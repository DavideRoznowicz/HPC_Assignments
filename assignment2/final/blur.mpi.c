/*
	File:     blur.mpi.c 
	Author:   Davide Roznowicz
	Version:  MPI
	Compiling (example):              openmpi/4.0.3/gnu/9.3.0
					  mpicc -O1 blur.mpi.c -o blur.mpi.x -lm
	Running (example):		  mpirun --mca btl ^openib -np ${procs} ./blur.mpi.x $kernel_type $kernel_size $ff $input_file $output_file 					  

		         
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <byteswap.h>
#include <math.h>
#include <sched.h>
#include <stdbool.h>
#include <mpi.h>


// options for kernel type
enum tipi{mean,weight,gaussian};

int binomialCoeff(int n,int k) { 
  int res=1; 
  if(k>n-k) 
    k=n-k; 
  for(int i=0;i<k;++i) { 
    res*=((n-i)/(i+1)); 
  } 
  return res; 
}


// accessory function for kernel
int* printPascal(int n) { 
  int* t=malloc(n*n*sizeof(int));
  n--;
  for (int i=0; i<= n;i++)
    t[i]=binomialCoeff(n,i);
  return t;
}

// accessory function for kernel
void print_kernel(float* matrix,int d) {
  for(int i=0;i<d;i++) {
    for(int j=0;j<d;j++)
      printf("%f ",matrix[i*d+j]);
    printf("\n");
  }
}




int main(int argc,char* argv[]) {



  
  int rank, tot_cores; // rank of each core/process ; total number of cores running
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD,&rank );
  MPI_Comm_size( MPI_COMM_WORLD,&tot_cores ); 



  

  FILE* fc; 			// declaring shared variable fc used at the end of it all
  float ff, w;  		// variables for weight kernel
  int* pascal;
  int type_kernel=weight; 	// options: mean, gaussian, weight 
  float acck=0; 		// accumulator for kernel matrix
  int d; 			// double offset
  char* input_file;  		// name of input file
  char* output_file;		// name of output file



  //+++++++++++START: parsing different command-line args according to assignment
  type_kernel=atoi(argv[1]); 		// avg, weight, gaussian -  0,1,2
  d=atoi(argv[2]);              	// dimension of kernel
  if (type_kernel==1) { 		// only in weight case
  	ff=atof(argv[3]);     	// ff only present for type_kernel==1
	input_file=malloc( (strlen(argv[4])+1) * sizeof(char) );
	input_file=argv[4];
	if (argc>5) {
		output_file=malloc( (strlen(argv[5])+1) * sizeof(char) );
		output_file=argv[5];
		printf("inp %s\n", (input_file) );
	}
  }
  else {  // argv[1]!=1
	input_file=malloc( (strlen(argv[3])+1) * sizeof(char) );
	input_file=argv[3];
	if (argc>4) {
		output_file=malloc( (strlen(argv[4])+1) * sizeof(char) );
		output_file=argv[4];
	}	
  }//++++++++++++END: parsing command-line args

	
	
	FILE* f;
	char c;
	char type[2];
	int stop=0;
	f=fopen(input_file,"r");   // opening file: "read mode"
	fscanf(f,"%s\n",type);  // parsing magic number
	
	if (rank==0) printf("type: %s\n", type);


	// +++++++++START PARSING COMMENTS
	while(!stop) { //entro nel loop
		fread(&c, sizeof(char),1, f);    // read one char
		fseek(f,-sizeof(char),SEEK_CUR); // go back one char
	while(c=='#') {   // if the char is #, read it all up to \n, otherwise exit
		while(c!='\n') { 
			fread(&c,sizeof(char),1,f);
		}
		fread(&c,sizeof(char),1,f);
		fseek(f,-sizeof(char),SEEK_CUR); // go back one char
	}
	stop=1;
	}
	// +++++++++END PARSING COMMENTS

	
	// parsing dimensions width and height: nx, ny
	int nx;
	int ny;
	fscanf(f,"%d %d\n",&nx,&ny); // scanning nx, ny

	// parsing maxval
	int maxval=0;
	fscanf(f,"%d\n",&maxval); //parsing maxval
	int offset=d/2;
        int nx_o=nx+d-1;
        int ny_o=ny+d-1;
        int size=d*d;

	
	if (rank==0){ // printing only for master
		printf("maxval: %d\n",maxval);
		printf("half size of kernel matrix: %d\n",offset);
		printf("nx %d ny %d\n",nx,ny);
	}
	

	
	// initial piece of matrix my_mwith offset (to be convoluted)
	unsigned short int* my_m;             // portion of matrix (with offset) a core is supposed to use
	unsigned short int* m_final;          // final piece of matrix that each core obtains after convolution
	unsigned short int* mylocal_buffer;   // storing piece of matrix just read by fread
	unsigned short int* garbage;          // store part of matrix got during reading but not used (before reaching the useful part)

	float* k=(float*)malloc(size*sizeof(float)); // kernel matrix ptr 
        float acc=0;                                 // accumulator for single point blurring
        unsigned short int tmp;
        int tot_pixels=nx*ny;
        int yfrom;
        int yto;
	
	// temporal vars
	int ktmp;
	int mtmp;
	int m_finaltmp;
	
	
	// time vars
	double tstart;
	double tend;
  if (rank==0) {
	tstart = MPI_Wtime(); // start time for MPI
  }

//****************************** Time starts
  
  // There is no reason to let one core read it all and then use scatter function. We can avoid a lot of
  // buffering and communication times just by letting each core open its own file and read its own piece.
  // It is almost sequential, but always much better than the other possibility just explained. Moreover,
  // while some cores are still busy reading, we let other cores take other workload to alleviate this
  // sequential behaviour.


  // How many lines should each core read? The optimal case would be that each core reads exactly 
  // the same number of lines. However, this is not always possible as ny might not be divisible by the number 
  // of cores. Then, we adopt the following choice: 
  // the first core gets ny/tot_cores (integer division) lines + one line. If there are still extra lines,
  // also core 1 gets the additional line; this continues until all extra lines have been assigned.
  // When this happens, each remaining core gets just ny/tot_cores (integer division).


  int extra_lines=ny % tot_cores;  // extra lines
  int base_lines=ny/tot_cores;
  int mylines=base_lines; // number of lines assigned to each core
  if (extra_lines!=0) {
  	if (rank<extra_lines) {
  		mylines++;  	
  	}
  }
  
  	
  

  

  // reading from file: each core reads its own piece of data in its own file: then initializes my_m;
  // the part of reading that is not useful for the computations of the specific cores is collected 
  // in garbage, then the pointer is freed.
  if (rank >= extra_lines) {
  	
  	my_m=(unsigned short int*)calloc(nx_o*(base_lines+d-1),sizeof(short int));
	m_final=(unsigned short int*)malloc(nx*(base_lines)*sizeof(short int));
	
	if (rank==0 && rank==tot_cores-1) { // edge case: last core is also the master ---> there's only one core
		yfrom=base_lines*rank+extra_lines;  // ==0
		yto=base_lines*rank+extra_lines+base_lines;  // ==base_lines

		fread(m_final, sizeof(short int), mylines*nx, f);
		for (int i=yfrom; i<yto; i++){ // i runs over rows
			for(int j=offset;j<nx_o-offset;j++) { // j runs over columns
				tmp=__bswap_16(m_final[(i-yfrom)*nx+(j-offset)]);
				my_m[(i+offset)*nx_o+j]=tmp;
			}
	  	}

	}

	else if (rank==tot_cores-1){  // last core
		yfrom=base_lines*rank+extra_lines-offset;
		yto=base_lines*rank+extra_lines+base_lines;
			
		// reading but ignoring piece of file that is not of interest for this core
		garbage=(unsigned short int*)malloc(yfrom*nx*sizeof(short int));
		fread(garbage, sizeof(short int), yfrom*nx, f);		


		// reading piece of interest
		mylocal_buffer=(unsigned short int*)malloc(nx*(mylines+offset)*sizeof(short int));
		fread(mylocal_buffer, sizeof(short int), (mylines+offset)*nx, f);
		for (int i=yfrom; i<yto; i++){
			for(int j=offset;j<nx_o-offset;j++) { // j runs over columns
				tmp=__bswap_16(mylocal_buffer[(i-yfrom)*nx+(j-offset)]);
				my_m[(i-yfrom)*nx_o+j]=tmp;
			}
		}
	
	
	}

	else if (rank==0 && rank!=tot_cores-1) { // first core, but there are also other cores this time
		yfrom=base_lines*rank+extra_lines;
                yto=base_lines*rank+extra_lines+base_lines+offset;
         
                // reading but ignoring piece of file that is not of interest for this core
                garbage=(unsigned short int*)malloc(yfrom*nx*sizeof(short int));
                fread(garbage, sizeof(short int), yfrom*nx, f);

         	
		// reading piece of interest
		mylocal_buffer=(unsigned short int*)malloc(nx*(mylines+offset)*sizeof(short int));
		fread(mylocal_buffer, sizeof(short int), (mylines+offset)*nx, f);
         	for (int i=yfrom; i<yto; i++){
         		for(int j=offset;j<nx_o-offset;j++) { // j runs over columns
         	                tmp=__bswap_16(mylocal_buffer[(i-yfrom)*nx+(j-offset)]);
         	                my_m[(i+offset)*nx_o+j]=tmp;
         	        }
         	}
         	
	}

	else {  // includes both perfect and non perfect divisibility
		yfrom=base_lines*rank+extra_lines-offset;
		yto=base_lines*rank+extra_lines+base_lines+offset;
		
		// reading but ignoring piece of file that is not of interest for this core
                garbage=(unsigned short int*)malloc(yfrom*nx*sizeof(short int));
                fread(garbage, sizeof(short int), yfrom*nx, f);
		


		// reading piece of interest
		mylocal_buffer=(unsigned short int*)malloc(nx*(mylines+d-1)*sizeof(short int));
		fread(mylocal_buffer, sizeof(short int), (mylines+d-1)*nx, f);
		for (int i=yfrom; i<yto; i++){
			for(int j=offset;j<nx_o-offset;j++) { // j runs over columns
				tmp=__bswap_16(mylocal_buffer[(i-yfrom)*nx+(j-offset)]);
				my_m[(i-yfrom)*nx_o+j]=tmp;
			}
		}
			
	}
  }
  else {    // rank < extra_lines
  	
  	my_m=(unsigned short int*)calloc(nx_o*(mylines+d-1),sizeof(short int));
  	m_final=(unsigned short int*)malloc(nx*(mylines)*sizeof(short int));
	
	
	if (rank==0) {
		yfrom=base_lines*rank+rank;
		yto=base_lines*rank+rank+mylines+offset;

		// reading but ignoring piece of file that is not of interest for this core
                garbage=(unsigned short int*)malloc(yfrom*nx*sizeof(short int));
                fread(garbage, sizeof(short int), yfrom*nx, f);

		
		// reading piece of interest
		mylocal_buffer=(unsigned short int*)malloc(nx*(mylines+offset)*sizeof(short int));
		fread(mylocal_buffer, sizeof(short int), (mylines+offset)*nx, f);
		for (int i=yfrom;i<yto;i++){ // i runs over rows
			for(int j=offset;j<nx_o-offset;j++) { // j runs over columns
				tmp=__bswap_16(mylocal_buffer[(i-yfrom)*nx+(j-offset)]);
				my_m[(i+offset)*nx_o+j]=tmp;
			}
		}
	
	
	
	
	}	
	else {
		yfrom=base_lines*rank+rank-offset;
		yto=base_lines*rank+rank+mylines+offset;
				
		// reading but ignoring piece of file that is not of interest for this core
                garbage=(unsigned short int*)malloc(yfrom*nx*sizeof(short int));
                fread(garbage, sizeof(short int), yfrom*nx, f);



		// reading piece of interest
		mylocal_buffer=(unsigned short int*)malloc(nx*(mylines+d-1)*sizeof(short int));
		fread(mylocal_buffer, sizeof(short int), (mylines+d-1)*nx, f);
		for (int i=yfrom;i<yto;i++){ // i runs over rows
			for(int j=offset;j<nx_o-offset;j++) { // j runs over columns
				tmp=__bswap_16(mylocal_buffer[(i-yfrom)*nx+(j-offset)]);
				my_m[(i-yfrom)*nx_o+j]=tmp;
			}
		}
		
  	}
  }


  if (tot_cores>1) {
	free(garbage); // free ptr garbage
	free(mylocal_buffer); // free ptr mylocal_buffer
  }

// master core declares and initializes the kernel matrix while other cores are still busy reading the file.
// As image matrixes tend to be large, while kernel matrixes are quite small, master core manages to finish 
// computing the kernel matrix before the reading is completed (reading from file is not parallel, but
// somehow sequential)


  if (rank==0) { // only master computes kernel
	// initializing kernel matrix
	switch(type_kernel) {
	
 	case mean:
		{
		float division_size=1.0/size;  // value of kernel matrix entries: division done only once
		for(int i=0;i<d;i++){
			for(int j=0; j<d; j++){
				k[i*d+j]=division_size;
			}
		}
		break;
 		}
 	case weight:
		{
		//ff=0.2;  // ff default
		w=(1-ff)/(size-1);
		for(int i=0;i<d;i++) // run over rows
			for (int j=0;j<d; j++) // run over columns
				k[i*d+j]=w;
		k[offset*d+offset]=ff;
		break;
 		}
	case gaussian:
		{
		pascal=malloc(size*sizeof(int));
		d--;
		for (int i=0; i<=d;i++)
			pascal[i]=binomialCoeff(d,i);
		for(int i=0;i<d;i++)
			for(int j=0;j<d;j++) {
				k[i*d+j]=pascal[i]*pascal[j];
				acck+=k[i*d+j];
			}
		float division_acck=1.0/acck;		 
		for(int i=0;i<d;i++)
			for(int j=0;j<d;j++)
				k[i*d+j]=k[i*d+j]*division_acck;
		break;
		}
	}
  
  } // end of kernel computation by master
  
  
  // sending kernel matrix to all cores: MPI_Bcast(buf, dim_of_buf, type, root, communicator)
  MPI_Bcast(k, size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  


  
  // ++++++++++++START: each core performs convolution on its own piece of matrix and saves the results into m_final
  for(int i=offset;i<offset+mylines;i++) {
	m_finaltmp=(i-offset)*nx-offset;
	for(int j=offset;j<nx_o-offset;j++) {
		acc=0;
		for(int ii=-offset;ii<=offset;ii++) {
			ktmp=(ii+offset)*d+offset; 
			mtmp=(i+ii)*nx_o+j;        
			for(int jj=-offset;jj<=offset;jj++) {
				acc+=k[ktmp+jj]*my_m[mtmp+jj];
			}
		}
		m_final[m_finaltmp+j]=__bswap_16((unsigned short int) acc);	
	}
  }

  free(my_m); // free ptr my_m  
  free(k);    // free ptr k  

  // ++++++++++++END


  // ++++++++++++++++START: Preparing input for Gatherv

  unsigned short int* m_write;
  // Allocating receive buffer to master core
  if (rank==0) {
  	m_write=(unsigned short int*)malloc(nx*ny*sizeof(short int));
  }


  // array containing the displacement in sendbuf
  int* displs = (int *)malloc(tot_cores*sizeof(int)); 
  
  // array cointaning the number of elements sent from each process
  int* rcounts = (int *)malloc(tot_cores*sizeof(int));
  
  if (rank==0) {
  	for (int i=0; i<tot_cores; i++) {
  		if (i<extra_lines) {
  			rcounts[i]=nx*(base_lines+1);
			displs[i]=nx*((base_lines+1)*i);
        	}
		else if (i==extra_lines) {
			rcounts[i]=nx*base_lines;
			displs[i]=nx*((base_lines+1)*extra_lines);

		}
        	else {   // i>extra_lines
     			rcounts[i]=nx*base_lines;
			displs[i]=nx*((base_lines+1)*extra_lines)+nx*(base_lines*(i-extra_lines));
     		}
  	}
  }
  
  int mypixels=nx*mylines;

  // ++++++++++++++++END: Preparing input for Gatherv


  // Sending m_final: each core sends its piece of blurred image to master
  MPI_Gatherv(m_final, mypixels, MPI_UNSIGNED_SHORT,
		m_write, rcounts, displs, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);	          
  MPI_Barrier(MPI_COMM_WORLD);
  
  free(m_final); // free ptr m_final



  // ++++++++++++START: producing the output name for the image (according to the assignment requests)
  int len_mystring=strlen(input_file)+ 23 + (int) (floor(log10(abs(d))) + 1) +(int) (floor(log10(abs(d))) + 1); 
  char mystring[len_mystring];
  if (rank==0){
	
	char kern_str[(int) (floor(log10(abs(d))) + 1) +1];
	char type_string[2];  // one more position for end of string
	sprintf(kern_str, "%d", d);
	sprintf(type_string, "%d", type_kernel );
	input_file[strlen(input_file)-4]='\0'; // wanna get rid of .pgm
	if (type_kernel==1) {
		snprintf(mystring, sizeof(mystring), "%s.b_%s_%sx%s_02.mpi.pgm", input_file, type_string, kern_str, kern_str);
		if (argc<=5){ // no output_file given in input
			output_file=malloc( (len_mystring)*sizeof(char) );
			strcpy(output_file, mystring);
		}
		printf("output default name for this specific image: %s", (mystring) ); 
	}
	else { // different kernels...ff not present
                snprintf(mystring, sizeof(mystring), "%s.b_%s_%sx%s.mpi.pgm", input_file, type_string, kern_str, kern_str);
                if (argc<=4){ // no output_file given in input
			output_file=malloc( (len_mystring)*sizeof(char) );
                        strcpy(output_file, mystring);
                }
                printf("output default name for this specific image: %s", (mystring) );			
	}
  }//+++++++++++++END: producing the output name
  

  //++++++++START WRITING master core opens the file and writes the header; then writes the whole final matrix
  if (rank==0) {

	fc=fopen(output_file, "w");
	fprintf(fc,"P5\n%d %d\n%d\n", nx, ny,maxval);
	fclose(fc);
	fc=fopen(output_file,"ab"); // open file for output at the end of file ( P5, nx, ny, maxval added before) : a=append, b=binary mode


  	// writing the image in the file
	fwrite(m_write,sizeof(short int),nx*ny,fc);
	fclose(fc);
	
  } //++++++++END WRITING 

  if (rank==0) {
  	tend = MPI_Wtime();          // end time for MPI
  	printf("\ntime for %d cores : %g sec\n\n", tot_cores, tend-tstart);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  



return 0;  


} // end of main
















