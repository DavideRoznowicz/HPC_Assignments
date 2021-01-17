/*
	File:     blur.omp.c 
	Author:   Davide Roznowicz
	Version:  OpenMP
	Compiling (example):          // load gnu>4.9.9 to use proc_bind clause
			 		module load gnu/9.3.0
		         		export OMP_NUM_THREADS=24; gcc -O3 -std=gnu99 -fopenmp blur.omp.c -o blur.omp.x -lm;
	Running (example):		./blur.omp.x $kernel_type $kernel_size $ff $input_file $output_file

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <byteswap.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <sched.h>






#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif
#define _GNU_SOURCE



#if defined(_OPENMP)
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)

#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec +	\
		     (double)myts.tv_nsec * 1e-9)

#else

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)

#define CPU_TIME_th (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)
#endif




// kind of kernels available
enum tipi{mean,weight,gaussian};


// acessory function to gaussian kernel
int binomialCoeff(int n,int k) { 
  int res=1; 
  if(k>n-k) 
    k=n-k; 
  for(int i=0;i<k;++i) { 
    res*=((n-i)/(i+1)); 
  } 
  return res; 
}

// accessory function to gaussian kernel
int* printPascal(int n) { 
  int* t=malloc(n*n*sizeof(int));
  n--;
  for (int i=0; i<= n;i++)
    t[i]=binomialCoeff(n,i);
  return t;
}



// print kernel matrix
void print_kernel(double* matrix,int d) {
  for(int i=0;i<d;i++) {
    for(int j=0;j<d;j++)
      printf("%f ",matrix[i*d+j]);
    printf("\n");
  }
}








int main(int argc,char* argv[]) {


  // first parsing is made by master thread
  float ff, w;                  // variables for weight kernel
  int type_kernel=weight;
  int d;			// kernel dimension
  char* input_file;             // name of input file
  char* output_file;            // name of output file



  //+++++++++++START: parsing different command-line args according to assignment
  type_kernel=atoi(argv[1]);            // avg, weight, gaussian -  0,1,2
  d=atoi(argv[2]);                      // dimension of kernel
  if (type_kernel==1) {                 // only in weight case
  	ff=atof(argv[3]);       	// ff only present for type_kernel==1
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
  




  FILE* f; // master thread opens file to start parsing
  char c;
  char type[2];
  int stop=0;
  f=fopen(input_file,"r");   // open file: "read mode"
  fscanf(f,"%s\n",type);  // parsing magic number
  printf("type: %s\n",type);

  // +++++++++START PARSING COMMENTS
  while(!stop) {
    fread(&c, sizeof(char),1, f); // read a char
    fseek(f,-sizeof(char),SEEK_CUR); // go back one char
    while(c=='#') { // if char is  #, reading it all up to \n, otherwise exiting
      while(c!='\n') { 
	fread(&c,sizeof(char),1,f);
      }
      fread(&c,sizeof(char),1,f);
      fseek(f,-sizeof(char),SEEK_CUR); // go back one char
    }
    stop=1;
  }
  // +++++++++END PARSING COMMENTS


  // parsing dimensions: nx, ny
  int nx; // width of the image
  int ny; // height of the image
  fscanf(f,"%d %d\n",&nx,&ny); // parsing of dimensions
  printf("nx %d ny %d\n",nx,ny); 


  //parsing maxval
  int maxval=0;
  fscanf(f,"%d\n",&maxval); // parsing maxval
  printf("maxval %d\n",maxval);
  printf("half-size %d\n",d/2);
  int offset=d/2;
  int nx_o=nx+d-1;
  int ny_o=ny+d-1;
  int size=d*d;



 
  
  // allocating matrix encircled by zeros 0 (considered an offset), but not initialized yet

  // malloc because we are paying attention to "first touch" policy
  unsigned short int* m=malloc(nx_o*ny_o*sizeof(short int)); // initial matrix m with offset (to be convoluted)
  unsigned short int* m1=malloc(nx*ny*sizeof(short int)); // final matrix m1 (obtained after convolution)
  float* k=malloc(d*d*sizeof(float)); // kernel matrix ptr
  unsigned short int* tmp_array; // temporary array for speeding up reading
  float acc=0; // accumulator for single point blurring 
  FILE* fc; // declaring shared variable fc: it will be used for writing in the end

  float acck=0; // accumulator for kernel matrix
  float division_size;
  float division_acck;
  int* pascal; 
  int tot_threads; // total number of spawned threads

  //various useful tmp vars
  int ktmp;
  int mtmp;
  int m1tmp;
  int ptmp;


  struct timespec ts; // for time storage
  double tend;
  double tstart = CPU_TIME; // starting time for parallel region


#pragma omp parallel shared(acc) proc_bind(close) // if we are among cores, close allows us to reduce time to retrieve data from memory
{	
	


//*************START   initialize zeros only on the "contour" of the matrix; enabling "first touch" for data in m


	// need to initialize to zero the first offset lines by first thread (master)
	#pragma omp master
	{
		for (int i=0; i<offset; i++){
			mtmp=i*nx_o;
			for (int j=0; j<nx_o; j++){
	        		m[mtmp+j]=0;
	        	}
		}
	}

	#pragma omp for schedule(static) nowait  // with static, equal-size chuncks are assigned to threads (n-th chunck to n-th thread)
        for (int i=offset;i<ny_o-offset;i++){ // i runs over rows
		mtmp=i*nx_o;
                for (int j=0; j<offset; j++){
                        m[mtmp+j]=0;
                }

                for (int j=nx_o-offset; j<nx_o; j++){
                        m[mtmp+j]=0;
                }
        } 

	if (omp_get_thread_num()==omp_get_num_threads()-1){ // need to initialize to zero the first offset lines by last thread
                for (int i=ny_o-offset; i<ny_o; i++){
			mtmp=i*nx_o;
                        for (int j=0; j<nx_o; j++){
                                m[mtmp+j]=0;
                        }
                }
        }




        #pragma omp single
        {
        tot_threads=omp_get_num_threads();
	tmp_array=(unsigned short int*)malloc(nx*sizeof(short int)); // allocating tmp_array
        }
        // reading from file: each thread reads its own piece of data and initializes a piece of m      
        unsigned short int tmp;
        #pragma omp for ordered schedule(static) private(tmp) nowait
        for(int i=offset;i<ny_o-offset;i++){ // i runs over rows
        	#pragma omp ordered  // forcing ordering as fread requires to be able to read data in order
        	{
		fread(tmp_array, sizeof(short int), nx, f);
		mtmp=i*nx_o;	
		for(int j=offset;j<nx_o-offset;j++) { // j runs over columns
        		tmp=__bswap_16(tmp_array[j-offset]);
        		m[mtmp+j]=tmp; // allocating to matrix m, enabling first touch policy
		}
                }
        }
        


//*************END   initialize zero matrix, enabling "first touch" for data in m





	// initializing kernel matrix
	switch(type_kernel) {
	
 	case mean:
		{
		#pragma omp single
		{
			division_size=1.0/size;
		}
		#pragma omp for schedule(dynamic)  // we don't care about the order here because the kernel matrix is used by any thread
		for(int i=0;i<d;i++){
			for(int j=0; j<d; j++){
				k[i*d+j]=division_size;
			}
		}
		break;
 		}
 	case weight:
		{
		#pragma omp single
		{
			w=(1-ff)/(size-1);
		}
		#pragma omp for schedule(dynamic)  // we don't care about the order here because the kernel matrix is used by any thread
		for(int i=0;i<d;i++) // run over rows
			for (int j=0;j<d; j++) // run over columns
				k[i*d+j]=w;
		#pragma omp single
		{
			k[offset*d+offset]=ff;
		}
		break;
 		}
	case gaussian:
		{
		#pragma omp single
		{
                	pascal=malloc(size*sizeof(int));
                	d--;
		}
		#pragma omp for schedule(dynamic)
		for (int i=0; i<=d;i++)
			pascal[i]=binomialCoeff(d,i);
		#pragma omp for reduction(+: acck) schedule(dynamic) 
		// we don't care about the order here because the kernel matrix is used by any thread
		for(int i=0;i<d;i++){
			for(int j=0;j<d;j++) {
				k[i*d+j]=ptmp*pascal[j];
				acck+=k[i*d+j];
			}
		}
		#pragma omp single
		{
			division_acck=1.0/acck;
		}
		#pragma omp for schedule(dynamic) // we don't care about the order here because the kernel matrix is used by any thread 
		for(int i=0;i<d;i++)
			for(int j=0;j<d;j++)
				k[i*d+j]=k[i*d+j]*division_acck;
		break;
		}
	}




	// ++++++++++++++++START OF COMPUTATIONS for blurring

	#pragma omp for reduction(+: acc) private(ktmp, mtmp, m1tmp) firstprivate(k)  // reduction in order to avoid data race on acc
	for(int i=offset;i<ny_o-offset;i++) {
		m1tmp=(i-offset)*nx-offset;
		for(int j=offset;j<nx_o-offset;j++) {
			acc=0;
			for(int ii=-offset;ii<=offset;ii++) {
				ktmp=(ii+offset)*d+offset;
				mtmp=(i+ii)*nx_o+j;
				for(int jj=-offset;jj<=offset;jj++) {
					acc+=k[ktmp+jj]*m[mtmp+jj];
				}
			}
			m1[m1tmp+j]=__bswap_16( (unsigned short int) acc );
			
		}
	}

	// ++++++++++++++++END OF COMPUTATIONS for blurring
	

 
	#pragma omp single // one thread that finished previous tasks earlier opens the file and writes the header
	{
	
	  free(m); // free ptr m
          free(k); // free ptr k


	  // ++++++++++++START: producing the output name for the image (according to the assignment requests)
	  int len_mystring=strlen(input_file)+ 23 + (int) (floor(log10(abs(d))) + 1) +(int) (floor(log10(abs(d))) + 1);
	  char mystring[len_mystring];

		   char kern_str[(int) (floor(log10(abs(d))) + 1) +1];
		   char type_string[2];  // one more position for end of string
		   sprintf(kern_str, "%d", d);
		   sprintf(type_string, "%d", type_kernel );
		   input_file[strlen(input_file)-4]='\0'; // wanna get rid of .pgm
		   if (type_kernel==1) {
		           snprintf(mystring, sizeof(mystring), "%s.b_%s_%sx%s_02.omp.pgm", input_file, type_string, kern_str, kern_str);
		           if (argc<=5){ // no output_file given in input
		                   output_file=malloc( (len_mystring)*sizeof(char) );
		                   strcpy(output_file, mystring);
		           }
		           printf("output default name for this specific image: %s", (mystring) );
		   }
		   else { // different kernels...ff not present
		           snprintf(mystring, sizeof(mystring), "%s.b_%s_%sx%s.omp.pgm", input_file, type_string, kern_str, kern_str);
		           if (argc<=4){ // no output_file given in input
		                   output_file=malloc( (len_mystring)*sizeof(char) );
		                   strcpy(output_file, mystring);
		           }
		           printf("output default name for this specific image: %s", (mystring) );
		   }
	  //+++++++++++++END: producing the output name









	// ++++++++++START OF WRITING 

	fc=fopen(output_file,"w");  // image name to write
	fprintf(fc,"P5\n%d %d\n%d\n", nx, ny,maxval);  // writing magic number, dimensions and maxval
	fclose(fc);  // closing file
	fc=fopen(output_file,"ab"); // reopen file for output, but in binary mode: a=append, b=binary mode
	}

	#pragma omp master
	{

	fwrite(m1,sizeof(short int),nx*ny,fc); // writing the whole matrix in just one shot by master thread
	fclose(fc);  

	}

	// ++++++++++END OF WRITING	


} // end of parallel region

  tend = CPU_TIME; // end time for parallel region, for master thread
  printf("\ntime for %d threads : %g sec\n\n", tot_threads, tend-tstart);
  
return 0;  


} // end of main
















