/*
	File:     myblur.c 
	Author:   Davide Roznowicz
	Version:  OpenMP
	Compilation (example):          # load gnu>4.9.9 to use proc_bind clause
			 		module load gnu/9.3.0
		         		export OMP_NUM_THREADS=24; gcc -std=gnu99 -fopenmp myblur.c -o myblur.x; ./myblur.x check_me.pgm 31

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



#define N_default 1000


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

#ifdef OUTPUT
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...)
#endif

#define CPU_ID_ENTRY_IN_PROCSTAT 39
#define HOSTNAME_MAX_LENGTH      200



int read_proc__self_stat ( int, int * );
int get_cpu_id           ( void       );






enum tipi{mean,weigth,gaussian};

int binomialCoeff(int n,int k) { 
  int res=1; 
  if(k>n-k) 
    k=n-k; 
  for(int i=0;i<k;++i) { 
    res*=((n-i)/(i+1)); 
  } 
  return res; 
}

int* printPascal(int n) { 
  int* t=malloc(n*n*sizeof(int));
  n--;
  for (int i=0; i<= n;i++)
    t[i]=binomialCoeff(n,i);
  return t;
}


double* kernel(int type,int d) {
  double* k=malloc(d*d*sizeof(double));
  int size=d*d;
  double ff=0;
  double w=0;
  int* pascal;
  double acck=0;
  switch(type) {
  case mean:
    for(int i=0;i<size;i++)
      k[i]=1.0/size;
    break;

  case weigth:
    ff=0.2;  // ff default
    w=(1-ff)/(d*d-1);
    for(int i=0;i<size;i++)
      k[i]=w;
    k[(d/2)*d+(d/2)]=ff;
    break;

  case gaussian:
    pascal=printPascal(d);
    for(int i=0;i<d;i++) 
      for(int j=0;j<d;j++) { 
	k[i*d+j]=pascal[i]*pascal[j];
	acck+=k[i*d+j];
      }
    for(int i=0;i<d;i++) 
      for(int j=0;j<d;j++) 
	k[i*d+j]=k[i*d+j]/(float)acck;
    break;
  }
  return k;
}

void print_kernel(double* matrix,int d) {
  for(int i=0;i<d;i++) {
    for(int j=0;j<d;j++)
      printf("%f ",matrix[i*d+j]);
    printf("\n");
  }
}


// print sum of kernel matrix entries
void print_sum(double* matrix,int d) {
  double acck=0;
  for(int i=0;i<d;i++)
    for(int j=0;j<d;j++)
      acck+=matrix[i*d+j];
  printf("%f\n",acck);
}









int main(int argc,char* argv[]) {


  // first parsing is made by master thread


  //sistemo gli input
  if(argc<2) {
    printf("ERRORE ERRORE MANCA IL NOME DEL FILE\n");
    return -1;
  }
  int d=3;
  if(argc>2) {
    d=atoi(argv[2]);
    if( d%2==0) {
      printf("NON VALIDO; DEVE ESSERE UN NUMERO DISPARI; UTILIZZO QUELLO DI DEF=3");
      d=3;
    }
    printf("DIM KERNEL D=%d\n",d);
  }
  else {
    printf("DIM KERNEL DEFAULT=3\n");
  }
  
  FILE* f;
  char c;
  char type[2];
  int stop=0;
  f=fopen(argv[1],"r");//parsing tipo immagine
  fscanf(f,"%s\n",type);  //parsing magic number
  printf("tipo: %s\n",type);

  //parsing commenti
  while(!stop) { //entro nel loop
    fread(&c, sizeof(char),1, f); //leggo un carattere
    fseek(f,-sizeof(char),SEEK_CUR); //torno indietro di uno
    while(c=='#') { //se il carattere è # leggo tutto fino a \n,altrimenti esco
      while(c!='\n') { 
	fread(&c,sizeof(char),1,f);
      }
      fread(&c,sizeof(char),1,f);
      fseek(f,-sizeof(char),SEEK_CUR); //torno indietro di uno
    }
    stop=1;
  }



  //parsing dimensioni: nx, ny
  int nx;
  int ny;
  fscanf(f,"%d %d\n",&nx,&ny);
  printf("nx %d ny %d\n",nx,ny);


  //parsing livello massimo: maxval
  int maxval=0;
  fscanf(f,"%d\n",&maxval);
  printf("maxval max %d\n",maxval);
  printf("SPESSORE BORDO ESTERNO %d\n",d/2);
  int offset=d/2;
  int nx_o=nx+d-1;
  int ny_o=ny+d-1;
  int size=d*d;



 
  
  //alloco matrice con spazio bordi e inizializzo a 0

  // malloc because we are paying attention "first touch" policy
  unsigned short int* m=malloc(nx_o*ny_o*sizeof(short int)); // initial matrix m wit offset (to be convoluted)
  unsigned short int* m1=malloc(nx_o*ny_o*sizeof(short int)); // final matrix m1 (obtained after convolution)
  float* k=malloc(d*d*sizeof(float)); // kernel matrix ptr 
  float acc=0; // accumulator for single point blurring 
  FILE* fc; // declaring shared variable fc 

  float ff, w;  // variables for weight kernel
  float acck=0; // accumulator for kernel matrix
  float myk;
  int* pascal;
  int type_kernel=mean; // options: mean, gaussian, weight 
  int tot_threads;



  struct  timespec ts; // for time storage
  double tend;
  double tstart = CPU_TIME; // start time for parallel region


#pragma omp parallel shared(acc) proc_bind(close) // if we are among cores, close allows us to reduce time to retrieve data from memory
{	
	
	#pragma omp single nowait
        {
        tot_threads=omp_get_num_threads();
        }
        // reading from file: each thread reads its own piece of data and initializes a piece of m      
                 int tmp;
                 #pragma omp for ordered schedule(static) private(tmp) nowait
                 for(int i=offset;i<ny_o-offset;i++){ // i runs over rows
                 	#pragma omp ordered  // forcing ordering as fread requires oredered data reading from opened file f
                        for(int j=offset;j<nx_o-offset;j++) { // j runs over columns
                        	fread(&tmp, sizeof(short int), 1, f);
                                tmp=__bswap_16(tmp);
                                m[i*nx_o+j]=tmp;
                	}
                 }
 




	/*
	#pragma omp single nowait   // on purpose I decide not to parallelize: this is because the previous for is oredered but nowait: the master should \
					finish first and then compute the kernel matrix, while most of the others are still reading: as the number of \
					cores increases, scalability should be enhanced.

	switch(type_kernel) {
         
        case mean:
           
          
                        myk=1.0/size;  // value of kernel matrix entries: division done only one
                 // we don't care about the order here because the kernel matrix is used by any thread
                for(int i=0;i<d;i++){
                        for(int j=0; j<d; j++){
                                k[i*d+j]=myk;
                        }
                }
                
		break;
  		
        case weigth:
                #pragma omp single nowait
                {
                        ff=0.2;  // ff default
                        w=(1-ff)/(size-1);
                
                // we don't care about the order here because the kernel matrix is used by any thread
                for(int i=0;i<d;i++) // run over rows
                        for (int j=0;j<d; j++) // run over columns
                                k[i*d+j]=w;
                
                
                        k[offset*d+offset]=ff;
                }
                break;
  		
        case gaussian:
                #pragma omp single nowait
                {
                        pascal=malloc(size*sizeof(int));
                        d--;
                
                
                for (int i=0; i<=d;i++)
                        pascal[i]=binomialCoeff(d,i);
                 // we don't care about the order here because the kernel matrix is used by any thread
                for(int i=0;i<d;i++)
                        for(int j=0;j<d;j++) {
                                k[i*d+j]=pascal[i]*pascal[j];

			}
		}
		break;
	}	

	*/




       

//*************START   initialize zero matrix, enabling "first touch" for data in m


	// need initialize first offset lines of zeros by first thread (master)
	#pragma omp master
	{
		for (int i=0; i<offset; i++){
			for (int j=0; j<nx_o; j++){
	        		m[i*nx_o+j]=0;
	        	}
		}
	}

	#pragma omp for schedule(static) nowait  // with static, equal-size chuncks are assigned to threads (n-th chunck to n-th thread)
        for (int i=offset;i<ny_o-offset;i++){ // i runs over rows

                for (int j=0; j<offset; j++){
                        m[i*nx_o+j]=0;
                }

                for (int j=nx_o-offset; j<nx_o; j++){
                        m[i*nx_o+j]=0;
                }
        }

	if (omp_get_thread_num()==omp_get_num_threads()-1){ // only last thread should initialize the last zero(offset) lines to enforce first touch
                for (int i=ny_o-offset; i<ny_o; i++){
                        for (int j=0; j<nx_o; j++){
                                m[i*nx_o+j]=0;
                        }
                }
        }
	/*
	#pragma omp single nowait
	{
	tot_threads=omp_get_num_threads();
	}
	// reading from file: each thread reads its own piece of data and initializes a piece of m 	
	//int tmp;
	#pragma omp for ordered schedule(static) private(tmp) //nowait
	for(int i=offset;i<ny_o-offset;i++){ // i runs over rows
		#pragma omp ordered  // forcing ordering as fread requires oredered data reading from opened file f
                for(int j=offset;j<nx_o-offset;j++) { // j runs over columns
			fread(&tmp, sizeof(short int), 1, f);
                	tmp=__bswap_16(tmp);
                	m[i*nx_o+j]=tmp;

		}

	}
*/

//*************END   initialize zero matrix, enabling "first touch" for data in m


	/*

	// initializing kernel matrix
	switch(type_kernel) {
	
 	case mean:
		#pragma omp single
		{
			myk=1.0/size;  // value of kernel matrix entries: division done only one
		}
		#pragma omp for schedule(dynamic) // we don't care about the order here because the kernel matrix is used by any thread
		for(int i=0;i<d;i++){
			for(int j=0; j<d; j++){
				k[i*d+j]=myk;
			}
		}
		break;
 
 	case weigth:
		#pragma omp single
		{
			ff=0.2;  // ff default
			w=(1-ff)/(size-1);
		}
		#pragma omp for schedule(dynamic) // we don't care about the order here because the kernel matrix is used by any thread
		for(int i=0;i<d;i++) // run over rows
			for (int j=0;j<d; j++) // run over columns
				k[i*d+j]=w;
		#pragma omp single
		{
			k[offset*d+offset]=ff;
		}
		break;
 
	case gaussian:
		#pragma omp single
		{
			pascal=malloc(size*sizeof(int));
			d--;
		}
		#pragma omp for schedule(dynamic)
		for (int i=0; i<=d;i++)
			pascal[i]=binomialCoeff(d,i);
		#pragma omp for reduction(+: acck) schedule(dynamic) // we don't care about the order here because the kernel matrix is used by any thread
		for(int i=0;i<d;i++)
			for(int j=0;j<d;j++) {
				k[i*d+j]=pascal[i]*pascal[j];
				acck+=k[i*d+j];
			}
		#pragma omp for schedule(dynamic) // we don't care about the order here because the kernel matrix is used by any thread 
		


*/

 
/* 
 
  
#ifdef DEBUG
  printf("DEBUG\n");
  k=kernel(mean, d);
  print_kernel(k, d);
  print_sum(k, d);
  printf("\n");
  k=kernel(gaussian,d);
  print_kiernel(k, d);
  print_sum(k, d);
  printf("\n");
  k=kernel(weigth,d);
  print_kernel(k, d);
  print_sum(k, d);
#endif

*/
	/*
	#pragma omp single
	{
		print_kernel(k,d);
	}	
	*/
	#pragma omp for reduction(+: acc) nowait firstprivate(k)// avoiding data race on acc
	for(int i=offset;i<ny_o-offset;i++) {
		for(int j=offset;j<nx_o-offset;j++) {
			acc=0;
			for(int ii=-offset;ii<=offset;ii++) {
				for(int jj=-offset;jj<=offset;jj++) {
					acc+=k[(ii+offset)*d+(jj+offset)]*m[(i+ii)*nx_o+(j+jj)];
				}
			}
			m1[i*nx_o+j]=(unsigned short int) acc;
			
		}
	}


	#pragma omp single // one thread that finished earlier previous tasks opens the file and writes the header
	{
	fc=fopen("omp_after.pgm","w");
	fprintf(fc,"P5\n%d %d\n%d\n", nx, ny,maxval);
	fclose(fc);
	fc=fopen("omp_after.pgm","ab"); // open file for output at the end of file ( P5, nx, ny, maxval added before) : a=append, b=binary mode
	}

	#pragma omp for ordered schedule(static)  private(tmp)  // all threads write in order each piece of the blurred image 
	for(int i=offset;i<ny_o-offset;i++) {
		#pragma omp ordered
		for(int j=offset;j<nx_o-offset;j++) {
			unsigned short int tmp=m1[i*nx_o+j];
			tmp=__bswap_16(tmp);
			fwrite(&tmp,sizeof(short int),1,fc);
		}
	} // implicit barrier: the file must be closed only after all threads have finished writing
	#pragma omp single nowait // one thread closes the file (might be different from the one that opened it)
	{
		fclose(fc);
	}

} // end of parallel

  tend = CPU_TIME; // end time for parallel region
  printf("\ntime for %d threads : %g sec\n\n", tot_threads, tend-tstart);
  
return 0;  


} // end of main
















