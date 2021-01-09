#include <stdio.h>
#include <stdlib.h>
#include <byteswap.h>
#include <math.h>
enum tipi{mean,weigth,gaussian};

int binomialCoeff(int n,int k) { 
  int res=1; 
  if(k>n-k) 
    k=n-k; 
  for(int i=0;i<k;++i) { 
    res*=(n-i); 
      res/=(i+1); 
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
// +++++++ now unsigned short int instead of short int
void write_img(unsigned short int* matrix,char* name,int nx,int ny,int lev,int d) {
  FILE* fc;
  fc=fopen(name,"w");
  fprintf(fc,"P5\n%d %d\n%d\n", nx, ny,lev);
  fclose(fc);
  fc=fopen(name,"ab"); // open file for output at the end of file ( P5, nx, ny, maxval added before) : a=append, b=binary mode
  int offset=d/2;
  int nx_o=nx+d-1;  // let's go up to d: result is line made of      offset+nx+offset
  int ny_o=ny+d-1;
  for(int i=offset;i<ny_o-offset;i++)
    for(int j=offset;j<nx_o-offset;j++) {
      unsigned short int tmp=matrix[i*nx_o+j];
      tmp=__bswap_16(tmp);
      fwrite(&tmp,sizeof(short int),1,fc);
    }
  fclose(fc);
}


double* kernel(int type,int d) {
  double* k=malloc(d*d*sizeof(double));
  int size=d*d;
  double ff=0;
  double w=0;
  int* pascal;
  double acc=0;
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
	acc+=k[i*d+j];
      }
    for(int i=0;i<d;i++) 
      for(int j=0;j<d;j++) 
	k[i*d+j]=k[i*d+j]/(float)acc;
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
  double acc=0;
  for(int i=0;i<d;i++)
    for(int j=0;j<d;j++)
      acc+=matrix[i*d+j];
  printf("%f\n",acc);
}

int main(int argc,char* argv[]) {
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
  int lev=0;
  fscanf(f,"%d\n",&lev);
  printf("lev max %d\n",lev); // lev is maxval
  printf("SPESSORE BORDO ESTERNO %d\n",d/2);
  int offset=d/2;
  int nx_o=nx+d-1;
  int ny_o=ny+d-1;
  int size=d*d;
  
  //alloco matrice con spazio bordi e inizializzo a 0
  unsigned short int* m=malloc(nx_o*ny_o*sizeof(short int));
  unsigned short int* m1=malloc(nx_o*ny_o*sizeof(short int));
  for(int i=0;i<nx_o*ny_o;i++)
    m[i]=0;
  int tmp;
  
  //leggo immagine
  for(int i=offset;i<ny_o-offset;i++)
    for(int j=offset;j<nx_o-offset;j++) {
      fread(&tmp, sizeof(short int), 1, f);
      tmp=__bswap_16(tmp); 
      m[i*nx_o+j]=tmp;
    }

  double* k;
#ifdef DEBUG
  printf("PUTTANATE DI DEBUG\n");
  k=kernel(mean, d);
  print_kernel(k, d);
  print_sum(k, d);
  printf("\n");
  k=kernel(gaussian,d);
  print_kernel(k, d);
  print_sum(k, d);
  printf("\n");
  k=kernel(weigth,d);
  print_kernel(k, d);
  print_sum(k, d);
#endif

  long int diff=0;
  k=kernel(weigth,d);
  print_kernel(k, d);
  for(int i=offset;i<ny_o-offset;i++) {
    for(int j=offset;j<nx_o-offset;j++) {
      double acc=0;
      for(int ii=-d/2;ii<=d/2;ii++) 
	for(int jj=-d/2;jj<=d/2;jj++)
	  acc+=k[(ii+d/2)*d+(jj+d/2)]*m[(i+ii)*nx_o+(j+jj)];
      m1[i*nx_o+j]=(unsigned short int)acc;
    }
  }
  


  write_img(m1,"copia.pgm", nx, ny, lev, d); // nx, ny, maxval, d==twice the offset     
  
  return 0;  
}
