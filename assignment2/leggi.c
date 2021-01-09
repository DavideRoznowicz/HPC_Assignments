#include <stdio.h>
#include <stdlib.h>
#include <byteswap.h>
int main(int argc,char* argv[]) {
  
  FILE* f;
  f=fopen(argv[1],"r");
  char c;
  char type[2];
  int stop=0;
  //read magic number
  fread(&c, sizeof(char),1, f);
  printf("%c",c); // P
  while(c!='\n') {
    fread(&c, sizeof(char),1, f);
    printf("%c",c); // 5
  }
  //parsing commenti una sola riga
  fread(&c,sizeof(char),1, f);
  //printf("%c",c); // # or if comments missing, then first char of width
  while(c=='#') { //se è un # inizia la linea di commento
    while(c!='\n') { //se è != da \n continuo a leggere
      fread(&c,sizeof(char),1, f);
      printf("%c",c); // text of comment or #
    }
    fread(&c,sizeof(char),1, f); //leggo un altro carattere e poi vedo se == #
    //printf("%c",c); // first character of width (line without #) 
  }
  //parsing dimensioni
  // width(nx) - height(ny)
  char dim[100];
  int count=1;
  dim[0]=c; // see some lines ago
  fread(&c,sizeof(char),1, f);
  dim[count]=c;
  //printf("%c",c);
  int base=0;
  while((c!='\n') && (c!=' ')) {
    count++;
    fread(&c,sizeof(char),1, f); // first elem was already taken before, but we also have final \n. 
    dim[count]=c;
    //printf("%c",c);
  }
  int len_xsize = count-base;  // len of xsize
  base=count;
  
  count++;
  fread(&c,sizeof(char),1, f);
  dim[count]=c;
  // height(ny)
  while((c!='\n') && (c!=' ')) {
	count++;	
	fread(&c,sizeof(char),1, f); // first elem was already taken before, but we also have final \n. 
	dim[count]=c;
	//printf("%c",c);
  }
  int len_ysize = count-base-1;  // len of ysize
  base=count;

  

  
  // prints gray max level, keep printing up to \n
  fread(&c,sizeof(char),1, f);
  count++;
  dim[count]=c;
  //printf("%c",c);  
  while(c!='\n') {
    count++;
    fread(&c,sizeof(char),1, f);
    dim[count]=c;
    //printf("%c",c);
  }

  int len_maxval = count-base-1;  // len of maxval
  //printf("lenx=%d, leny=%d\nlen_maxval=%d\n", len_xsize, len_ysize, len_maxval);

  
  int kk=0;
  char* mychar=malloc(len_xsize*sizeof(char));
  
  while (dim[kk]!=' ' && dim[kk]!='\n') 
  {
	mychar[kk]=dim[kk];
	kk++;
  }
  kk++;
  int nx=atoi(mychar);

  int mkk=0;
  mychar=malloc(len_ysize*sizeof(char));
  while (dim[kk]!=' ' && dim[kk]!='\n')
  {
	mychar[mkk]=dim[kk];
	kk++;
	mkk++;
  }
  kk++;
  int ny=atoi(mychar);


  mkk=0;
  mychar=malloc(len_maxval*sizeof(char));
  while (dim[kk]!=' ' && dim[kk]!='\n')
  {
        mychar[mkk]=dim[kk];
        kk++;
	mkk++;
  }
  kk++;
  int maxval=atoi(mychar);

  printf("%d %d\n%d\n", nx, ny, maxval);


  //for (int jj=0; jj<len_xsize; ++jj)  printf("mychar[%d] = %c\n", jj, mychar[jj] );
  //printf("nnx = %d, nny=%d\nmaxval=%d\n", nx, ny, maxval);



  //dim=(const char*) dim;
  //sscanf(&dim, "%d%*c%d%*c%d", xsize, ysize, maxval);
  //printf("%d, %d, %d\n", *xsize, *ysize, *maxval);
   

 // xsize=(unsigned short int*) dim;
 // printf("%hn", *xsize);
  
  //for (int jj=0; jj<count; ++jj)  printf("dim[%d] = %c\n", jj, dim[jj] );

  short int n;
  // while(fread(&n,sizeof(short int),1,f)==1) { // up to EOF
  //  printf("%d ",__bswap_16(n));
  //}
  
  int read_from=atoi(argv[2]);  // def 0 
  int read_up_to=atoi(argv[3]);  // def ny
  // print clean matrix of gray density
  unsigned short int* nn = malloc(nx*ny*sizeof(short int)); // allocate nx*ny*sizeof(short int) bytes
  fread(nn, sizeof(short int), nx*ny, f); // just one fread for the whole clean matrix
  for (int ii=read_from; ii<read_up_to; ++ii){  //rows
	for (int jj=nx-100; jj<nx; ++jj){
		printf("%d ", __bswap_16(nn[ii*nx+jj])); // print pixel degree of gray one by one

  	}
	printf("\n");
  }  

  printf("ny=%d", ny);





  
  fclose(f);
  return 0;
}
