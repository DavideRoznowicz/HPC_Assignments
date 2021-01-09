#include <stdio.h>
#include <stdlib.h>
#include <byteswap.h>
int main() {
  
	
  int nx=15;
  int ny=15;
  int lev=65535;
short int* m=malloc(nx*ny*sizeof(short int));
 int count=0;
  for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++) {
      count++;
      m[j*nx+i]=count;
    }
  }
  
  FILE* fc;
  fc=fopen("sfumatura_big_endian.pgm","w");
  fprintf(fc,"P5\n%d %d\n%d\n", nx, ny,lev);
  fclose(fc);
  fc=fopen("sfumatura_big_endian.pgm","ab");
  for(long int j=0;j<ny*nx;j++) {
    short int tmp=m[j];

	tmp=__bswap_16(tmp);

	fwrite(&tmp,sizeof(short int),1,fc);
  }
  fclose(fc);
  

  fc=fopen("sfumatura_little_endian.pgm","w");
  fprintf(fc,"P5\n%d %d\n%d\n-1\n", nx, ny,lev);
  fclose(fc);
  fc=fopen("sfumatura_little_endian.pgm","ab");
  for(long int j=0;j<ny*nx;j++) {
    short int tmp=m[j];

	//tmp=__bswap_16(tmp);

	fwrite(&tmp,sizeof(short int),1,fc);
  }
  fclose(fc);
   

  fc=fopen("sfumatura_sbagliata.pgm","w");
  fprintf(fc,"P5\n%d %d\n%d\n", nx, ny,lev);
  fclose(fc);
  fc=fopen("sfumatura_sbagliata.pgm","ab");
  for(long int j=0;j<ny*nx;j++) {
    short int tmp=m[j];

	//tmp=__bswap_16(tmp);

	fwrite(&tmp,sizeof(short int),1,fc);
  }
  fclose(fc);
}
