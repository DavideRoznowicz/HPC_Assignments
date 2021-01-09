#include <stdio.h>
#include <stdlib.h>
#include <byteswap.h>
int main(int argc,char* argv[]) {
  
  FILE* f;
  f=fopen(argv[1],"r");
  unsigned short int *c;
  char type[2];
  int stop=0;
  int k=0;
  c=(unsigned short int*)calloc( (150*150), sizeof(unsigned short int) );
  // while(k<2000) {
  fread(c, sizeof(unsigned short int),(150*150) , f);
  for (int ss; ss<(150*150); ++ss)	printf("c[%d] = %hd\n",ss,*(c+ss));
  //	k++;
  // }
  fclose(f);
  return 0;
}
