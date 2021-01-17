#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <byteswap.h>
#include <math.h>
#include <sched.h>
#include <stdbool.h>
#include <math.h>

/*
	compiling: 		gcc create_image.c -o create_image.x -lm
	running  : 		./create_image.x 2 100 100

	reading matrix:     ./leggi.x wblack.pgm
	viewing matrix:     eog wblack.pgm
	view header   :	hexdump -C wblack.pgm |less

*/


void create_image(int ncores, int nx, int ny) {
	
	
	int maxval=65535;
	int old_pixels=nx*ny;
		
	int new_nx=ceil(sqrt(ncores*nx*ny));
	int new_ny=floor(sqrt(ncores*nx*ny));
	
	int new_pixels=new_nx*new_ny;
	unsigned short int* m=(unsigned short int*)malloc(new_pixels*sizeof(short int));
	
	
	FILE* fc;
	fc=fopen("wblack.pgm","w");
	fprintf(fc,"P5\n%d %d\n%d\n", new_nx, new_ny, maxval);
	fclose(fc);
	fc=fopen("wblack.pgm","ab");
	
	for (int i=0; i<new_ny; i++){
		for (int j=0; j<new_nx; j++){
			m[i*new_nx+j]=__bswap_16(floor((maxval-1)*drand48())); // random image
		}
	}
	
	fwrite(m, sizeof(short int), new_pixels, fc);
	fclose(fc);





}


int main(int argc,char* argv[]) {

	// if ncores==1, then the image generated is the one with nx, ny given
	// otherwise the dimensions are s.t. workload is proportional for the number of cores

	int ncores=atoi(argv[1]);    // number of cores
	int nx=atoi(argv[2]);        // nx for image for one core
	int ny=atoi(argv[3]);        // ny for image for one core
	
	
	create_image(ncores, nx, ny);  // the image is generated



return 0;
}
