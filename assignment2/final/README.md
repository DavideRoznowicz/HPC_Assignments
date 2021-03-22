# Directory organization

File organization in the folder: The files that were specifically requested in [Assignment2.pdf](../Assignment2.pdf) and [Final details for Assignment 2.pdf](../Final\details\for\Assignment\2.pdf) (with a specific name required) are present with the right name.
Some additional files are present, together with some others that were requested but without any name indication:

- create_image.c		file used to create a new image of specific size, proportionate to the number of cores/threads (used for weak scalability)
- strong_scal_mpi.sh		script for strong scalability for blur.mpi.c
- strong_scal_omp.sh		script for strong scalability for blur.omp.c
- weak_scal_mpi.sh		script for weak scalability for blur.mpi.c
- weak_scal_omp.sh		script for weak scalability for blur.omp.c
- ss_mpi.csv			file containing the results of the strong scalability for mpi
- ss_omp.csv			file containing the results of the strong scalability for omp
- ww_mpi.csv			file containing the results of the weak scalability for mpi
- ww_omp.csv			file containing the results of the weak scalability for omp






## Example (bash script) for compiling and then running MPI/OpenMP:


procs=4				# number of cores\
kernel_type=1			# type of kernel: average, weight, gaussian - 0,1,2\
kernel_size=101			# dimensions of the square kernel\
ff='0.2'			# ff should be specified only if kernel_type=1, otherwise it must not\
image=test_picture.pgm		# input image\
output_file="image_after.pgm"	# if output_file is specified, then the name for the image in output will be exactly <output_file>\
				# otherwise the name is univocally determined according to the assignment instructions


## Compiling/Running MPI:

module load openmpi/4.0.3/gnu/9.3.0\
mpicc -O1 blur.mpi.c -o blur.mpi -lm\
mpirun --mca btl ^openib -np ${procs} ./blur.mpi $kernel_type $kernel_size $ff $image #$output_file



## Compiling/Running OpenMP:

module load openmpi/4.0.3/gnu/9.3.0\
export OMP_NUM_THREADS=5\
gcc -O3 -std=gnu99 -fopenmp blur.omp.c -o blur.omp -lm\
./blur.omp $kernel_type $kernel_size $ff $image #$output_file

















