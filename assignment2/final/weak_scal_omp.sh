# -l nodes=1:ppn=24
# PBS -l walltime=08:00:00




### PLEASE NOTICE: this script was used for scalability, but has different input command-line args to the ones
###                currently supported by the blurring file. In fact, these inputs were adapted to meet the 
###                specific requirements for the assignment.














cd /u/dssc/s275995/Davide/mynewgit/assignment2


# load gnu>4.9.9 to use proc_bind clause
module load gnu/9.3.0

echo THIS IS OPENMP RUNNING...



image_name=wblack.pgm
mpi_cfile=omp_blur_scal  		# C file to launch: without .c or .x
output_file=ww_omp.csv  		# where results should be written
max_cores=24				# max cores to perform scalability on
creation_file=create_image              # file that generates the suitable image for n cores: the image is almost a square
default_left=4409                       # length of the generated image for one core
default_right=4409                      # height of the generated image for one core



#export OMP_NUM_THREADS=3


# compiling in openmp 
#gcc -O3 -std=gnu99 -fopenmp ${mpi_cfile}.c -o ${mpi_cfile}.x

#compiling create_image file
gcc -O3 ${creation_file}.c -o ${creation_file}.x -lm;  # linking math library


#printf "kernel : ${kernel_dim}\nimage : ${image_name}\ncfile : ${mpi_cfile}.c\n" 2>>${output_file} 1>>${output_file}



#delete previous version of the file
#rm ${output_file}



# compiling C file
gcc -O3 -std=gnu99 -fopenmp ${mpi_cfile}.c -o ${mpi_cfile}.x

printf "#threads,kernel_dim,iter1,iter2,iter3\n" 2>>${output_file} 1>>${output_file}  # Header .csv
#echo Using a kernel of dimension ${kernel_dim} 2>>${output_file} 1>>${output_file}
for kernel_dim in 11 101; do
	for procs in $(seq 1 ${max_cores}); do  # loop over threads
		./${creation_file}.x ${procs} ${default_left} ${default_right}  # need to generate a new image every time procs changes
		export OMP_NUM_THREADS=${procs}
		printf "${procs},${kernel_dim}"  2>>${output_file} 1>>${output_file}
		for i in 1 2 3; do
			printf "," 2>>${output_file} 1>>${output_file}
			./${mpi_cfile}.x ${image_name} ${kernel_dim} 2>>${output_file}  1>>${output_file}
			
		done
		printf "\n" 2>>${output_file} 1>>${output_file}
	done
done


#printf "\n\n" 2>>${output_file}  1>>${output_file}
echo $PBS_NODEFILE 2>>${output_file}  1>>${output_file}












