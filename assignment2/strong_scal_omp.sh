# -l nodes=1:ppn=24
# PBS -l walltime=23:00:00

cd /u/dssc/s275995/Davide/mynewgit/assignment2


# load gnu>4.9.9 to use proc_bind clause
module load gnu/9.3.0

echo THIS IS OPENMP RUNNING...



image_name=earth-large.pgm
mpi_cfile=omp_blur_scal  		# C file to launch: without .c or .x
output_file=ss_omp.csv  		# where results should be written
max_cores=24				# max cores to perform scalability on


#export OMP_NUM_THREADS=3


# compiling in openmp 
#gcc -O3 -std=gnu99 -fopenmp ${mpi_cfile}.c -o ${mpi_cfile}.x


#printf "kernel : ${kernel_dim}\nimage : ${image_name}\ncfile : ${mpi_cfile}.c\n" 2>>${output_file} 1>>${output_file}



#delete previous version of the file
rm ${output_file}



printf "#threads,kernel_dim,iter1,iter2,iter3\n" 2>>${output_file} 1>>${output_file}  # Header .csv
#echo Using a kernel of dimension ${kernel_dim} 2>>${output_file} 1>>${output_file}
for kernel_dim in 11 101; do
	for procs in $(seq 1 ${max_cores}); do
		export OMP_NUM_THREADS=${procs}
		gcc -O3 -std=gnu99 -fopenmp ${mpi_cfile}.c -o ${mpi_cfile}.x
		printf "${procs},${kernel_dim}"  2>>${output_file} 1>>${output_file}
		for i in 1 2 3; do
			printf "," 2>>${output_file} 1>>${output_file}
			./${mpi_cfile}.x ${image_name} ${kernel_dim} 2>>${output_file}  1>>${output_file}
			
		done
		printf "\n" 2>>${output_file} 1>>${output_file}
	done
done


#printf "\n\n" 2>>${output_file}  1>>${output_file}
#lscpu 2>>${output_file}  1>>${output_file}












