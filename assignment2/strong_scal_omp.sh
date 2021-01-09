# -l nodes=1:ppn=24
# PBS -l walltime=24:00:00

cd $PBS_O_WORKDIR

# load gnu>4.9.9 to use proc_bind clause
module load gnu/9.3.0

echo THIS IS OPENMP RUNNING...



#export OMP_NUM_THREADS=3


export image_name=check_me.pgm
export kernel_dim=101
export final_image_name=omp_after.pgm


echo number of threads=$OMP_NUM_THREADS
echo image name=$image_name
echo kernel=$kernel_dim





# compiling in openmp 
gcc -O3 -std=gnu99 -fopenmp omp_blur_scal.c -o omp_blur_scal.x
./omp_blur_scal.x $image_name $kernel_dim





image_name=earth-large.pgm
mpi_cfile=omp_blur_scal  # C file to launch: without .c or .x
output_file=ww_omp.csv  # where results should be written
max_cores=24  # max cores to perform scalability on



#printf "kernel : ${kernel_dim}\nimage : ${image_name}\ncfile : ${mpi_cfile}.c\n" 2>>${output_file} 1>>${output_file}



#delete previous version of the file
#rm ${output_file}

#script for blurring via mpi 
echo THIS IS MPI RUNNING ...   #2>>${output_file} 1>>${output_file}


# compiling with -O1 in mpi
mpicc -O1 ${mpi_cfile}.c -o ${mpi_cfile}.x


printf "#P,kernel_dim,iter1,iter2,iter3\n" 2>>${output_file} 1>>${output_file}  # Header .csv
#echo Using a kernel of dimension ${kernel_dim} 2>>${output_file} 1>>${output_file}
for kernel_dim in 11 101; do
	for procs in $(seq 1 ${max_cores}); do
		printf "${procs},${kernel_dim}"  2>>${output_file} 1>>${output_file}
		for i in 1 2 3; do
			printf "," 2>>${output_file} 1>>${output_file}
			mpirun  --mca btl '^openib' -np  ${procs}  ./${mpi_cfile}.x ${image_name} ${kernel_dim} 2>>${output_file}  1>>${output_file}
			
		done
		printf "\n" 2>>${output_file} 1>>${output_file}
	done
done


printf "\n\n"
lscpu 2>>${output_file}  1>>${output_file}














