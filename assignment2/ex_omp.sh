# .sh script for fast processing blurring with openmp and loading
# the resulting image on github to view it in eog on my pc


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

# upload on my_git 
git add $final_image_name 
git commit -m "checking omp blurring $final_image_name" 
git push







# compiling in openmp
gcc -O3 -std=gnu99 -fopenmp omp_blur_scal.c -o omp_blur_scal.x
./omp_blur_scal.x $image_name $kernel_dim


# upload on my_git
git add $final_image_name
git commit -m "checking omp blurring $final_image_name"
git push











