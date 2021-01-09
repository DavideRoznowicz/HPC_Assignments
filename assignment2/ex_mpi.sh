# exec script for blurring via mpi



module load openmpi/4.0.3/gnu/9.3.0



echo THIS IS MPI RUNNING ...

procs=5
image=check_me.pgm


mpicc -O1 mpi_blur_scal.c -o mpi_blur_scal.x
mpirun --mca btl ^openib -np ${procs} ./mpi_blur_scal.x ${image} 101

# git upload
git add mpi_after.pgm; git commit -m "mpi_after.pgm"; git push





echo NOW mpi_blur_scl.c
mpicc -O1 mpi_blur_scal.c -o mpi_blur_scal.x
mpirun --mca btl ^openib -np ${procs} ./mpi_blur_scal.x ${image} 101


# git upload
git add mpi_after.pgm; git commit -m "mpi_after.pgm"; git push








