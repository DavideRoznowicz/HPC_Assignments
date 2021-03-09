#!/bin/bash
#PBS -l nodes=1:ppn=48
#PBS -l walltime=24:00:00

cd $PBS_O_WORKDIR
module load   openmpi/4.0.3/gnu/9.3.0
mpicc mpi_pi.c -o mpi_pi.x


for MOVES in 100000000 1000000000 10000000000 100000000000; do
 for i in 1 2 3; do
  for procs in 1 4 8 12 16 20 24 28 32 36 40 44 48 ; do
   echo "executing on ", ${procs}, "  processors"
   /usr/bin/time -f " ~%U usr %S system %E total %P CPU \n\n\n" mpirun  --mca btl '^openib' -np ${procs} ./mpi_pi.x  ${MOVES} 2>>ss_${MOVES}.${procs} 1>>ss_${MOVES}.${procs}
  done
 done
done
