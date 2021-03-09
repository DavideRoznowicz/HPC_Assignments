#!/bin/bash
#PBS -l nodes=1:ppn=48
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR
gcc pi.c -o pi.x



for MOVES in 1000000000 10000000000 100000000000; do 
  for i in 1 2 3; do
    /usr/bin/time -f " ~%U usr %S system %E total %P CPU \n\n\n" ./pi.x  ${MOVES} 2>>single_core_${MOVES}.txt 1>>single_core_${MOVES}.txt
  done
done
