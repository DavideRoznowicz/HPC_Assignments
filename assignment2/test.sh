#!/bin/bash
#PBS -M davide.roznowicz@gmail.com
#PBS -q dssc
#PBS -l walltime=01:00:00



cd /u/dssc/s275995/Davide/mynewgit/assignment2 


echo finally 1>>test.txt

for i in $(seq 1 1000000); do
	printf "ciao: ${i}\n" 1>>test.txt 2>>test.txt

done

echo $PBS_NODEFILE 2>>test.txt 1>>test.txt

