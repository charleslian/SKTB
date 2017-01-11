#!/bin/sh
#PBS -N vasp
#PBS -q vasp
#PBS -l nodes=1:ppn=6

# Change into the working directory
cd $PBS_O_WORKDIR    

# Calculate the number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
slater-koster > TBresult
