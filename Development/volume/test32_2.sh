#!/bin/bash
#PBS -l select=4:ncpus=16:mem=32gb
#PBS -l walltime=1:00:00
#PBS -q short_cpuQ
module load mpich-3.2
cd $PBS_O_WORKDIR
mpirun.actual -n 32 ./parallelOMP 2