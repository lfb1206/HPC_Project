#!/bin/bash
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -l walltime=4:00:00
#PBS -q short_cpuQ
module load mpich-3.2
cd $PBS_O_WORKDIR
mpirun.actual -n 1 ./parallelOMP 1