#!/bin/bash
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=2:00:00
#PBS -q short_cpuQ
module load mpich-3.2
cd $PBS_O_WORKDIR
mpirun.actual -n 1 ./parallelOMP 4
