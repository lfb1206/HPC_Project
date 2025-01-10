#!/bin/bash
#PBS -l select=2:ncpus=8:mem=32gb
#PBS -l walltime=2:00:00
#PBS -q short_cpuQ
module load mpich-3.2
cd $PBS_O_WORKDIR
mpirun.actual -n 4 ./parallelOMP 4