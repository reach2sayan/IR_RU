#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --mem=48g

robustrelax_vasp -id -c 0.05 srun --mpi=pmi2 -n 64 
