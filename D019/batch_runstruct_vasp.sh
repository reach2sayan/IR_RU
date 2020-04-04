#!/bin/bash

#SBATCH --mem=100G

runstruct_vasp -lu -w vaspf.wrap srun --mpi=pmi2 -n 240
