#!/bin/bash

#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=126g

fitfc -si=str_relax.out -f -frnn=1.5