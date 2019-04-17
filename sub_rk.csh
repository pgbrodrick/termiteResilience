#!/bin/bash
#SBATCH --job-name=term_rk
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH -o logs/o
#SBATCH -e logs/e

#module load r/3.5.1

Rscript calc_rk.r ${1} ${2} ${3}
