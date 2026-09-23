#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --job-name="PAnal"
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=outfiles/PAnal.out
#SBATCH --cpus-per-task=100

module load stack/2024-06  gcc/12.2.0
module load r/4.4.0

export OMP_NUM_THREADS=1

Rscript --vanilla R/validationOfPvalues.R