#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --job-name="fit"
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=outfiles/fit.out
#SBATCH --cpus-per-task=100

module load stack/2024-06  gcc/12.2.0
module load r/4.4.0

export OMP_NUM_THREADS=1

Rscript --vanilla R/drugInteraction.R