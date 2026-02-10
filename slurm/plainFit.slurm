#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --job-name="plainFit"
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=outfiles/plainFit.out
#SBATCH --cpus-per-task=100

module load stack/2024-06  gcc/12.2.0
module load r/4.4.0

export OMP_NUM_THREADS=1

Rscript --vanilla R/plainEstimation.R
