#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name="prepare data"
#SBATCH --mem-per-cpu=50GB
#SBATCH --output=outfiles/prep.out
#SBATCH --cpus-per-task=1

module load stack/2024-06  gcc/12.2.0
module load r/4.4.0

export OMP_NUM_THREADS=1
# stop at the first failing script
set -e

Rscript --vanilla R/preprocessing/data_preprocessing.R
Rscript --vanilla R/preprocessing/data_preparation.R
Rscript --vanilla R/preprocessing/lagged_time.R
