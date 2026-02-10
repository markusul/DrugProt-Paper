#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --job-name="opt"
#SBATCH --mem-per-cpu=50GB
#SBATCH --output=outfiles/opt.out
#SBATCH --cpus-per-task=1

module load stack/2024-06  gcc/12.2.0
module load r/4.4.0

export OMP_NUM_THREADS=1

Rscript --vanilla R/anchorG_opt_res.R