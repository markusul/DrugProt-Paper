module load stack/2024-06  gcc/12.2.0
module load r/4.4.0
export OMP_NUM_THREADS=1

sbatch --time=5:00:00 --job-name="Z 1" --mem-per-cpu=40GB --output=outfiles/Z_1.out --cpus-per-task=1 --wrap "Rscript --vanilla R/getZ.R 1"
sbatch --time=155:00:00 --job-name="Z 2" --mem-per-cpu=200GB --output=outfiles/Z_2.out --cpus-per-task=1 --wrap "Rscript --vanilla R/getZ.R 2"
sbatch --time=125:00:00 --job-name="Z 3" --mem-per-cpu=200GB --output=outfiles/Z_3.out --cpus-per-task=1 --wrap "Rscript --vanilla R/getZ.R 3"
