module load stack/2024-06  gcc/12.2.0
module load r/4.4.0
export OMP_NUM_THREADS=1


for i in {1..19}
do 
  sbatch --time=150:00:00 --job-name="anchorG $i" --mem-per-cpu=2GB --output=outfiles/anchorG_$i.out --cpus-per-task=100 --wrap "Rscript --vanilla R/anchorG_CV.R $i"
done