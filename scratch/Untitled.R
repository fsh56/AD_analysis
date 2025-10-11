#SBATCH --partition=tier2q
#SBATCH --time=6:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --error=job_errors_%A_%a.log
#SBATCH --output=job_outputs_%A_%a.log
#SBATCH --array=574,575,578-1155
#SBATCH --nodes=1
module load gcc
module load R
Rscript code/R2_pm_clumpdf_fusios.R $disease param_0.csv $((SLURM_ARRAY_TASK_ID + 0)) GenomicSEM-ChildrenFathered-MatPatOff-cleaned.txt 0.001 200 0.1