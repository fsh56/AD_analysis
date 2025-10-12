#!/bin/bash
#SBATCH --job-name=harmonize_linear
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/harmonize_linear%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/harmonize_linear%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --partition=tier1q
#SBATCH --cpus-per-task=1
#SBATCH --array=0-2

# Create log directory if it doesn't exist
mkdir -p /gpfs/data/gao-lab/people/Sihao/ad_analysis/logs

# Load required modules
module load gcc/12.1.0
module load R/4.5.1

# define GWAS files array
GWAS_FILES=(
    "amyloid.assoc.linear.with_rsid.txt.gz"
    "tangles.assoc.linear.with_rsid.txt.gz"
    "gpath.assoc.linear.with_rsid.txt.gz"
)

# define tissue
TISSUE="Brain_Amygdala"

# get current file based on array task ID
GWAS_FILE=${GWAS_FILES[$SLURM_ARRAY_TASK_ID]}

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "GWAS file: $GWAS_FILE"
echo "Tissue: $TISSUE"
echo ""

# run R script
Rscript harm_linear_with_rsid.R ${GWAS_FILE} ${TISSUE}

echo "Done!"