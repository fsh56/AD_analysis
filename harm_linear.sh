#!/bin/bash
#SBATCH --job-name=harmonize_linear
#SBATCH --output=/scratch/sfeng56/logs/harmonize_linear%j.out
#SBATCH --error=/scratch/sfeng56/logs/harmonize_linear%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-4

# Load required modules
module load gcc/12.1.0
module load R/4.5.1

BASE_DIR="/gpfs/data/gao-lab/people/Sihao"
R_SCRIPT="${BASE_DIR}/scripts/harm_linear.R"

# define GWAS files array
GWAS_FILES=(
    "amyloid.assoc.linear.with_rsid.txt.gz"
    "tangles.assoc.linear.with_rsid.txt.gz"
    "gpath.assoc.linear.with_rsid.txt.gz"
)

# define tissue
TISSUE="Brain_Anterior_cingulate_cortex_BA24"

# get current file based on array task ID
GWAS_FILE=${GWAS_FILES[$SLURM_ARRAY_TASK_ID]}

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "GWAS file: $GWAS_FILE"
echo "Tissue: $TISSUE"
echo ""

# run R script
Rscript ${R_SCRIPT} ${GWAS_FILE} ${TISSUE}

echo "Done!"