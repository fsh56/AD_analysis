#!/bin/bash
#SBATCH --job-name=harmonize_gwas_eqtl
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/logistic_harmonize_%A_%a.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/logistic_harmonize_%A_%a.err
#SBATCH --time=5:00:00
#SBATCH --mem=32G
#SBATCH --partition=tier1q
#SBATCH --array=1-4

# make log directory
mkdir -p /gpfs/data/gao-lab/people/Sihao/ad_analysis/logs

# load modules
module load gcc/12.1.0
module load R/4.5.1

# define GWAS files array
GWAS_FILES=(
    "cogdx_ad.assoc.logistic.gz"
    "hspath_typ.assoc.logistic.gz"
    "tdp_st4_binary.assoc.logistic.gz"
    "dlbany.assoc.logistic.gz"
    "dcfdx_ad.assoc.logistic.gz"
)

# define tissue
TISSUE="Brain_Amygdala"

# get current file based on array task ID
GWAS_FILE=${GWAS_FILES[$SLURM_ARRAY_TASK_ID]}

echo "Processing: ${GWAS_FILE}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Tissue: ${TISSUE}"

# run R file
Rscript harm_logistic.R ${GWAS_FILE} ${TISSUE}

echo "Completed: ${GWAS_FILE}"