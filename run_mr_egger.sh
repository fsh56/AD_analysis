#!/bin/bash
#SBATCH --job-name=mr_analysis
#SBATCH --output=/scratch/sfeng56/ad_analysis/logs/mr_%A_%a.out
#SBATCH --error=/scratch/sfeng56/ad_analysis/logs/mr_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-3

INPUT_DIR="/scratch/sfeng56/ad_analysis/clumped_data2"
OUTPUT_DIR="/scratch/sfeng56/ad_analysis/mr_results_multiP"
R_SCRIPT="/scratch/sfeng56/draft/mr_egger_ivw.R"


# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Get list of input files
FILES=(${INPUT_DIR}/clumped_eqtl_dlbany_Brain_Frontal_Cortex_BA9_1_20kb_r20.01_p*.txt.gz)

# Get the current file based on array task ID
CURRENT_FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}

# Extract base name without extension for output prefix
BASENAME=$(basename ${CURRENT_FILE} .txt.gz)
OUTPUT_PREFIX="${OUTPUT_DIR}/${BASENAME}"

echo "=========================================="
echo "Processing file: ${CURRENT_FILE}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "=========================================="

# Run the R script
Rscript mr_egger_ivw.R ${CURRENT_FILE} ${OUTPUT_PREFIX}

echo "=========================================="
echo "Task completed for: ${BASENAME}"
echo "=========================================="