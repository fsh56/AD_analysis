#!/bin/bash
#SBATCH --job-name=ld_clump
#SBATCH --output=logs/ld_clump_%A_%a.out
#SBATCH --error=logs/ld_clump_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-4

# LD Clumping Pipeline with ieugwasr
# This script runs LD clumping with different window sizes and p-value thresholds
# Based on original clumping_params.R approach
# Author: Updated script
# Date: 2025-10-21

# Set variables
INPUT_FILE="/scratch/sfeng56/data/harmonized_data_Brain_Amygdala/eqtl_amyloid_Brain_Amygdala_1.txt.gz"
OUTPUT_DIR="/scratch/sfeng56/data/harmonized_data_Brain_Amygdala/ld_clumping_results"
R_SCRIPT="ld_clumping_ieugwasr.R"
LOG_DIR="logs"

# Create directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Define parameter combinations
# Format: "WINDOW_KB P_VALUE"
PARAMS=(
    "20 0.001"
    "20 0.0001"
    "50 0.001"
    "50 0.0001"
)

# Get parameters for this array job
PARAM_SET=${PARAMS[$SLURM_ARRAY_TASK_ID-1]}
WINDOW_KB=$(echo $PARAM_SET | awk '{print $1}')
PVAL=$(echo $PARAM_SET | awk '{print $2}')

echo "=================================================="
echo "LD Clumping Job Started"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=================================================="
echo ""
echo "Parameters:"
echo "  Window size: ${WINDOW_KB} kb"
echo "  P-value threshold: ${PVAL}"
echo "  R2 threshold: 0.1 (fixed)"
echo ""

# Check if input file exists
if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: Input file not found: ${INPUT_FILE}"
    exit 1
fi

# Check if R script exists
if [ ! -f "${R_SCRIPT}" ]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    exit 1
fi

# Load required modules (adjust based on your cluster)
# module load R/4.2.0
# module load plink/1.90

echo "Starting LD clumping..."
echo "Start time: $(date)"
echo ""

# Run R script
Rscript ${R_SCRIPT} ${INPUT_FILE} ${OUTPUT_DIR} ${WINDOW_KB} ${PVAL}

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo "SUCCESS: Completed for Window=${WINDOW_KB}kb, P-value=${PVAL}"
    
    # Print output file info
    OUTPUT_FILE="${OUTPUT_DIR}/clumped_w${WINDOW_KB}kb_p${PVAL}.txt"
    if [ -f "${OUTPUT_FILE}" ]; then
        LINE_COUNT=$(wc -l < ${OUTPUT_FILE})
        FILE_SIZE=$(du -h ${OUTPUT_FILE} | cut -f1)
        echo "Output file: ${OUTPUT_FILE}"
        echo "File size: ${FILE_SIZE}"
        echo "Number of lines: ${LINE_COUNT}"
    fi
else
    echo ""
    echo "ERROR: Failed for Window=${WINDOW_KB}kb, P-value=${PVAL}"
    exit 1
fi

echo ""
echo "End time: $(date)"
echo "=================================================="