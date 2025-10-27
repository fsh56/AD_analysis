#!/bin/bash
#SBATCH --job-name=ld_pergene
#SBATCH --output=/scratch/sfeng56/logs/ld_pergene_ba9_%A_%a.out
#SBATCH --error=/scratch/sfeng56/logs/ld_pergene_ba9_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G

# Load required modules
module load gcc/12.1.0
module load R/4.5.1
module load lapack/3.11.0
module load atlas/3.10.3
module load plink/1.9

# Set variables
INPUT_FILE="/gpfs/data/gao-lab/people/Sihao/data/harmonized_data_Brain_Frontal_Cortex_BA9/eqtl_amyloid_Brain_Frontal_Cortex_BA9_1.txt.gz"
OUTPUT_DIR="/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/ldByGene"
R_SCRIPT="/scratch/sfeng56/draft/ld_per_gene.R"
LOG_DIR="/scratch/sfeng56/logs"

# Create directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Define parameter combinations
# Format: "WINDOW_KB P_VALUE"
PARAMS=(
    "50 0.0005"
    "50 0.001"
)

# Get parameters for this array job
PARAM_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
PARAM_SET=${PARAMS[$PARAM_INDEX]}
WINDOW_KB=$(echo $PARAM_SET | awk '{print $1}')
PVAL=$(echo $PARAM_SET | awk '{print $2}')

echo "=================================================="
echo "LD Clumping Job Started (Per Gene Mode)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=================================================="
echo ""
echo "Parameters:"
echo "  Window size: ${WINDOW_KB} kb"
echo "  P-value threshold: ${PVAL}"
echo "  R2 threshold: 0.1 (fixed)"
echo "  Mode: Per-gene clumping"
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

echo "Starting per-gene LD clumping..."
echo "Start time: $(date)"
echo ""

# Run R script
Rscript ${R_SCRIPT} ${INPUT_FILE} ${OUTPUT_DIR} ${WINDOW_KB} ${PVAL}

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo "SUCCESS: Completed for Window=${WINDOW_KB}kb, P-value=${PVAL}"
    
    # Print output file info
    OUTPUT_FILE="${OUTPUT_DIR}/BA9_clumped_w${WINDOW_KB}kb_p${PVAL}.txt"
    SUMMARY_FILE="${OUTPUT_DIR}/BA9_gene_summary_w${WINDOW_KB}kb_p${PVAL}.txt"
    LOG_FILE="${OUTPUT_DIR}/BA9_clumping_log_w${WINDOW_KB}kb_p${PVAL}.txt"
    
    if [ -f "${OUTPUT_FILE}" ]; then
        LINE_COUNT=$(wc -l < ${OUTPUT_FILE})
        FILE_SIZE=$(du -h ${OUTPUT_FILE} | cut -f1)
        echo "Main output file: ${OUTPUT_FILE}"
        echo "  File size: ${FILE_SIZE}"
        echo "  Number of lines: ${LINE_COUNT}"
    fi
    
    if [ -f "${SUMMARY_FILE}" ]; then
        LINE_COUNT=$(wc -l < ${SUMMARY_FILE})
        echo "Gene summary file: ${SUMMARY_FILE}"
        echo "  Number of genes: $((LINE_COUNT - 1))"
    fi
    
    if [ -f "${LOG_FILE}" ]; then
        echo "Log file: ${LOG_FILE}"
        echo ""
        echo "=== Preview of Results ==="
        head -20 ${LOG_FILE}
    fi
else
    echo ""
    echo "ERROR: Failed for Window=${WINDOW_KB}kb, P-value=${PVAL}"
    exit 1
fi

echo ""
echo "End time: $(date)"
echo "=================================================="