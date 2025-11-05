#!/bin/bash
#SBATCH --job-name=ld_clump_ba9
#SBATCH --output=/scratch/sfeng56/logs/ld_%x_chr%a.out
#SBATCH --error=/scratch/sfeng56/logs/ld_%x_chr%a.err
#SBATCH --time=3:00:00
#SBATCH --mem=4G
#SBATCH --array=1-22

# ==========================================================
# LD Clumping for Brain_Amygdala tissue - One GWAS per job
# Each array task processes one chromosome (1-22)
# ==========================================================

# Load required modules
module load gcc/12.1.0
module load R/4.5.1
module load lapack/3.11.0
module load atlas/3.10.3
module load plink/1.9

# ==================== CONFIGURATION ====================
# MODIFY THIS SECTION FOR EACH GWAS
GWAS_NAME="dlbany"  # Change this for each GWAS submission
# =======================================================

# Fixed parameters
TISSUE_NAME="Brain_Frontal_Cortex_BA9"
CLUMP_KB=50
CLUMP_R2=0.1
CLUMP_P=0.001

# Directories
BASE_DIR="/gpfs/data/gao-lab/people/Sihao"
INPUT_DIR="${BASE_DIR}/data/harmonized_data_Brain_Frontal_Cortex_BA9"
OUTPUT_DIR="${BASE_DIR}/AD/ba9/clumped_data"
R_SCRIPT="${BASE_DIR}/scripts/ld_clumping_amygdala.R"
LOG_DIR="/scratch/sfeng56/logs"

# Create directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Get chromosome number from array task ID
CHR=${SLURM_ARRAY_TASK_ID}

echo "=================================================="
echo "LD Clumping Job Started (Per Gene Mode)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=================================================="
echo ""
echo "Configuration:"
echo "  GWAS: ${GWAS_NAME}"
echo "  Tissue: ${TISSUE_NAME}"
echo "  Chromosome: ${CHR}"
echo ""
echo "Parameters:"
echo "  Window size: ${CLUMP_KB} kb"
echo "  R2 threshold: ${CLUMP_R2}"
echo "  P-value threshold: ${CLUMP_P}"
echo "  Mode: Per-gene clumping"
echo ""

# Define input file
INPUT_FILE="${INPUT_DIR}/eqtl_${GWAS_NAME}_${TISSUE_NAME}_${CHR}.txt.gz"

# Check if input file exists
if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: Input file not found: ${INPUT_FILE}"
    echo "Please check:"
    echo "  1. GWAS name is correct: ${GWAS_NAME}"
    echo "  2. File naming format: eqtl_{gwas}_{tissue}_{chr}.txt.gz"
    echo "  3. Input directory: ${INPUT_DIR}"
    exit 1
fi

# Check if R script exists
if [ ! -f "${R_SCRIPT}" ]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    exit 1
fi

echo "Input file: ${INPUT_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Starting per-gene LD clumping..."
echo "Start time: $(date)"
echo ""

# Run R script
Rscript ${R_SCRIPT} ${INPUT_FILE} ${OUTPUT_DIR} ${CLUMP_KB} ${CLUMP_R2} ${CLUMP_P}

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "SUCCESS: Completed for ${GWAS_NAME} chr${CHR}"
    echo "=================================================="
    
    # Format parameter strings for filenames
    PVAL_STR=$(printf "%s" ${CLUMP_P})
    R2_STR=$(printf "%s" ${CLUMP_R2})
    
    # Expected output files
    OUTPUT_FILE="${OUTPUT_DIR}/clumped_${GWAS_NAME}_${TISSUE_NAME}_chr${CHR}_w${CLUMP_KB}kb_p${PVAL_STR}_r${R2_STR}.txt.gz"
    SUMMARY_FILE="${OUTPUT_DIR}/gene_summary_${GWAS_NAME}_${TISSUE_NAME}_chr${CHR}_w${CLUMP_KB}kb_p${PVAL_STR}_r${R2_STR}.txt"
    
    # Print output file info
    if [ -f "${OUTPUT_FILE}" ]; then
        FILE_SIZE=$(du -h ${OUTPUT_FILE} | cut -f1)
        echo ""
        echo "Clumped data file: ${OUTPUT_FILE}"
        echo "  File size: ${FILE_SIZE}"
        
        # Check line count (for gzipped file)
        LINE_COUNT=$(zcat ${OUTPUT_FILE} | wc -l)
        echo "  Number of lines: ${LINE_COUNT}"
    else
        echo ""
        echo "WARNING: Clumped data file not found"
    fi
    
    if [ -f "${SUMMARY_FILE}" ]; then
        LINE_COUNT=$(wc -l < ${SUMMARY_FILE})
        echo ""
        echo "Gene summary file: ${SUMMARY_FILE}"
        echo "  Number of genes: $((LINE_COUNT - 1))"
        echo ""
        echo "=== Preview of Gene Summary ==="
        head -10 ${SUMMARY_FILE}
    else
        echo ""
        echo "WARNING: Gene summary file not found"
    fi
else
    echo ""
    echo "=================================================="
    echo "ERROR: Failed for ${GWAS_NAME} chr${CHR}"
    echo "=================================================="
    exit 1
fi

echo ""
echo "End time: $(date)"
echo "=================================================="

