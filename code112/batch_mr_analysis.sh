#!/bin/bash
#SBATCH --job-name=MR_amy
#SBATCH --output=/scratch/sfeng56/logs/MR_%x.out
#SBATCH --error=/scratch/sfeng56/logs/MR_%x.err
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# =============================================================================
# MR Analysis for Brain_Amygdala - One GWAS per job
# This script performs MR analysis on clumped data
# Input: All chromosomes for one GWAS
# Output: Three method-specific result files (IVW, Egger, WeightedMedian)
# =============================================================================

# Load required modules
module load gcc/12.1.0
module load R/4.5.1

# ==================== CONFIGURATION ====================
# MODIFY THIS SECTION FOR EACH GWAS
GWAS_NAME="dlbany"  # Change this for each GWAS submission
# =======================================================

# Fixed parameters
TISSUE_NAME="Brain_Amygdala"
CLUMP_KB=50
CLUMP_R2=0.1
CLUMP_P=0.001

# Directories
BASE_DIR="/gpfs/data/gao-lab/people/Sihao"
CLUMP_DIR="${BASE_DIR}/AD/amygdala/clumped_data"
OUTPUT_DIR="${BASE_DIR}/AD/amygdala/mr_results"
R_SCRIPT="${BASE_DIR}/scripts/mr_analysis_amygdala.R"
LOG_DIR="/scratch/sfeng56/logs"

# Create directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

echo "=========================================="
echo "MR Analysis Job Started"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date: $(date)"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  GWAS: ${GWAS_NAME}"
echo "  Tissue: ${TISSUE_NAME}"
echo ""
echo "Parameters:"
echo "  Window size: ${CLUMP_KB} kb"
echo "  R2 threshold: ${CLUMP_R2}"
echo "  P-value threshold: ${CLUMP_P}"
echo ""
echo "Directories:"
echo "  Clumped data: ${CLUMP_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo ""

# Check if R script exists
if [ ! -f "${R_SCRIPT}" ]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    exit 1
fi

# Check if clumped data directory exists
if [ ! -d "${CLUMP_DIR}" ]; then
    echo "ERROR: Clumped data directory not found: ${CLUMP_DIR}"
    exit 1
fi

# Format parameter strings
PVAL_STR=$(printf "%s" ${CLUMP_P})
R2_STR=$(printf "%s" ${CLUMP_R2})

# Check if at least one chromosome file exists
SAMPLE_FILE="${CLUMP_DIR}/clumped_${GWAS_NAME}_${TISSUE_NAME}_chr1_w${CLUMP_KB}kb_p${PVAL_STR}_r${R2_STR}.txt.gz"

if [ ! -f "${SAMPLE_FILE}" ]; then
    echo "ERROR: No clumped data files found for ${GWAS_NAME}"
    echo "Expected format: clumped_${GWAS_NAME}_${TISSUE_NAME}_chr*_w${CLUMP_KB}kb_p${PVAL_STR}_r${R2_STR}.txt.gz"
    echo "Sample file checked: ${SAMPLE_FILE}"
    echo ""
    echo "Available files in directory:"
    ls -1 ${CLUMP_DIR}/clumped_${GWAS_NAME}_*.txt.gz 2>/dev/null | head -5
    exit 1
fi

echo "Clumped data files found. Starting MR analysis..."
echo "Start time: $(date)"
echo ""

# Run the R script
Rscript ${R_SCRIPT} \
    "${GWAS_NAME}" \
    "${TISSUE_NAME}" \
    "${CLUMP_DIR}" \
    "${OUTPUT_DIR}" \
    ${CLUMP_KB} \
    ${CLUMP_P} \
    ${CLUMP_R2}

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "SUCCESS: MR analysis completed for ${GWAS_NAME}"
    echo "=========================================="
    
    # Expected output files
    IVW_FILE="${OUTPUT_DIR}/IVW_${GWAS_NAME}_${TISSUE_NAME}_w${CLUMP_KB}kb_p${PVAL_STR}_r${R2_STR}.csv"
    EGGER_FILE="${OUTPUT_DIR}/Egger_${GWAS_NAME}_${TISSUE_NAME}_w${CLUMP_KB}kb_p${PVAL_STR}_r${R2_STR}.csv"
    WM_FILE="${OUTPUT_DIR}/WeightedMedian_${GWAS_NAME}_${TISSUE_NAME}_w${CLUMP_KB}kb_p${PVAL_STR}_r${R2_STR}.csv"
    
    echo ""
    echo "Output files:"
    
    # Check IVW results
    if [ -f "${IVW_FILE}" ]; then
        LINE_COUNT=$(wc -l < ${IVW_FILE})
        FILE_SIZE=$(du -h ${IVW_FILE} | cut -f1)
        echo "  ✓ IVW: ${IVW_FILE}"
        echo "    Genes: $((LINE_COUNT - 1)), Size: ${FILE_SIZE}"
    else
        echo "  ✗ IVW: File not created"
    fi
    
    # Check Egger results
    if [ -f "${EGGER_FILE}" ]; then
        LINE_COUNT=$(wc -l < ${EGGER_FILE})
        FILE_SIZE=$(du -h ${EGGER_FILE} | cut -f1)
        echo "  ✓ Egger: ${EGGER_FILE}"
        echo "    Genes: $((LINE_COUNT - 1)), Size: ${FILE_SIZE}"
    else
        echo "  ✗ Egger: File not created"
    fi
    
    # Check Weighted Median results
    if [ -f "${WM_FILE}" ]; then
        LINE_COUNT=$(wc -l < ${WM_FILE})
        FILE_SIZE=$(du -h ${WM_FILE} | cut -f1)
        echo "  ✓ Weighted Median: ${WM_FILE}"
        echo "    Genes: $((LINE_COUNT - 1)), Size: ${FILE_SIZE}"
    else
        echo "  ✗ Weighted Median: File not created"
    fi
    
else
    echo ""
    echo "=========================================="
    echo "ERROR: MR analysis failed for ${GWAS_NAME}"
    echo "=========================================="
    echo "Please check the log file for details"
    exit 1
fi

echo ""
echo "End time: $(date)"
echo "=========================================="
