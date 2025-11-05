#!/bin/bash
# =============================================================================
# Calculate Genomic Inflation Factor (Lambda) for MR Results
# =============================================================================

# Load R module
module load gcc/12.1.0
module load R/4.5.1

# Directories
BASE_DIR="/gpfs/data/gao-lab/people/Sihao"
MR_DIR="${BASE_DIR}/AD/ba9/mr_results"
R_SCRIPT="${BASE_DIR}/scripts/calculate_lambda.R"

# Output file (optional - if not specified, will save to mr_results/lambda_summary.csv)
OUTPUT_FILE="${MR_DIR}/lambdas/lambda_summary.csv"

echo "=========================================="
echo "Calculating Genomic Inflation Factors"
echo "Date: $(date)"
echo "=========================================="
echo ""
echo "MR results directory: ${MR_DIR}"
echo "R script: ${R_SCRIPT}"
echo "Output file: ${OUTPUT_FILE}"
echo ""

# Check if R script exists
if [ ! -f "${R_SCRIPT}" ]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    exit 1
fi

# Check if MR directory exists
if [ ! -d "${MR_DIR}" ]; then
    echo "ERROR: MR results directory not found: ${MR_DIR}"
    exit 1
fi

# Count MR result files
MR_FILES=$(ls ${MR_DIR}/*_*_Brain_Amygdala_w*.csv 2>/dev/null | grep -v "all_gwas" | wc -l)
echo "Found ${MR_FILES} MR result files"
echo ""

if [ ${MR_FILES} -eq 0 ]; then
    echo "ERROR: No MR result files found in ${MR_DIR}"
    echo "Expected format: {method}_{gwas}_{tissue}_w{kb}kb_p{pval}_r{r2}.csv"
    exit 1
fi

# Run R script
echo "Running lambda calculation..."
echo ""

Rscript ${R_SCRIPT} ${MR_DIR} ${OUTPUT_FILE}

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "SUCCESS: Lambda calculation completed"
    echo "=========================================="
    echo ""
    echo "Output file: ${OUTPUT_FILE}"
    
    if [ -f "${OUTPUT_FILE}" ]; then
        echo ""
        echo "Preview of results:"
        head -10 ${OUTPUT_FILE} | column -t -s','
        echo ""
        echo "Total rows: $(tail -n +2 ${OUTPUT_FILE} | wc -l)"
    fi
else
    echo ""
    echo "=========================================="
    echo "ERROR: Lambda calculation failed"
    echo "=========================================="
    exit 1
fi

echo ""
echo "Completed at $(date)"
echo "=========================================="