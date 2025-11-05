#!/bin/bash
# Extract significant genes from MR results

# Load R
module load gcc/12.1.0
module load R/4.5.1

# Directories
BASE_DIR="/gpfs/data/gao-lab/people/Sihao"
MR_DIR="${BASE_DIR}/AD/amygdala/mr_results"
OUTPUT_DIR="${BASE_DIR}/AD/amygdala/gene_list"
R_SCRIPT="${BASE_DIR}/scripts/extract_significant_genes.R"

echo "=========================================="
echo "Extracting Significant Genes (p < 0.05)"
echo "=========================================="
echo ""
echo "MR results: ${MR_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# Run R script
Rscript ${R_SCRIPT} ${MR_DIR} ${OUTPUT_DIR}

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "SUCCESS"
    echo "=========================================="
    echo ""
    echo "Output files in: ${OUTPUT_DIR}"
    ls -lh ${OUTPUT_DIR}
else
    echo ""
    echo "ERROR: Script failed"
    exit 1
fi
