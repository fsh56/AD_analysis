#!/bin/bash
#SBATCH --job-name=mr_fusiomr_test
#SBATCH --output=/scratch/sfeng56/ad_analysis/logs/mr_by_gene_seso_%j.out
#SBATCH --error=/scratch/sfeng56/ad_analysis/logs/mr_by_gene_seso_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

# Input/Output directories
INPUT_FILE="/scratch/sfeng56/ad_analysis/clumped_data/clumped_eqtl_amyloid_Brain_Frontal_Cortex_BA9_1.txt.gz_20kb_p0.001.txt.gz"
OUTPUT_BASE="/scratch/sfeng56/data/mr_results_fusiomr"
R_SCRIPT="/scratch/sfeng56/draft/mr_by_gene_seso.R"

# Create directories
mkdir -p "$OUTPUT_BASE"
mkdir -p "/scratch/sfeng56/ad_analysis/logs"

# Load modules
module load gcc/12.1.0
module load R/4.5.1

echo "========================================"
echo "FusioMR-SESO Test Run"
echo "========================================"
echo "Start time: $(date)"
echo ""

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    echo "Please check the file path and name."
    exit 1
fi

BASENAME=$(basename "$INPUT_FILE" .txt.gz)
TISSUE_NAME="amyloid_Brain_Frontal_Cortex_BA9_1"

OUTPUT_DIR="${OUTPUT_BASE}/${TISSUE_NAME}"
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "Processing file: $BASENAME"
echo "Tissue: $TISSUE_NAME"
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "========================================"
echo ""

# Run the R script
Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_DIR"

echo ""
echo "========================================"
echo "Checking results..."
echo "========================================"

# Check if any FusioMR results exist
FUSIO_COUNT=$(find "$OUTPUT_DIR" -name "*_fusiomr.csv" | wc -l)

if [ $FUSIO_COUNT -eq 0 ]; then
    echo "Warning: No FusioMR results found!"
    echo "Check the error log for details."
    exit 1
else
    echo "Success! Found $FUSIO_COUNT FusioMR result files"
    
    # Create summary file for this tissue
    find "$OUTPUT_DIR" -name "*_fusiomr.csv" -exec head -1 {} \; -quit > "${OUTPUT_DIR}/summary_fusiomr.csv"
    find "$OUTPUT_DIR" -name "*_fusiomr.csv" -exec tail -n +2 {} \; >> "${OUTPUT_DIR}/summary_fusiomr.csv"
    
    # Count total genes
    TOTAL_GENES=$(tail -n +2 "${OUTPUT_DIR}/summary_fusiomr.csv" | wc -l)
    echo "Total genes analyzed: $TOTAL_GENES"
    
    # Print summary statistics
    echo ""
    echo "Summary Statistics:"
    echo "-------------------"
    
    # Count significant results (p < 0.05)
    SIG_COUNT=$(tail -n +2 "${OUTPUT_DIR}/summary_fusiomr.csv" | awk -F',' '$6 < 0.05' | wc -l)
    echo "Significant results (p < 0.05): $SIG_COUNT"
    
    # Show first few results
    echo ""
    echo "First 5 results:"
    head -6 "${OUTPUT_DIR}/summary_fusiomr.csv" | column -t -s','
fi

echo ""
echo "========================================"
echo "Test Run Complete!"
echo "End time: $(date)"
echo "========================================"
echo "Results saved in: $OUTPUT_DIR"
echo "Summary file: ${OUTPUT_DIR}/summary_fusiomr.csv"