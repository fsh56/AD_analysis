#!/bin/bash
# =============================================================================
# Single Chromosome Test Script
# =============================================================================
# Use this script to test the analysis on a single chromosome before
# submitting the full array job
# =============================================================================

# Configuration
BASE_DIR="/gpfs/data/gao-lab/people/Sihao/amyloid_ba9"
R_SCRIPT="/scratch/sfeng56/draft/mr_analysis_ba9.R"
CHR=1  # Test with chromosome 1

# Load R module (adjust based on your system)
module load gcc/12.1.0
module load R/4.5.1

echo "=========================================="
echo "Testing MR Analysis on Chromosome ${CHR}"
echo "=========================================="
echo ""

# Check if input file exists
INPUT_FILE="${BASE_DIR}/ldByGene${CHR}/BA9_clumped_w50kb_p0.001.txt"
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    echo "Please check the BASE_DIR and chromosome number."
    exit 1
fi

echo "✓ Input file found: $INPUT_FILE"
echo ""

# Run the analysis
echo "Starting analysis..."
Rscript ${R_SCRIPT} ${CHR} ${BASE_DIR}

# Check results
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Test completed successfully!"
    echo "=========================================="
    echo ""
    echo "Output files should be in:"
    echo "  ${BASE_DIR}/ldByGene${CHR}/"
    echo ""
    echo "Expected files:"
    echo "  - BA9_chr${CHR}_IVW_results.csv"
    echo "  - BA9_chr${CHR}_Egger_results.csv"
    echo "  - BA9_chr${CHR}_WeightedMedian_results.csv"
    echo ""
    
    # Check if output files exist
    OUTPUT_DIR="${BASE_DIR}/ldByGene${CHR}"
    echo "Checking output files..."
    for method in IVW Egger WeightedMedian; do
        FILE="${OUTPUT_DIR}/BA9_chr${CHR}_${method}_results.csv"
        if [ -f "$FILE" ]; then
            LINES=$(wc -l < "$FILE")
            echo "  ✓ $method: $FILE ($LINES lines)"
        else
            echo "  ✗ $method: File not found"
        fi
    done
    echo ""
    echo "If the test looks good, submit the full array job with:"
    echo "  sbatch run_mr_array.sh"
else
    echo ""
    echo "=========================================="
    echo "Test failed! Please check the errors above."
    echo "=========================================="
    exit 1
fi