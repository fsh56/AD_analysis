#!/bin/bash
#SBATCH --job-name=p0001
#SBATCH --output=/scratch/sfeng56/ad_analysis/logs/mr_p001_%j.out
#SBATCH --error=/scratch/sfeng56/ad_analysis/logs/mr_p001_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# Directories and files
INPUT_FILE="/scratch/sfeng56/ad_analysis/clumped_data/clumped_eqtl_dlbany_Brain_Frontal_Cortex_BA9_1_20kb_r20.01_p0.001.txt.gz"
OUTPUT_DIR="/scratch/sfeng56/data/mr_results_p001"
R_SCRIPT="/scratch/sfeng56/draft/mr.R"

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "/scratch/sfeng56/ad_analysis/logs"

# Load modules
module load gcc/12.1.0
module load R/4.5.1

echo "======================================================================"
echo "Starting MR analysis for p=0.00x1"
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Date: $(date)"
echo "======================================================================"
echo ""

# Check if file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found!"
    exit 1
fi

# Run R script
Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_DIR"

EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "======================================================================"
    echo "✓ Analysis completed successfully!"
    echo "Results saved in: $OUTPUT_DIR"
    echo "End time: $(date)"
    echo "======================================================================"
    echo ""
    echo "Output files:"
    ls -lh "${OUTPUT_DIR}"/*.csv 2>/dev/null || echo "No CSV files found"
else
    echo ""
    echo "======================================================================"
    echo "✗ Analysis failed with exit code: $EXIT_CODE"
    echo "Check error log for details"
    echo "======================================================================"
    exit $EXIT_CODE
fi
