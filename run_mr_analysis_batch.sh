#!/bin/bash
#SBATCH --job-name=mr_analysis_batch
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/mr_analysis_batch_%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/mr_analysis_batch_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

# Define directories
INPUT_DIR="/gpfs/data/gao-lab/people/Sihao/ad_analysis/clumped_data"
OUTPUT_DIR="/gpfs/data/gao-lab/people/Sihao/data/mr_results_draft"
R_SCRIPT="/gpfs/data/gao-lab/people/Sihao/draft/mr_analysis.R"
LOG_DIR="/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create output directory: $OUTPUT_DIR"
    exit 1
fi

# Check if R script exists
if [ ! -f "$R_SCRIPT" ]; then
    echo "ERROR: R script not found: $R_SCRIPT"
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Load modules
module load gcc/12.1.0
module load R/4.5.1

# Initialize counters
TOTAL_FILES=0
SUCCESS_COUNT=0
FAILED_COUNT=0

echo "========================================="
echo "MR Analysis Batch Processing"
echo "========================================="
echo "Start time: $(date)"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "R script: $R_SCRIPT"
echo "========================================="
echo ""

# Process all .txt.gz files in the input directory
for INPUT_FILE in ${INPUT_DIR}/*.txt.gz; do
    
    # Check if any files match the pattern
    if [ ! -f "$INPUT_FILE" ]; then
        echo "WARNING: No .txt.gz files found in $INPUT_DIR"
        break
    fi
    
    TOTAL_FILES=$((TOTAL_FILES + 1))
    
    # Extract filename
    FILENAME=$(basename "$INPUT_FILE")
    
    # Extract gene ID
    GENE_ID=$(echo "$FILENAME" | grep -oP 'ENSG\d+\.\d+' | head -1)
    
    # If gene ID extraction fails, try alternative pattern without version number
    if [ -z "$GENE_ID" ]; then
        GENE_ID=$(echo "$FILENAME" | grep -oP 'ENSG\d+' | head -1)
    fi
    
    # If still no gene ID found, use filename without extension
    if [ -z "$GENE_ID" ]; then
        GENE_ID=$(echo "$FILENAME" | sed 's/\.txt\.gz$//' | sed 's/^clumped_//')
        echo "WARNING: Could not extract ENSG ID from $FILENAME, using: $GENE_ID"
    fi
    
    # Output prefix
    OUTPUT_PREFIX="${OUTPUT_DIR}/${GENE_ID}"
    
    echo "----------------------------------------"
    echo "Processing file $TOTAL_FILES: $FILENAME"
    echo "Gene ID: $GENE_ID"
    echo "Input: $INPUT_FILE"
    echo "Output prefix: $OUTPUT_PREFIX"
    echo "Time: $(date)"
    
    # Run R script and capture output
    Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_PREFIX" 2>&1
    
    # Check exit status
    if [ $? -eq 0 ]; then
        echo "✓ Successfully completed: $GENE_ID"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "✗ Failed: $GENE_ID"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
    echo ""
done

echo "Results saved in: $OUTPUT_DIR"
echo ""
echo "All Done!"