#!/bin/bash
#SBATCH --job-name=mr_analysis_batch
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/mr_analysis_batch_%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/mr_analysis_batch_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G

INPUT_DIR="/gpfs/data/gao-lab/people/Sihao/ad_analysis/clumped_data/single_genes"
OUTPUT_DIR="/gpfs/data/gao-lab/people/Sihao/data/mr_results_draft"

R_SCRIPT="/gpfs/data/gao-lab/people/Sihao/draft/mr_analysis.R"

# Create output directory
mkdir -p "$OUTPUT_DIR"
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create output directory: $OUTPUT_DIR"
    exit 1
fi

# Check if R script exists
if [ ! -f "$R_SCRIPT" ]; then
    echo "ERROR: R script not found: $R_SCRIPT"
    exit 1
fi

# Load modules
module load gcc/12.1.0
module load R/4.5.1

# Process all files matching the pattern
for INPUT_FILE in ${INPUT_DIR}/clumped_cogdx_ad_Brain_Amygdala_9_ENSG*.txt.gz; do
    
    # Check if file exists
    if [ ! -f "$INPUT_FILE" ]; then
        echo "WARNING: File not found: $INPUT_FILE"
        continue
    fi
    
    # Extract gene name from filename
    FILENAME=$(basename "$INPUT_FILE")
    GENE_NAME=$(echo "$FILENAME" | grep -oP 'ENSG\d+\.\d+')
    
    # Output prefix
    OUTPUT_PREFIX="${OUTPUT_DIR}/${GENE_NAME}"
    
    echo "processing gene: $GENE_NAME"
    echo "Input file: $INPUT_FILE"
    
    # Run R script
    Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_PREFIX"
    
    if [ $? -eq 0 ]; then
        echo "Successfully completed: $GENE_NAME"
    else
        echo "Failed: $GENE_NAME"
    fi
    echo ""
done

echo "Results saved in: $OUTPUT_DIR"
echo "All Done!"
