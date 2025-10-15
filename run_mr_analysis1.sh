#!/bin/bash
#SBATCH --job-name=mr_analysis1
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/mr_analysis1_%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/mr_analysis1_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G

# Run MR
INPUT_FILE="/gpfs/data/gao-lab/people/Sihao/ad_analysis/clumped_data/single_genes/clumped_cogdx_ad_Brain_Amygdala_9_ENSG00000106714.18_20kb_p0.01.txt.gz"

# Extract gene name from filename for output prefix
FILENAME=$(basename "$INPUT_FILE")
GENE_NAME=$(echo "$FILENAME" | grep -oP 'ENSG\d+\.\d+')

# Output directory and file prefix
OUTPUT_DIR="/gpfs/data/gao-lab/people/Sihao/data/mr_results_draft"
OUTPUT_PREFIX="${OUTPUT_DIR}/${GENE_NAME}"

# R script location
R_SCRIPT="/gpfs/data/gao-lab/people/Sihao/draft/mr_analysis.R"

# create output directory
mkdir -p "$OUTPUT_DIR"
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create output directory: $OUTPUT_DIR"
    exit 1
fi
# check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi
# check if R script exists
if [ ! -f "$R_SCRIPT" ]; then
    echo "ERROR: R script not found: $R_SCRIPT"
    exit 1
fi

echo "running MR analysis for gene: $GENE_NAME..."

# Load modules
module load gcc/12.1.0
module load R/4.5.1

# Run R script
Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_PREFIX"

echo "All results saved in: $OUTPUT_PREFIX"
echo "All Done"
echo ""