#!/bin/bash
#SBATCH --job-name=ld_clump_gene
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/ld_clump_gene_%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/ld_clump_gene_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --partition=tier2q

# Create log directory if it doesn't exist
mkdir -p /gpfs/data/gao-lab/people/Sihao/ad_analysis/logs

# Load required modules
module load gcc/12.1.0
module load R/4.5.1
module load lapack/3.11.0
module load atlas/3.10.3
module load plink/1.9

# ========================================
# Parameters (EDIT HERE)
# ========================================
GENE_ID="ENSG00000173253.16"
#GENE_ID="ENSG00000130559.19"
CLUMP_KB=20 
CLUMP_R2=0.1
CLUMP_P=0.01

# Define paths
INPUT_DIR="/gpfs/data/gao-lab/people/Sihao/data/harmonized_data_Brain_Amygdala"
OUTPUT_DIR="/gpfs/data/gao-lab/people/Sihao/ad_analysis/clumped_data/single_genes"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Input and output files
INPUT_FILE="${INPUT_DIR}/eqtl_cogdx_ad_Brain_Amygdala_9.txt.gz"
OUTPUT_FILE="${OUTPUT_DIR}/clumped_cogdx_ad_Brain_Amygdala_9_${GENE_ID}_${CLUMP_KB}kb_p${CLUMP_P}.txt.gz"

echo "Job ID: $SLURM_JOB_ID"
echo "Input file: $INPUT_FILE"
echo "Gene ID: $GENE_ID"
echo "Output file: $OUTPUT_FILE"
echo ""
echo "LD Clumping Parameters:"
echo "  clump_kb: $CLUMP_KB"
echo "  clump_r2: $CLUMP_R2"
echo "  clump_p: $CLUMP_P"
echo ""

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    exit 1
fi

# Run R script with parameters
Rscript clumping_per_gene.R ${INPUT_FILE} ${GENE_ID} ${OUTPUT_FILE} ${CLUMP_KB} ${CLUMP_R2} ${CLUMP_P}

echo "All done!"
echo ""
