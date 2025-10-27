#!/bin/bash
#SBATCH --job-name=ld_clump_genes_parallel
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/ld_clump_genes_parallel_%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/ld_clump_genes_parallel_%j.err
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
# Define multiple gene IDs (space-separated)
GENE_IDS=(
  "ENSG00000130559.19" # include most SNPs
  "ENSG00000173253.16"
  "ENSG00000106714.18"
  "ENSG00000226007.3" # include least SNPs
  "ENSG00000288062.1" # sed -n 300p
)

CLUMP_KB=20 
CLUMP_R2=0.1
CLUMP_P=0.01

# Number of parallel jobs (should match --cpus-per-task)
MAX_PARALLEL=10

# Define paths
INPUT_DIR="/gpfs/data/gao-lab/people/Sihao/data/harmonized_data_Brain_Amygdala"
OUTPUT_DIR="/gpfs/data/gao-lab/people/Sihao/ad_analysis/clumped_data/single_genes"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Input file
INPUT_FILE="${INPUT_DIR}/eqtl_cogdx_ad_Brain_Amygdala_9.txt.gz"

echo "Job ID: $SLURM_JOB_ID"
echo "Input file: $INPUT_FILE"
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

# Function to process a single gene
process_gene() {
    local GENE_ID=$1
    local OUTPUT_FILE="${OUTPUT_DIR}/clumped_cogdx_ad_Brain_Amygdala_9_${GENE_ID}_${CLUMP_KB}kb_p${CLUMP_P}.txt.gz"
    
    echo "[$(date)] Starting: $GENE_ID"
    Rscript clumping_per_gene.R ${INPUT_FILE} ${GENE_ID} ${OUTPUT_FILE} ${CLUMP_KB} ${CLUMP_R2} ${CLUMP_P} 2>&1 | sed "s/^/[$GENE_ID] /"
    
    local EXIT_CODE=$?
    if [ $EXIT_CODE -eq 0 ]; then
        echo "[$(date)] ✓ Completed: $GENE_ID"
    else
        echo "[$(date)] ✗ Failed: $GENE_ID (exit code: $EXIT_CODE)"
    fi
    
    return $EXIT_CODE
}

# Export function and variables for parallel execution
export -f process_gene
export INPUT_FILE OUTPUT_DIR CLUMP_KB CLUMP_R2 CLUMP_P

# Process genes in parallel using GNU parallel or xargs
# Using background jobs with job control
echo "Starting parallel processing..."
echo ""

job_count=0
for GENE_ID in "${GENE_IDS[@]}"; do
    # Wait if we've reached max parallel jobs
    while [ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]; do
        sleep 1
    done
    
    # Start job in background
    process_gene "$GENE_ID" &
    ((job_count++))
done

# Wait for all background jobs to complete
echo ""
echo "Waiting for all jobs to complete..."
wait


echo "All genes processed!"
