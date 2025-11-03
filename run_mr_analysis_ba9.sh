#!/bin/bash
#SBATCH --job-name=MRamy
#SBATCH --output=/scratch/sfeng56/logs/MR_amy_chr%a.out
#SBATCH --error=/scratch/sfeng56/logs/MR_amy_chr%a.err
#SBATCH --array=1-22
#SBATCH --time=4:00:00
#SBATCH --mem=16G

# =============================================================================
# Shell Script for Running MR Analysis on All Chromosomes
# =============================================================================
# This script submits array jobs for chromosomes 1-22
# Each chromosome is processed independently in parallel
# =============================================================================


BASE_DIR="/gpfs/data/gao-lab/people/Sihao/amyloid_amydgala"
R_SCRIPT="/scratch/sfeng56/draft/mr_analysis_ba9.R"

# Load R module (adjust based on your system)
module load gcc/12.1.0
module load R/4.5.1

# Get chromosome number from SLURM array task ID
CHR=${SLURM_ARRAY_TASK_ID}

# Print job information
echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Chromosome: ${CHR}"
echo "Start Time: $(date)"
echo "=========================================="
echo ""

# Check if input file exists
INPUT_FILE="${BASE_DIR}/ldByGene${CHR}/BA9_clumped_w50kb_p0.001.txt"
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

echo "Input file found: $INPUT_FILE"
echo ""

# Run the R script
echo "Starting MR analysis for chromosome ${CHR}..."
Rscript ${R_SCRIPT} ${CHR} ${BASE_DIR}

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "SUCCESS: Chromosome ${CHR} completed"
    echo "End Time: $(date)"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "ERROR: Chromosome ${CHR} failed"
    echo "End Time: $(date)"
    echo "=========================================="
    exit 1
fi