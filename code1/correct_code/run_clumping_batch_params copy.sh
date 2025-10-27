#!/bin/bash
#SBATCH --job-name=ld_clump_batch
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/ld_clump_batch_%A_%a.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/ld_clump_batch_%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --partition=tier2q
#SBATCH --array=1-176%10

# Create log directory if it doesn't exist
mkdir -p /gpfs/data/gao-lab/people/Sihao/ad_analysis/logs

# Load required modules
module load gcc/12.1.0
module load R/4.5.1
module load lapack/3.11.0
module load atlas/3.10.3
module load plink/1.9

# Define paths
INPUT_DIR="/gpfs/data/gao-lab/people/Sihao/ad_analysis/harmonized_data"
OUTPUT_DIR="/gpfs/data/gao-lab/people/Sihao/ad_analysis/clumped_data"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

CLUMP_KB=20 
CLUMP_R2=0.1
CLUMP_P=0.001

FILES=(${INPUT_DIR}/eqtl_*.txt.gz)
INPUT_FILE="${FILES[$SLURM_ARRAY_TASK_ID-1]}"
BASENAME=$(basename ${INPUT_FILE})
OUTPUT_FILE="${OUTPUT_DIR}/clumped_${BASENAME}_${CLUMP_KB}kb_p${CLUMP_P}.txt.gz"

# Print job information
echo "Job ID: $SLURM_ARRAY_JOB_ID"
echo "Input file: $INPUT_FILE"
echo "clump_kb: $CLUMP_KB"
echo "clump_r2: $CLUMP_R2"
echo "clump_p: $CLUMP_P"
echo ""

if [ ! -f "$INPUT_FILE" ]; then
echo "Error: Input file not found: $INPUT_FILE"
exit 1
fi
if [ -f "$OUTPUT_FILE" ]; then
echo "Output file already exists, skipping: $OUTPUT_FILE"
exit 0
fi

# Run R script
echo "Starting clumping for: $(basename $INPUT_FILE)"
Rscript /gpfs/data/gao-lab/people/Sihao/draft/clumping_params.R ${INPUT_FILE} ${OUTPUT_FILE} ${CLUMP_KB} ${CLUMP_R2} ${CLUMP_P}

# Check if R script succeeded
if [ $? -eq 0 ]; then
echo "Successfully completed: $(basename $INPUT_FILE)"
else
  echo "Error processing: $(basename $INPUT_FILE)"
exit 1
fi

echo "All output files are saved to $OUTPUT_FILE"
echo "Clumping done!"
