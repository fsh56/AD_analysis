#!/bin/bash
#SBATCH --job-name=seso
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/logs/mr_seso%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/logs/mr_seso%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=32G

# ===== 添加这部分：设置临时目录 =====
export TMPDIR=/scratch/sfeng56/tmp
export TEMP=$TMPDIR
export TMP=$TMPDIR
mkdir -p $TMPDIR
echo "Using temporary directory: $TMPDIR"
# =====================================

# Directories and files
INPUT_DIR="/scratch/sfeng56/ad_analysis/amyloid/BA9"
OUTPUT_BASE="/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/mr_results"
R_SCRIPT="/scratch/sfeng56/draft/mr_seso.R"

# Create output directories
mkdir -p "$OUTPUT_BASE"
mkdir -p "/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/logs/"

# Load modules
module load gcc/12.1.0
module load R/4.5.1

echo "======================================================================"
echo "Starting MR SESO analysis for amyloid ba9 data"
echo "Date: $(date)"
echo "======================================================================"
echo ""

# Define the input files
declare -a INPUT_FILES=(
  "${INPUT_DIR}/clumped_amyloid_BA9_1_w20kb_p0.0001.txt"
  "${INPUT_DIR}/clumped_amyloid_BA9_1_w20kb_p0.001.txt"
)

# Process each file
for INPUT_FILE in "${INPUT_FILES[@]}"; do
  if [ ! -f "$INPUT_FILE" ]; then
    echo "WARNING: File not found: $INPUT_FILE"
    echo "Skipping..."
    echo ""
    continue
  fi
  
  BASENAME=$(basename "$INPUT_FILE" .txt)
  
  echo "----------------------------------------------------------------------"
  echo "Processing: $BASENAME"
  echo "Input file: $INPUT_FILE"
  echo "Start time: $(date)"
  echo "----------------------------------------------------------------------"
  
  # Run R script
  Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_BASE"
  
  EXIT_CODE=$?
  
  if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Successfully completed"
  else
    echo "✗ Failed with exit code: $EXIT_CODE"
  fi
  
  echo "End time: $(date)"
  echo ""
done

echo "======================================================================"
echo "Combining all results across p-value thresholds..."
echo "======================================================================"
echo ""

# Combine seso results
echo "Combining seso results..."
SESO_FILES=("${OUTPUT_BASE}"/amyloid_seso_p*.csv)

if [ ${#SESO_FILES[@]} -gt 0 ] && [ -f "${SESO_FILES[0]}" ]; then
  head -1 "${SESO_FILES[0]}" > "${OUTPUT_BASE}/all_seso_combined.csv"
  for file in "${SESO_FILES[@]}"; do
    tail -n +2 "$file" >> "${OUTPUT_BASE}/all_seso_combined.csv"
  done
  echo "✓ SESO combined file created: all_seso_combined.csv"
  echo "  Total rows: $(wc -l < "${OUTPUT_BASE}/all_seso_combined.csv")"
else
  echo "✗ No SESO files found"
fi

echo ""
echo "======================================================================"
echo "MR SESO analysis completed"
echo "End time: $(date)"
echo "======================================================================"