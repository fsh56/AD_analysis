#!/bin/bash
#SBATCH --job-name=mr_amyloid
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/amyloid_amydgala/logs/mr_%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/amyloid_amydgala/logs/mr_%j.err
#SBATCH --time=3:00:00
#SBATCH --partition=tier2q
#SBATCH --mem=16G

# Directories and files
INPUT_DIR="/scratch/sfeng56/ad_analysis/clump1020/Brain_Amygdala"
OUTPUT_BASE="/gpfs/data/gao-lab/people/Sihao/amyloid_amydgala/mr_results0.0005"
R_SCRIPT="/scratch/sfeng56/draft/mr3.R"

# Create output directories
mkdir -p "$OUTPUT_BASE"
mkdir -p "/gpfs/data/gao-lab/people/Sihao/amyloid_amydgala/logs"

# Load modules
module load gcc/12.1.0
module load R/4.5.1

echo "======================================================================"
echo "Starting MR analysis for amyloid amydgala 1 data"
echo "Date: $(date)"
echo "======================================================================"
echo ""

# Define the three input files
declare -a INPUT_FILES=(
  "${INPUT_DIR}/clumped_w50kb_p0.0005.txt"
  "${INPUT_DIR}/clumped_w20kb_p0.0005.txt "
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

# Combine IVW results
echo "Combining IVW results..."
IVW_FILES=("${OUTPUT_BASE}"/amyloid_IVW_p*.csv)
if [ ${#IVW_FILES[@]} -gt 0 ] && [ -f "${IVW_FILES[0]}" ]; then
  head -1 "${IVW_FILES[0]}" > "${OUTPUT_BASE}/all_IVW_combined.csv"
  for file in "${IVW_FILES[@]}"; do
    tail -n +2 "$file" >> "${OUTPUT_BASE}/all_IVW_combined.csv"
  done
  echo "✓ IVW combined file created: all_IVW_combined.csv"
  echo "  Total rows: $(wc -l < "${OUTPUT_BASE}/all_IVW_combined.csv")"
else
  echo "✗ No IVW files found"
fi

echo ""

# Combine Egger results
echo "Combining Egger results..."
EGGER_FILES=("${OUTPUT_BASE}"/amyloid_Egger_p*.csv)
if [ ${#EGGER_FILES[@]} -gt 0 ] && [ -f "${EGGER_FILES[0]}" ]; then
  head -1 "${EGGER_FILES[0]}" > "${OUTPUT_BASE}/all_Egger_combined.csv"
  for file in "${EGGER_FILES[@]}"; do
    tail -n +2 "$file" >> "${OUTPUT_BASE}/all_Egger_combined.csv"
  done
  echo "✓ Egger combined file created: all_Egger_combined.csv"
  echo "  Total rows: $(wc -l < "${OUTPUT_BASE}/all_Egger_combined.csv")"
else
  echo "✗ No Egger files found"
fi

echo ""

# Combine Weighted Median results (新增)
echo "Combining Weighted Median results..."
WM_FILES=("${OUTPUT_BASE}"/amyloid_WeightedMedian_p*.csv)
if [ ${#WM_FILES[@]} -gt 0 ] && [ -f "${WM_FILES[0]}" ]; then
  head -1 "${WM_FILES[0]}" > "${OUTPUT_BASE}/all_WeightedMedian_combined.csv"
  for file in "${WM_FILES[@]}"; do
    tail -n +2 "$file" >> "${OUTPUT_BASE}/all_WeightedMedian_combined.csv"
  done
  echo "✓ Weighted Median combined file created: all_WeightedMedian_combined.csv"
  echo "  Total rows: $(wc -l < "${OUTPUT_BASE}/all_WeightedMedian_combined.csv")"
else
  echo "✗ No Weighted Median files found"
fi

echo ""

# Create master combined file (更新)
echo "Creating master combined file..."
MASTER_FILE="${OUTPUT_BASE}/all_methods_combined.csv"
TEMP_CREATED=false

if [ -f "${OUTPUT_BASE}/all_IVW_combined.csv" ]; then
  cat "${OUTPUT_BASE}/all_IVW_combined.csv" > "$MASTER_FILE"
  TEMP_CREATED=true
  echo "  Added IVW results"
fi

if [ -f "${OUTPUT_BASE}/all_Egger_combined.csv" ]; then
  if [ "$TEMP_CREATED" = true ]; then
    tail -n +2 "${OUTPUT_BASE}/all_Egger_combined.csv" >> "$MASTER_FILE"
  else
    cat "${OUTPUT_BASE}/all_Egger_combined.csv" > "$MASTER_FILE"
    TEMP_CREATED=true
  fi
  echo "  Added Egger results"
fi

if [ -f "${OUTPUT_BASE}/all_WeightedMedian_combined.csv" ]; then
  if [ "$TEMP_CREATED" = true ]; then
    tail -n +2 "${OUTPUT_BASE}/all_WeightedMedian_combined.csv" >> "$MASTER_FILE"
  else
    cat "${OUTPUT_BASE}/all_WeightedMedian_combined.csv" > "$MASTER_FILE"
    TEMP_CREATED=true
  fi
  echo "  Added Weighted Median results"
fi

if [ "$TEMP_CREATED" = true ]; then
  echo "✓ Master combined file created: all_methods_combined.csv"
  echo "  Total rows: $(wc -l < "$MASTER_FILE")"
else
  echo "✗ Cannot create master file - no result files found"
fi

echo ""
echo "======================================================================"
echo "Analysis complete!"
echo "Results directory: $OUTPUT_BASE"
echo "End time: $(date)"
echo "======================================================================"
echo ""
echo "Output files:"
ls -lh "${OUTPUT_BASE}"/*.csv 2>/dev/null || echo "No CSV files found"

