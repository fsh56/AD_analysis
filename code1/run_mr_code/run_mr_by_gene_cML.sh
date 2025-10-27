#!/bin/bash
#SBATCH --job-name=mr_cML
#SBATCH --output=/scratch/sfeng56/ad_analysis/logs/mr_by_gene_cML_%j.out
#SBATCH --error=/scratch/sfeng56/ad_analysis/logs/mr_by_gene_cML_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

INPUT_DIR="/scratch/sfeng56/ad_analysis/clumped_data"
OUTPUT_BASE="/scratch/sfeng56/data/mr_results_by_gene"
R_SCRIPT="/scratch/sfeng56/draft/mr_by_gene_cML.R"

mkdir -p "$OUTPUT_BASE"
mkdir -p "/scratch/sfeng56/ad_analysis/logs"

module load gcc/12.1.0
module load R/4.5.1

# Process each clumped file
for INPUT_FILE in ${INPUT_DIR}/*.txt.gz; do
    
    BASENAME=$(basename "$INPUT_FILE" .txt.gz)
    # Remove the _20kb_p0.001.txt.gz suffix to get a cleaner name
    TISSUE_NAME=$(echo "$BASENAME" | sed 's/_20kb_p0\.001//g' | sed 's/clumped_eqtl_//g')
    
    OUTPUT_DIR="${OUTPUT_BASE}/${TISSUE_NAME}"
    mkdir -p "$OUTPUT_DIR"
    
    echo "========================================"
    echo "Processing file: $BASENAME"
    echo "Output directory: $OUTPUT_DIR"
    echo "========================================"
    
    Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_DIR"
    
    echo ""
done

echo "combining all cML results..."

# Check if any cML results exist
CML_COUNT=$(find "$OUTPUT_BASE" -name "*_cML.csv" | wc -l)

if [ $CML_COUNT -eq 0 ]; then
    echo "Warning: No cML results found!"
    exit 1
else
    echo "Found $CML_COUNT cML result files"
    
    # Combine all cML results
    find "$OUTPUT_BASE" -name "*_cML.csv" -exec head -1 {} \; -quit > "${OUTPUT_BASE}/all_cML_combined.csv"
    find "$OUTPUT_BASE" -name "*_cML.csv" -exec tail -n +2 {} \; >> "${OUTPUT_BASE}/all_cML_combined.csv"
    
    # Count total genes
    TOTAL_GENES=$(tail -n +2 "${OUTPUT_BASE}/all_cML_combined.csv" | wc -l)
    echo "Total genes analyzed: $TOTAL_GENES"
fi

echo ""
echo "Done! Results saved in: $OUTPUT_BASE"
echo "Summary file: all_cML_combined.csv"