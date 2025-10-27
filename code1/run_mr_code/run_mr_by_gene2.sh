#!/bin/bash
#SBATCH --job-name=mr_by_gene
#SBATCH --output=/scratch/sfeng56/ad_analysis/logs/mr_by_gene2_%j.out
#SBATCH --error=/scratch/sfeng56/ad_analysis/logs/mr_by_gene2_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

INPUT_DIR="/scratch/sfeng56/ad_analysis/clumped_data"
OUTPUT_BASE="/scratch/sfeng56/data/mr_results_by_gene"
R_SCRIPT="/scratch/sfeng56/draft/mr_by_gene2.R"

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
    
    echo "processing file: $BASENAME"
    echo "output directory: $OUTPUT_DIR"
    Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_DIR"
    
    echo ""
done

echo "combining all results..."
echo ""

# Combine all IVW results
echo "Combining IVW results..."
find "$OUTPUT_BASE" -name "*_IVW.csv" -exec head -1 {} \; -quit > "${OUTPUT_BASE}/all_IVW_combined.csv"
find "$OUTPUT_BASE" -name "*_IVW.csv" -exec tail -n +2 {} \; >> "${OUTPUT_BASE}/all_IVW_combined.csv"

# Combine all Egger results
echo "Combining Egger results..."
find "$OUTPUT_BASE" -name "*_Egger.csv" -exec head -1 {} \; -quit > "${OUTPUT_BASE}/all_Egger_combined.csv"
find "$OUTPUT_BASE" -name "*_Egger.csv" -exec tail -n +2 {} \; >> "${OUTPUT_BASE}/all_Egger_combined.csv"

# Create a master combined file with all methods
echo "Creating master combined file..."
cat "${OUTPUT_BASE}/all_IVW_combined.csv" > "${OUTPUT_BASE}/all_methods_combined.csv"
tail -n +2 "${OUTPUT_BASE}/all_Egger_combined.csv" >> "${OUTPUT_BASE}/all_methods_combined.csv"

echo "Done! Results saved in: $OUTPUT_BASE"
echo "Summary files:"
echo "  - all_IVW_combined.csv"
echo "  - all_Egger_combined.csv"
echo "  - all_methods_combined.csv"