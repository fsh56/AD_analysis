#!/bin/bash
#SBATCH --job-name=mr_analysis_sep
#SBATCH --output=/scratch/sfeng56/ad_analysis/logs/mr_analysis_separate_%j.out
#SBATCH --error=/scratch/sfeng56/ad_analysis/logs/mr_analysis_separate_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G

INPUT_DIR="/scratch/sfeng56/ad_analysis/clumped_data"
OUTPUT_DIR="/scratch/sfeng56/data/mr_results2"
R_SCRIPT="/scratch/sfeng56/draft/mr_analysis_separate.R"

mkdir -p "$OUTPUT_DIR"

module load gcc/12.1.0
module load R/4.5.1

for INPUT_FILE in ${INPUT_DIR}/*.txt.gz; do
    
    GENE_ID=$(basename "$INPUT_FILE" | grep -oP 'ENSG\d+\.\d+')
    OUTPUT_PREFIX="${OUTPUT_DIR}/${GENE_ID}"
    
    echo "Processing: $GENE_ID"
    Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_PREFIX"
    
done

cat ${OUTPUT_DIR}/*_IVW_summary.csv | head -1 > ${OUTPUT_DIR}/all_methods_combined.csv
for file in ${OUTPUT_DIR}/*_summary.csv; do
    tail -n +2 "$file" >> ${OUTPUT_DIR}/all_methods_combined.csv
done

echo "Done! Results in: $OUTPUT_DIR"
