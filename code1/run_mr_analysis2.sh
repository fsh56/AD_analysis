#!/bin/bash
#SBATCH --job-name=mr_analysis_batch
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/mr_analysis_batch_%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/mr_analysis_batch_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G

INPUT_DIR="/gpfs/data/gao-lab/people/Sihao/ad_analysis/clumped_data"
OUTPUT_DIR="/gpfs/data/gao-lab/people/Sihao/data/mr_results2"
R_SCRIPT="/gpfs/data/gao-lab/people/Sihao/draft/mr_analysis2.R"

mkdir -p "$OUTPUT_DIR"

module load gcc/12.1.0
module load R/4.5.1

for INPUT_FILE in ${INPUT_DIR}/*.txt.gz; do
    
    GENE_ID=$(basename "$INPUT_FILE" | grep -oP 'ENSG\d+\.\d+')
    OUTPUT_PREFIX="${OUTPUT_DIR}/${GENE_ID}"
    
    echo "Processing: $GENE_ID"
    Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_PREFIX"
    
done

echo "Done! Results in: $OUTPUT_DIR"