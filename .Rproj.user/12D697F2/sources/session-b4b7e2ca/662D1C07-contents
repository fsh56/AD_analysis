#!/bin/bash
#SBATCH --job-name=harmonize_gwas_eqtl
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/harmonize_%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/harmonize_%j.err
#SBATCH --time=5:00:00
#SBATCH --mem=32G
#SBATCH --partition=tier1q

# make log directory
mkdir -p /gpfs/data/gao-lab/people/Sihao/ad_analysis/logs

# load modules
module load gcc/12.1.0
module load R/4.5.1

# define args
GWAS_FILE="amyloid.assoc.linear.gz"
TISSUE="Brain_Amygdala"

# run R file
Rscript harm_linear.R ${GWAS_FILE} ${TISSUE}
