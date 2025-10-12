#!/bin/bash
#SBATCH --job-name=gwas_map_rsid
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/gwas_map_rsid_%A_%a.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/gwas_map_rsid_%A_%a.err
#SBATCH --array=1-8
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=tier2q

# make log directory
mkdir -p /gpfs/data/gao-lab/people/Sihao/ad_analysis/logs

# load modules
module load gcc/12.1.0
module load R/4.5.1

# Define GWAS files array
FILES=(
  "amyloid.assoc.linear.gz"
  "dcfdx_ad.assoc.logistic.gz"
  "gpath.assoc.linear.gz"
  "cogdx_ad.assoc.logistic.gz"
  "dlbany.assoc.logistic.gz"
  "hspath_typ.assoc.logistic.gz"
  "tdp_st4_binary.assoc.logistic.gz"
  "tangles.assoc.linear.gz"
)

FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}
Rscript map_rsid_gwas.R $FILE
echo "Completed processing: $FILE"