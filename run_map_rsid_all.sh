#!/bin/bash
#SBATCH --job-name=map_rsid
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/map_rsid_%A_%a.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/map_rsid_%A_%a.err
#SBATCH --array=2-13
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --partition=tier1q

# make log directory
mkdir -p /gpfs/data/gao-lab/people/Sihao/ad_analysis/logs

# load modules
module load gcc/12.1.0
module load R/4.5.1

# define tissues
tissues=(
  "Brain_Amygdala"
  "Brain_Anterior_cingulate_cortex_BA24"
  "Brain_Caudate_basal_ganglia"
  "Brain_Cerebellar_Hemisphere"
  "Brain_Cerebellum"
  "Brain_Cortex"
  "Brain_Frontal_Cortex_BA9"
  "Brain_Hippocampus"
  "Brain_Hypothalamus"
  "Brain_Nucleus_accumbens_basal_ganglia"
  "Brain_Putamen_basal_ganglia"
  "Brain_Spinal_cord_cervical_c-1"
  "Brain_Substantia_nigra"
)

# 获取当前任务的组织
tissue=${tissues[$SLURM_ARRAY_TASK_ID-1]}

echo "processing tissue: $tissue"

Rscript map_rsid_all.R $tissue