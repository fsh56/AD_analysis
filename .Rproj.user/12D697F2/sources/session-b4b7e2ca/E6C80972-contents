#!/bin/bash
#SBATCH --job-name=map_rsid
#SBATCH --output=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/map_rsid_%j.out
#SBATCH --error=/gpfs/data/gao-lab/people/Sihao/ad_analysis/logs/map_rsid_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --partition=tier2q

# make log directory
mkdir -p /gpfs/data/gao-lab/people/Sihao/ad_analysis/logs

# load modules
module load gcc/12.1.0
module load R/4.5.1

Rscript map_rsid.R