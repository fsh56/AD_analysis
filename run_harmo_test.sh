#!/bin/bash
#SBATCH --job-name=harmo_test
#SBATCH --output=/scratch/sfeng56/draft/code1027/logs/harmo_test%j.out
#SBATCH --error=/scratch/sfeng56/draft/code1027/logs/harmo_test%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G

# Load modules
module load gcc/12.1.0
module load R/4.5.1

# Directories and files
R_SCRIPT="/scratch/sfeng56/draft/code1027/harmo_test.R"

Rscript "$R_SCRIPT"