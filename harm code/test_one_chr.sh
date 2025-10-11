#!/bin/bash
#SBATCH --job-name=test_one_chr
#SBATCH --time=6:00:00
#SBATCH --error=test_chr1.err
#SBATCH --output=test_chr1.out
#SBATCH --mem=128G
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load R/4.5.1

Rscript process_eqtl.R

