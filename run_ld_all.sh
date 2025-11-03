#!/bin/bash
#SBATCH --job-name=ld_pergene
#SBATCH --output=/scratch/sfeng56/logs/ld_ba9_chr%a.out
#SBATCH --error=/scratch/sfeng56/logs/ld_ba9_chr%a.err
#SBATCH --array=1-22
#SBATCH --time=4:00:00
#SBATCH --mem=8G

CHR=${SLURM_ARRAY_TASK_ID}

# Load required modules
module load gcc/12.1.0
module load R/4.5.1
module load lapack/3.11.0
module load atlas/3.10.3
module load plink/1.9

# Set variables (使用 CHR 变量)
INPUT_FILE="/gpfs/data/gao-lab/people/Sihao/data/harmonized_data_Brain_Frontal_Cortex_BA9/eqtl_amyloid_Brain_Frontal_Cortex_BA9_${CHR}.txt.gz"
OUTPUT_DIR="/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/ldByGene${CHR}"
R_SCRIPT="/gpfs/data/gao-lab/people/Sihao/correct_code/ld_code/ld_per_gene.R"
LOG_DIR="/scratch/sfeng56/logs"

# Create directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Fixed parameters
WINDOW_KB=50
PVAL=0.0001

echo "=================================================="
echo "LD Clumping Job Started"
echo "Chromosome: chr${CHR}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=================================================="

# 检查输入文件
if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: Input file not found: ${INPUT_FILE}"
    exit 1
fi

# 运行分析
echo "Processing chr${CHR}..."
Rscript ${R_SCRIPT} ${INPUT_FILE} ${OUTPUT_DIR} ${WINDOW_KB} ${PVAL}

if [ $? -eq 0 ]; then
    echo "SUCCESS: Completed chr${CHR}"
else
    echo "ERROR: Failed chr${CHR}"
    exit 1
fi

echo "End time: $(date)"
echo "=================================================="