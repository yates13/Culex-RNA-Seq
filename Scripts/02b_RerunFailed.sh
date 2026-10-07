#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=HisatRerun
#SBATCH --output=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog/%x.%A-%a.log
#SBATCH --error=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog/%x.%A-%a.err
#SBATCH --time=20:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shawn.yates@colostate.edu

echo "[$0] $SLURM_JOB_NAME $@"
module purge

SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found"; exit 1; }
conda activate ${SCRATCH}/conda_envs/hisat2_env || { echo "ERROR: conda activate failed"; exit 1; }
date

BASE_DIR="/scratch/alpine/c832500103@colostate.edu/RNA_Seq"
TRIM_DIR="${BASE_DIR}/Data/slurm_script/trimmed"
OUT_DIR="${BASE_DIR}/Data/02_RNA_AedesHisatArray/CxtHisat"
HISAT_INDEX="${BASE_DIR}/Genome/Index_Genome/hisatIndex/CtarK1_index"

declare -a R1_FILES=(
#  "trimmed.BF-RNA-D1-R1-Cxt_R1_001.fastq.gz"
  "trimmed.ZH-RNA-D3-R5-Cxt_R1_001.fastq.gz"
)

line="${R1_FILES[$SLURM_ARRAY_TASK_ID]}"
TWO=${line/_R1_/_R2_}
OUT=${TWO%R2_001.fastq.gz}
#OUT=${line%R1_001.fastq.gz}

echo "R1: ${line}"
echo "R2: ${TWO}"
echo "Output prefix: ${OUT}"

#hisat2 --phred33 --rna-strandness RF -p 16 \
#    -x "${HISAT_INDEX}" \
#    -1 "${TRIM_DIR}/${line}" -2 "${TRIM_DIR}/${TWO}" \
#    | samtools sort -@ 4 -o "${OUT_DIR}/${OUT}.bam" -

#samtools index "${OUT_DIR}/${OUT}.bam"
hisat2 --phred33 --rna-strandness RF -p 16 \
    -x "${HISAT_INDEX}" \
    -1 "${TRIM_DIR}/${line}" -2 "${TRIM_DIR}/${TWO}" \
    -S "${OUT_DIR}/${OUT}.sam"
