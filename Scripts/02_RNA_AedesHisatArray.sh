#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=HisatArray
#SBATCH --output=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog/%x.%j.out
#SBATCH --error=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog/%x.%j.err
#SBATCH --time=2:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shawn.yates@colostate.edu

#Run this first to generate a list of files (from inside the trimmed dir produced by 01b):
#ls | grep -v / | grep ".fastq.gz" | grep "R1" | grep "trimmed" > CulexTrimmed.txt
#wc -l CulexTrimmed.txt
#sbatch --array=0-<N> <script.sh>
# example for this subset      sbatch --array=0-5 02_RNA_AedesHisatArray.sh
#Where N is the number of lines in your AedesTrimmed.txt file minus 1. You can get that easily using the following command
# to run this script use sbatch --array=0-5 02_RNA_AedesHisatArray.sh

echo "[$0] $SLURM_JOB_NAME $@" # log the command line
module purge

# --- Load conda (avoids @ symbol issues) ---
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3"; exit 1; }
#conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/hisat2_env || { echo "ERROR: conda activate hisat2_env failed"; exit 1; }
conda activate ${SCRATCH}/conda_envs/hisat2_env
date # timestamp

# --- Path setup: Data/<script_name>/ output, reading from prior stage's Data folder ---
BASE_DIR="/scratch/alpine/c832500103@colostate.edu/RNA_Seq"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}" .sh)"
TRIM_DIR="${BASE_DIR}/Data/slurm_script/trimmed"
OUT_DIR="${BASE_DIR}/Data/${SCRIPT_NAME}/CxtHisat"
HISAT_INDEX="${BASE_DIR}/Genome/Index_Genome/hisatIndex/CtarK1_index"


mkdir -vp "${OUT_DIR}"

filename="${TRIM_DIR}/CulexTrimmed.txt" # $1
linenum=0
while read -r line
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then
      TWO=${line/_R1_/_R2_}
      OUT=${TWO%R2_001.fastq.gz} #Remove end.
      echo ${TWO}
      echo $OUT
      hisat2 --phred33 --rna-strandness RF -p 16 \
        -x "${HISAT_INDEX}" \
        -1 "${TRIM_DIR}/${line}" -2 "${TRIM_DIR}/${TWO}" \
        -S "${OUT_DIR}/${OUT}.sam"
    fi
    linenum=$((linenum + 1))
done < "$filename"
