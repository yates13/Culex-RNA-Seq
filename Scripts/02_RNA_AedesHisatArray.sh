#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=HisatArray
#SBATCH --output=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog/%x.%j.out
#SBATCH --error=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog/%x.%j.err
#SBATCH --time=18:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shawn.yates@colostate.edu

# =============================================================================
# 02_RNA_AedesHisatArray.sh   (name is legacy — actually runs on Culex tarsalis)
#
# Purpose:
#   SLURM array job. Each array task aligns one paired-end RNA-seq sample
#   (trimmed FASTQs from step 01b) to the CtarK1 HISAT2 index, producing a
#   SAM file per sample.
#
# Input:
#   ${TRIM_DIR}/CulexTrimmed.txt   - list of R1 filenames, one per line,
#                                     line number (0-indexed) = array task ID
#   ${TRIM_DIR}/<sample>_R1_...fastq.gz  and matching _R2_ file
#
# Output:
#   ${BASE_DIR}/Data/02_RNA_AedesHisatArray/CxtHisat/<sample>.sam
#
# Setup (run once before first submission):
#   cd <trimmed dir produced by 01b>
#   ls | grep -v / | grep "_R1_001.fastq.gz" | grep "trimmed" > CulexTrimmed.txt
#   wc -l CulexTrimmed.txt
#
# Usage:
#   sbatch --array=0-<N> 02_RNA_AedesHisatArray.sh
#   where N = (number of lines in CulexTrimmed.txt) - 1
#   example (6 samples): sbatch --array=0-5 02_RNA_AedesHisatArray.sh
#
# Notes:
#   - HISAT_INDEX below points directly at Index_Genome/ (NOT a "hisatIndex"
#     subfolder) — confirmed via `ls` that the .ht2 index files live there
#     directly. If you rebuild the index in a different location later,
#     update this path to match.
#   - conda activate error-check is enabled (uncommented) below so the job
#     fails loudly instead of silently continuing with no environment if
#     activation breaks.
# =============================================================================

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

# --- Sanity checks before doing any work -------------------------------------
filename="${TRIM_DIR}/CulexTrimmed.txt"

if [ ! -d "${TRIM_DIR}" ]; then
    echo "ERROR: TRIM_DIR does not exist: ${TRIM_DIR}"
    exit 1
fi

if [ ! -f "${filename}" ]; then
    echo "ERROR: Sample list not found: ${filename}"
    echo "  Generate it first (see header comment) before submitting this array job."
    exit 1
fi

# Confirm the HISAT2 index actually exists at this prefix before wasting
# array-task compute time on a broken -x path.
if [ ! -f "${HISAT_INDEX}.1.ht2" ] && [ ! -f "${HISAT_INDEX}.1.ht2l" ]; then
    echo "ERROR: HISAT2 index files not found at prefix: ${HISAT_INDEX}"
    exit 1
fi

# --- Create output directory --------------------------------------------------
mkdir -vp "${OUT_DIR}"

#filename="${TRIM_DIR}/CulexTrimmed.txt" # $1
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
