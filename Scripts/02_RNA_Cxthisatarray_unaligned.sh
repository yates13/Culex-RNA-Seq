#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=HisatCxtRVFV
#SBATCH --output=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog/%x.%A_%a.out
#SBATCH --error=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog/%x.%A_%a.err
#SBATCH --time=18:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shawn.yates@colostate.edu

# =============================================================================
# 02_RNA_CxtRVFV_HisatArray.sh
#
# Purpose:
#   SLURM array job. Each task processes one trimmed paired-end RNA-seq sample
#   in two sequential HISAT2 stages and keeps THREE datasets:
#
#   Stage 1: trimmed reads  -> Cxt GENOME index
#       OUTPUT 1  Cxt-aligned reads (SAM)
#                 Cxt-unaligned pairs (FASTQ.gz)  -> input to stage 2
#   Stage 2: Cxt-unaligned pairs -> RVFV index
#       OUTPUT 2  RVFV-aligned reads (SAM)
#       OUTPUT 3  pairs that aligned to NEITHER Cxt nor RVFV (FASTQ.gz)
#                 -> ssRNA+ re-alignment / BLAST vs insect-specific viruses
#
# Input:
#   ${TRIM_DIR}/CulexTrimmed.txt  - list of R1 filenames, one per line
#                                   (line number, 0-indexed = array task ID)
#
# Output (under ${BASE_DIR}/Data/<script_name>/):
#   1_CxtAligned/<sample>.sam
#   1_CxtUnaligned/<sample>unaligned_R1/R2.fastq.gz     (intermediate)
#   2_RVFVAligned/<sample>.sam
#   3_NeitherAligned/<sample>unaligned_R1/R2.fastq.gz   (final "neither" set)
#   Summaries/<sample>cxt_summary.txt, <sample>rvfv_summary.txt
#
# Usage:
#   sbatch --array=0-<N> 02_RNA_CxtRVFV_HisatArray.sh
#   N = (lines in CulexTrimmed.txt) - 1
#
# Notes:
#   - The RVFV index must contain all three genome segments (L, M, S) in one
#     FASTA before running hisat2-build. See HISAT_RVFV below.
#   - Stage 2 only sees reads that did not align concordantly to Cxt.
# =============================================================================

echo "[$0] $SLURM_JOB_NAME $@"
module purge

# --- Load conda ---
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3"; exit 1; }
conda activate ${SCRATCH}/conda_envs/hisat2_env || { echo "ERROR: conda activate hisat2_env failed"; exit 1; }
date

# --- Paths ---
BASE_DIR="/scratch/alpine/c832500103@colostate.edu/RNA_Seq"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}" .sh)"
TRIM_DIR="${BASE_DIR}/Data/slurm_script/trimmed"
RUN_DIR="${BASE_DIR}/Data/${SCRIPT_NAME}"

CXT_SAM_DIR="${RUN_DIR}/1_CxtAligned"
CXT_UNAL_DIR="${RUN_DIR}/1_CxtUnaligned"
RVFV_SAM_DIR="${RUN_DIR}/2_RVFVAligned"
NEITHER_DIR="${RUN_DIR}/3_NeitherAligned"
SUM_DIR="${RUN_DIR}/Summaries"

# Cxt GENOME index prefix (same path as the previous script)
HISAT_CXT="${BASE_DIR}/Genome/Index_Genome/hisatIndex/CtarK1_index"
# RVFV index prefix -- CONFIRM this path once the RVFV index is built/uploaded
HISAT_RVFV="${BASE_DIR}/Genome/Index_RVFV/RVFV_index"

THREADS=${SLURM_CPUS_PER_TASK:-16}

# --- Sanity checks ---
filename="${TRIM_DIR}/CulexTrimmed.txt"

[ -d "${TRIM_DIR}" ] || { echo "ERROR: TRIM_DIR does not exist: ${TRIM_DIR}"; exit 1; }
[ -f "${filename}" ] || { echo "ERROR: Sample list not found: ${filename}"; exit 1; }

for IDX in "${HISAT_CXT}" "${HISAT_RVFV}"; do
    if [ ! -f "${IDX}.1.ht2" ] && [ ! -f "${IDX}.1.ht2l" ]; then
        echo "ERROR: HISAT2 index files not found at prefix: ${IDX}"
        exit 1
    fi
done

mkdir -vp "${CXT_SAM_DIR}" "${CXT_UNAL_DIR}" "${RVFV_SAM_DIR}" "${NEITHER_DIR}" "${SUM_DIR}"

linenum=0
while read -r line
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then
      TWO=${line/_R1_/_R2_}
      OUT=${TWO%R2_001.fastq.gz}   # sample prefix (keeps trailing "_")
      echo ${TWO}
      echo ${OUT}

      # ---------------- Stage 1: Cxt genome ----------------
      echo "=== Stage 1: Cxt genome alignment ==="
      hisat2 --phred33 --rna-strandness RF -p ${THREADS} \
        -x "${HISAT_CXT}" \
        -1 "${TRIM_DIR}/${line}" -2 "${TRIM_DIR}/${TWO}" \
        --un-conc-gz "${CXT_UNAL_DIR}/${OUT}unaligned_R%.fastq.gz" \
        --summary-file "${SUM_DIR}/${OUT}cxt_summary.txt" \
        -S "${CXT_SAM_DIR}/${OUT}.sam" \
        || { echo "ERROR: Cxt hisat2 failed for ${OUT}"; exit 1; }
      date

      CXT_R1="${CXT_UNAL_DIR}/${OUT}unaligned_R1.fastq.gz"
      CXT_R2="${CXT_UNAL_DIR}/${OUT}unaligned_R2.fastq.gz"
      if [ ! -s "${CXT_R1}" ] || [ ! -s "${CXT_R2}" ]; then
          echo "ERROR: no Cxt-unaligned FASTQs for ${OUT} (every pair aligned to Cxt?)"
          exit 1
      fi

      # ---------------- Stage 2: RVFV ----------------
      echo "=== Stage 2: RVFV alignment of Cxt-unaligned reads ==="
      hisat2 --phred33 --rna-strandness RF -p ${THREADS} \
        --no-spliced-alignment \
        -x "${HISAT_RVFV}" \
        -1 "${CXT_R1}" -2 "${CXT_R2}" \
        --un-conc-gz "${NEITHER_DIR}/${OUT}unaligned_R%.fastq.gz" \
        --summary-file "${SUM_DIR}/${OUT}rvfv_summary.txt" \
        -S "${RVFV_SAM_DIR}/${OUT}.sam" \
        || { echo "ERROR: RVFV hisat2 failed for ${OUT}"; exit 1; }
      date
    fi
    linenum=$((linenum + 1))
done < "$filename"
