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
# One-time setup: build the RVFV HISAT2 index (run once before first submission)
#   Reference: RefSeq segments NC_014397 (L), NC_014396 (M), NC_014395 (S).
#   RVFV has three segments, so all three go into ONE fasta.
#
#   mkdir -p ${BASE_DIR}/Genome/Index_RVFV
#   cd ${BASE_DIR}/Genome/Index_RVFV
#
#   # 1. Download all three segments as a single fasta (works from login node)
#   curl -o RVFV_all.fasta \
#     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_014397,NC_014396,NC_014395&rettype=fasta&retmode=text"
#
#   # 2. Check: 3 sequences, approx L 6.4 kb, M 3.9 kb, S 1.7 kb
#   grep ">" RVFV_all.fasta
#   awk '/^>/{if(n)print n,l; n=$1; l=0; next}{l+=length($0)}END{print n,l}' RVFV_all.fasta
#
#   # 3. Build the index (hisat2_env active). Takes seconds; ~11 kb genome.
#   hisat2-build -p 4 RVFV_all.fasta RVFV_index
#   ls RVFV_index.*.ht2      # expect 8 files, .1.ht2 through .8.ht2
#
#   Result: index prefix = ${BASE_DIR}/Genome/Index_RVFV/RVFV_index (HISAT_RVFV below)
#
# Inputs required BEFORE running (the script checks 1-4 and exits if missing):
#   1. Trimmed paired-end FASTQs in ${TRIM_DIR}
#      (${BASE_DIR}/Data/slurm_script/trimmed). Each sample needs an R1 and R2:
#      <sample>_R1_001.fastq.gz and <sample>_R2_001.fastq.gz, with "trimmed"
#      in the name. Example: trimmed.BF-RNA-D1-R1-Cxt_R1_001.fastq.gz
#   2. Sample list ${TRIM_DIR}/CulexTrimmed.txt (see "How the array is set up")
#   3. Cxt GENOME HISAT2 index: ${BASE_DIR}/Genome/Index_Genome/hisatIndex/
#      CtarK1_index.*.ht2 (8 files; the prefix is CtarK1_index)
#   4. RVFV HISAT2 index: ${BASE_DIR}/Genome/Index_RVFV/RVFV_index.*.ht2
#      (8 files; built with the one-time setup below)
#   5. The hisat2_env conda environment (conda_envs/hisat2_env)
#   6. SLURM log folders slurmlogs/outlog and slurmlogs/errlog (if they are
#      missing, SLURM cannot write the logs and the job can fail silently)
#
# How the array is set up with CulexTrimmed.txt:
#   - CulexTrimmed.txt lists ONE R1 filename per line, one line per sample.
#     Build it in the trimmed folder after trimming is finished:
#       cd ${TRIM_DIR}
#       ls | grep -v / | grep "_R1_001.fastq.gz" | grep "trimmed" > CulexTrimmed.txt
#       wc -l CulexTrimmed.txt
#   - The array task ID is the 0-indexed line number: task 0 runs line 1,
#     task 1 runs line 2, and so on. So the array range is 0 to N-1, where
#     N is the number of lines. For 24 samples: --array=0-23.
#   - Test one sample first (--array=0), check the outputs, then run the
#     rest (--array=1-23).
#   - Regenerate CulexTrimmed.txt whenever trimmed files are added or removed.
#     The list is sorted alphabetically, so adding samples can shift which task
#     ID maps to which sample. Do not mix task IDs from runs that used
#     different versions of the list.
#
# Notes:
#   - Stage 2 only sees reads that did not align concordantly to Cxt, so pairs
#     with one host-aligned mate can pass through to the RVFV/neither sets.
#   - When run through sbatch, SCRIPT_NAME resolves to "slurm_script", so the
#     outputs land in ${BASE_DIR}/Data/slurm_script/, not a folder named after
#     this file.
#   - 1_CxtUnaligned/ is kept after the run (it is the full "not Cxt" set);
#     3_NeitherAligned/ is a subset of it. Delete only after counts are checked.
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
