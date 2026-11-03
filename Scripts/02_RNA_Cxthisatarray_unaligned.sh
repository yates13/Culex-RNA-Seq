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
# 02_RNA_Cxthisatarray_unaligned.sh
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
#   "Aligned" = at least one mate mapped (concordant, discordant, or one-mate
#   pairs, with both records kept in the SAM). "Unaligned" = BOTH mates
#   unmapped. Each pair lands in exactly one output per stage.
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
#   7. samtools >= 1.10 in hisat2_env (the script checks this and exits if
#      it is missing or too old)
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
#   - HISAT2 writes a raw SAM (including unmapped records); split_sam() then
#     separates it by SAM flags. Aligned SAM = everything except pairs with
#     both mates unmapped (flag 12). Unaligned FASTQs = pairs with both mates
#     unmapped. --un-conc-gz is NOT used, because it also captured discordant
#     and one-mate-aligned pairs, which put host-mapped reads into the
#     "unaligned" sets.
#   - Stage 2 therefore only sees pairs with no Cxt alignment at all. A viral
#     pair with one spuriously host-mapped mate will stay in the Cxt SAM and
#     won't reach the RVFV/neither sets.
#   - The raw SAM is deleted after splitting; budget scratch space for one
#     uncompressed raw SAM per task while it runs.
#   - Requires samtools >= 1.10 in hisat2_env (view -G).
#   - Count check per sample: unique read names in 1_CxtAligned SAM + pairs in
#     1_CxtUnaligned FASTQ = total input pairs. The same holds for stage 2
#     (RVFV SAM + Neither FASTQ = Cxt-unaligned pairs).
#   - When run through sbatch, SCRIPT_NAME resolves to "slurm_script", so the
#     outputs land in ${BASE_DIR}/Data/slurm_script/, not a folder named after
#     this file.
#   - 1_CxtUnaligned/ is kept after the run (full "not Cxt" set);
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

command -v samtools >/dev/null || { echo "ERROR: samtools not in hisat2_env"; exit 1; }
samtools view --help 2>&1 | grep -q -- '-G' || { echo "ERROR: samtools too old (need >= 1.10 for view -G)"; exit 1; }

mkdir -vp "${CXT_SAM_DIR}" "${CXT_UNAL_DIR}" "${RVFV_SAM_DIR}" "${NEITHER_DIR}" "${SUM_DIR}"

# Split a raw HISAT2 SAM into:
#   ALN  = all records except pairs with BOTH mates unmapped (concordant + discordant + one-mate)
#   U1/U2 = FASTQ.gz of pairs with BOTH mates unmapped
split_sam () {
    set -o pipefail
    local RAW="$1" ALN="$2" U1="$3" U2="$4"

    # -G 12 excludes only records with BOTH 0x4 and 0x8 set
    # (needs samtools >= 1.10; the unmapped mate of a one-mate-mapped pair is kept)
    samtools view -h -G 12 "${RAW}" > "${ALN}" || return 1

    samtools view -u -f 12 -F 256 "${RAW}" \
      | samtools sort -n -u -@ "${THREADS}" -T "${RUN_DIR}/${OUT}sorttmp" - \
      | samtools fastq -n -1 "${U1}" -2 "${U2}" -0 /dev/null -s /dev/null - \
      || return 1
}

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
      RAW1="${CXT_SAM_DIR}/${OUT}.raw.sam"
      hisat2 --phred33 --rna-strandness RF -p ${THREADS} \
        -x "${HISAT_CXT}" \
        -1 "${TRIM_DIR}/${line}" -2 "${TRIM_DIR}/${TWO}" \
        --summary-file "${SUM_DIR}/${OUT}cxt_summary.txt" \
        -S "${RAW1}" \
        || { echo "ERROR: Cxt hisat2 failed for ${OUT}"; exit 1; }

      CXT_R1="${CXT_UNAL_DIR}/${OUT}unaligned_R1.fastq.gz"
      CXT_R2="${CXT_UNAL_DIR}/${OUT}unaligned_R2.fastq.gz"
      split_sam "${RAW1}" "${CXT_SAM_DIR}/${OUT}.sam" "${CXT_R1}" "${CXT_R2}" \
        || { echo "ERROR: Cxt split failed for ${OUT}"; exit 1; }
      [ -f "${CXT_R1}" ] && [ -f "${CXT_R2}" ] \
        || { echo "ERROR: missing Cxt-unaligned FASTQs for ${OUT}"; exit 1; }
      rm -f "${RAW1}"

      # ---------------- Stage 2: RVFV ----------------
      RAW2="${RVFV_SAM_DIR}/${OUT}.raw.sam"
      hisat2 --phred33 --rna-strandness RF -p ${THREADS} \
        --no-spliced-alignment \
        -x "${HISAT_RVFV}" \
        -1 "${CXT_R1}" -2 "${CXT_R2}" \
        --summary-file "${SUM_DIR}/${OUT}rvfv_summary.txt" \
        -S "${RAW2}" \
        || { echo "ERROR: RVFV hisat2 failed for ${OUT}"; exit 1; }

      split_sam "${RAW2}" "${RVFV_SAM_DIR}/${OUT}.sam" \
        "${NEITHER_DIR}/${OUT}unaligned_R1.fastq.gz" "${NEITHER_DIR}/${OUT}unaligned_R2.fastq.gz" \
        || { echo "ERROR: RVFV split failed for ${OUT}"; exit 1; }
      [ -f "${NEITHER_DIR}/${OUT}unaligned_R1.fastq.gz" ] && [ -f "${NEITHER_DIR}/${OUT}unaligned_R2.fastq.gz" ] \
        || { echo "ERROR: missing Neither FASTQs for ${OUT}"; exit 1; }
      rm -f "${RAW2}"

      date
    fi
    linenum=$((linenum + 1))
done < "$filename"
