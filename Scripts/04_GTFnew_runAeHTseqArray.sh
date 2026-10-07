#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=htseq
#SBATCH --time=18:00:00
#SBATCH --mem=32G
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shawn.yates@colostate.edu
#SBATCH --output=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog/%x.%A-%a.log
#SBATCH --error=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog/%x.%A-%a.err

# =============================================================================
# 04_GTFnew_runAeHTseqArray.sh
#
# Purpose:
#   SLURM array job. Each array task runs htseq-count on one aligned SAM
#   file (from step 02, HISAT2) against the updated gene-only GTF
#   (basefeatures_updated_high_confidence.with_notes.gtf), producing a
#   per-gene read count file.
#
# NOTE ON GTF: this GTF contains ONLY "gene" feature rows (no exon/CDS/
#   transcript rows), confirmed via `cut -f3 | sort | uniq -c`. Because of
#   this, --type=gene --idattr=gene_id is used instead of htseq-count's
#   default (--type=exon), so reads are counted against the FULL gene body
#   span rather than exon-only regions. This means intronic/UTR reads are
#   now included in gene counts — a reasonable tradeoff given the available
#   annotation, but worth stating explicitly in methods.
#
# Input:
#   ${HISAT_DIR}/CxtAligned.txt   - list of .sam filenames, one per line,
#                                    line number (0-indexed) = array task ID
#   ${HISAT_DIR}/<sample>.sam     - aligned reads from step 02
#   GTF_FILE (see below)          - gene-body annotation for counting
#
# Output:
#   ${BASE_DIR}/Data/04_GTFnew_runAeHTseqArray/hsCountsCxt_NewGTF/<sample>.HSCounts.txt
#
# Setup (run once before first submission, from inside the CxtHisat dir
# produced by step 02 — regenerate this fresh any time CxtHisat changes):
#   ls | grep -v / | grep ".sam" | grep "Cxt_" > CxtAligned.txt
#   wc -l CxtAligned.txt
#
# Usage:
#   sbatch --array=0-<N> 04_GTFnew_runAeHTseqArray.sh
#   where N = (number of lines in CxtAligned.txt) - 1
# =============================================================================


echo "[$0] $SLURM_JOB_NAME $@" # log the command line

# --- Pre-flight: ensure SLURM log directories exist -------------------------
mkdir -vp /scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog
mkdir -vp /scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog

# Clear existing modules and initialize Conda
module purge

# --- Load conda (avoids @ symbol issues) ---
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3"; exit 1; }
conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/rnaPseudo_clean2 || { echo "ERROR: conda activate rnaPseudo_clean failed"; exit 1; }

# Print job timestamp
date

# --- Path setup for new directory structure ---
BASE_DIR="/scratch/alpine/c832500103@colostate.edu/RNA_Seq"
HISAT_DIR="${BASE_DIR}/Data/slurm_script/CxtHisat"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}" .sh)"
OUT_DIR="${BASE_DIR}/Data/${SCRIPT_NAME}/hsCountsCxt_NewGTF"
GTF_FILE="${BASE_DIR}/Genome/Updated_GTF_With_Notes/basefeatures_updated_high_confidence.with_notes.gtf"

filename="${HISAT_DIR}/CxtAligned.txt"

# --- Sanity checks before doing any work -------------------------------------
if [ ! -d "${HISAT_DIR}" ]; then
    echo "ERROR: HISAT_DIR does not exist: ${HISAT_DIR}"
    exit 1
fi

if [ ! -f "${filename}" ]; then
    echo "ERROR: Sample list not found: ${filename}"
    echo "  Generate it first (see header comment) before submitting this array job."
    exit 1
fi

if [ ! -f "${GTF_FILE}" ]; then
    echo "ERROR: GTF file not found: ${GTF_FILE}"
    exit 1
fi

# Ensure the output directory exists before running htseq-count
mkdir -vp "${OUT_DIR}"

linenum=0
while IFS= read -r line || [ -n "$line" ]; do
    if [ "$SLURM_ARRAY_TASK_ID" -eq "$linenum" ]; then
      # Build output filename from inpt .sam basename
      base="$(basename "${line}" .sam)"
      suff="${OUT_DIR}/${base}.HSCounts.txt"

      # Confirm the SAM file actually exists before running htseq-count on it
      if [ ! -f "${HISAT_DIR}/${line}" ]; then
          echo "ERROR: SAM file not found for task ${SLURM_ARRAY_TASK_ID}: ${HISAT_DIR}/${line}"
          exit 1
      fi

      echo "Processing file: ${line}"
      echo "Output path:    ${suff}"

      htseq-count --stranded=reverse --type=gene --idattr=gene_id \
        "${HISAT_DIR}/${line}" \
        "${GTF_FILE}" \
        > "$suff" || { echo "ERROR: htseq-count failed for ${line}"; exit 1; }

      break # Exit loop early once matching task ID is processed
    fi
    linenum=$((linenum + 1))
done < "$filename"
