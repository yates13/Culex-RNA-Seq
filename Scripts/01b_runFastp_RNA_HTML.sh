#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=fastpLoop
#SBATCH --output=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog/%x.%j.out
#SBATCH --error=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog/%x.%j.err
#SBATCH --time=10:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=17
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shawn.yates@colostate.edu

# =============================================================================
# 01b_runFastp_RNA_HTML.sh
#
# Purpose:
#   Quality-trims paired-end RNA-seq reads with fastp for every R1/R2 pair
#   found in RAW_DIR. Produces trimmed FASTQs plus per-sample HTML/JSON QC
#   reports.
#
# Input:
#   ${RAW_DIR}/*_R1_001.fastq.gz  (and matching *_R2_001.fastq.gz)
#
# Output:
#   ${BASE_DIR}/Data/01b_runFastp_RNA_HTML/trimmed/
#     trimmed.<sample>_R1_001.fastq.gz
#     trimmed.<sample>_R2_001.fastq.gz
#     htmls/<sample>_R2_001.fastq.gz.html   (fastp QC report)
#     htmls/<sample>_R2_001.fastq.gz.json   (fastp QC data)
#
# Usage:
#   sbatch 01b_runFastp_RNA_HTML.sh
#
# Notes:
#   - Output folder name is derived automatically from this script's own
#     filename (via SCRIPT_NAME), so renaming the script changes where
#     output is written.
#   - SLURM writes its own --output/--error log files BEFORE any commands
#     in this script run, so slurmlogs/outlog and slurmlogs/errlog must
#     already exist or the job will fail at submission/startup with no
#     useful error message. The pre-flight check below creates them if
#     missing, but note that SLURM itself reads the #SBATCH paths above
#     before this script body executes — if this is truly the FIRST run
#     and those folders don't exist yet, create them manually once before
#     submitting:
#       mkdir -p /scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog
#       mkdir -p /scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog
# =============================================================================

echo "[$0] $SLURM_JOB_NAME $@"  # log the command line for this run
date

# --- Pre-flight: ensure SLURM log directories exist -------------------------
# (Belt-and-suspenders — see note above. If these are missing on a FIRST
#  submission, SLURM may already have failed before reaching this line;
#  this mainly guards against them being deleted between runs.)
mkdir -vp /scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog
mkdir -vp /scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog

# --- Load conda environment --------------------------------------------------
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3"; exit 1; }
conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/cellSquito || { echo "ERROR: conda activate cellSquito failed"; exit 1; }

# Fix for fastp's libisal.so.2 dependency not being found on the default library path
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

# --- Path setup ---------------------------------------------------------------
BASE_DIR="/scratch/alpine/c832500103@colostate.edu/RNA_Seq"
RAW_DIR="${BASE_DIR}/Raw_Files"

# Derive this script's own name (e.g. "01b_runFastp_RNA_HTML") for the output folder
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}" .sh)"
OUT_DIR="${BASE_DIR}/Data/${SCRIPT_NAME}/trimmed"

# --- Sanity check: confirm raw input directory exists and has R1 files -------
if [ ! -d "${RAW_DIR}" ]; then
    echo "ERROR: RAW_DIR does not exist: ${RAW_DIR}"
    exit 1
fi

shopt -s nullglob
r1_files=( "${RAW_DIR}"/*R1_001.fastq.gz )
shopt -u nullglob

if [ ${#r1_files[@]} -eq 0 ]; then
    echo "ERROR: No *_R1_001.fastq.gz files found in ${RAW_DIR}"
    exit 1
fi

# --- Create output directories (trimmed FASTQs + HTML/JSON QC reports) ------
mkdir -vp "${OUT_DIR}/htmls"

for FILE in "${RAW_DIR}"/*R1_001.fastq.gz; do
    echo "Processing R1: ${FILE}"
    TWO="${FILE/_R1/_R2}"
    echo "Processing R2: ${TWO}"

    # Confirm the matching R2 file actually exists before running fastp on it
    if [ ! -f "${TWO}" ]; then
        echo "ERROR: Matching R2 file not found for ${FILE} (expected: ${TWO}) — skipping this pair."
        continue
    fi

    TRIM1="${OUT_DIR}/trimmed.$(basename "${FILE}")"
    TRIM2="${OUT_DIR}/trimmed.$(basename "${TWO}")"

    cmd="fastp -i ${FILE} -I ${TWO} -o ${TRIM1} -O ${TRIM2} -h ${OUT_DIR}/htmls/$(basename "${TWO}").html -j ${OUT_DIR}/htmls/$(basename "${TWO}").json -w $((SLURM_NTASKS-1))"
    echo "Running: ${cmd}"
    eval "${cmd}"
done
