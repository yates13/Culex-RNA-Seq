#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=fastpLoop
#SBATCH --output=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/outlog/%x.%j.out
#SBATCH --error=/scratch/alpine/c832500103@colostate.edu/RNA_Seq/slurmlogs/errlog/%x.%j.err
#SBATCH --time=5:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=17
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shawn.yates@colostate.edu

SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3"; exit 1; }
conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/cellSquito || { echo "ERROR: conda activate cellSquito failed"; exit 1; }

# Fix for fastp's libisal.so.2 dependency not being found on the default library path
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

BASE_DIR="/scratch/alpine/c832500103@colostate.edu/RNA_Seq"
RAW_DIR="${BASE_DIR}/Raw_Files"

# Derive this script's own name (e.g. "01b_runFastp_RNA_HTML") for the output folder
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}" .sh)"
OUT_DIR="${BASE_DIR}/Data/${SCRIPT_NAME}/trimmed"

mkdir -vp "${OUT_DIR}/htmls"

for FILE in "${RAW_DIR}"/*R1_001.fastq.gz; do
    echo "Processing R1: ${FILE}"
    TWO="${FILE/_R1/_R2}"
    echo "Processing R2: ${TWO}"

    TRIM1="${OUT_DIR}/trimmed.$(basename "${FILE}")"
    TRIM2="${OUT_DIR}/trimmed.$(basename "${TWO}")"

    cmd="fastp -i ${FILE} -I ${TWO} -o ${TRIM1} -O ${TRIM2} -h ${OUT_DIR}/htmls/$(basename "${TWO}").html -j ${OUT_DIR}/htmls/$(basename "${TWO}").json -w $((SLURM_NTASKS-1))"
    echo "Running: ${cmd}"
    eval "${cmd}"
done
