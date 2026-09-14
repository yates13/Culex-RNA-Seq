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

echo "[$0] $SLURM_JOB_NAME $@" # log the command line

# Clear existing modules and initialize Conda
module purge
# --- Load conda (avoids @ symbol issues) ---
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3"; exit 1; }
conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/rnaPseudo_clean || { echo "ERROR: conda activate rnaPseudo_clean failed"; exit 1; }

# Print job timestamp
date

#Run this first to generate a list of files (from inside the CxtHisat dir produced by 02):
#ls | grep -v / | grep ".sam" | grep "Cxt_"  > CxtAligned.txt
#wc -l CxtAligned.txt
#sbatch --array=0-N 04_GTFnew_runAeHTseqArray.sh
#Where N is the number of lines in your CxtAligned.txt file minus 1.

# --- Path setup for new directory structure ---
BASE_DIR="/scratch/alpine/c832500103@colostate.edu/RNA_Seq"
HISAT_DIR="${BASE_DIR}/Data/slurm_script/CxtHisat"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}" .sh)"
OUT_DIR="${BASE_DIR}/Data/${SCRIPT_NAME}/hsCountsCxt_NewGTF"


filename="${HISAT_DIR}/CxtAligned.txt"

# Ensure the output directory exists before running htseq-count
mkdir -vp "${OUT_DIR}"

linenum=0
while IFS= read -r line || [ -n "$line" ]; do
    if [ "$SLURM_ARRAY_TASK_ID" -eq "$linenum" ]; then
      # Build output filename from inpt .sam basename

      base="$(basename "${line}" .sam)"
      suff="${OUT_DIR}/${base}.HSCounts.txt"


      echo "Processing file: ${line}"
      echo "Output path:    ${suff}"

      htseq-count --stranded=reverse --type=gene --idattr=gene_id \
        "${HISAT_DIR}/${line}" \
        /scratch/alpine/c832500103@colostate.edu/RNA_Seq/Genome/Updated_GTF_With_Notes/basefeatures_updated_high_confidence.with_notes.gtf \
        > "$suff" || { echo "ERROR: htseq-count failed for ${line}"; exit 1; }
#      htseq-count --stranded=reverse \
#        "${HISAT_DIR}/${line}" \
#        /scratch/alpine/c832500103@colostate.edu/RNA_Seq/Genome/GTF/Culex-tarsalis_knwr_BASEFEATURES_CtarK1.1_filtered_WNV_rep_fixed.gtf \
#        > "$suff" || { echo "ERROR: htseq-count failed for ${line}"; exit 1; }

#      htseq-count --stranded=reverse \
#        "${HISAT_DIR}/${line}" \
#        /scratch/alpine/camcorey@colostate.edu/CR_RNASeq/index/Culex-tarsalis_knwr_BASEFEATURES_CtarK1.1_filtered_WNV_rep_fixed.gtf \
#        > "$suff"
      break # Exit loop early once matching task ID is processed
    fi
    linenum=$((linenum + 1))
done < "$filename"
