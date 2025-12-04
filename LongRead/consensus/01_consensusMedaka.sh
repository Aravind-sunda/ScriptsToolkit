#!/bin/bash



# Author: Aravind Sundaravadivelu
# Created: 2025-12-04
# Description: Merge ONT FASTQs and run Medaka consensus.
# Version: 0.1

usage() {
  cat <<EOF
  
Usage: $(basename "$0") -i INPUTDIR -r DRAFT_REF -o OUTDIR -t THREADS [options]

Required:
  -i  Input directory containing FASTQ files, can be mutliple files of same sample which will be merged
  -r  Draft reference FASTA file
  -o  Output directory
  -t  Number of threads for medaka

Optional:
    -h  Show this help message and exit

Behavior:
  - This script uses a mamba/conda environment named "Medaka".
  - If it does not already exist, it will be created automatically.

Example:
  $(basename "$0") -i /data/fastq -r ref.fa -o /data/medaka_out -t 16

EOF
}
# ===== Parse flags instead of positional args =====
# -i : INPUTDIR
# -r : DRAFT_REF
# -o : OUTDIR
# -t : THREADS
while getopts ":i:r:o:t:h" opt; do
  case "$opt" in
    i) INPUTDIR="$OPTARG" ;;
    r) DRAFT_REF="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h)
      usage
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      exit 1
      ;;
  esac
done

# Basic check that all required flags were passed
if [ -z "${INPUTDIR:-}" ] || [ -z "${DRAFT_REF:-}" ] || \
   [ -z "${OUTDIR:-}" ]   || [ -z "${THREADS:-}" ]; then
  echo "Error: -i INPUTDIR -r DRAFT_REF -o OUTDIR -t THREADS are required." >&2
  echo "Usage: $0 -i INPUTDIR -r DRAFT_REF -o OUTDIR -t THREADS" >&2
  exit 1
fi


# STEP 1: Creating a conda environment and installing the yaml file if the environment does not exist
module load mamba
mamba activate

if mamba info --envs | grep -q Medaka; then
  echo "[INFO] Medaka environment already exists"
else
  mamba create -n Medaka -c conda-forge -c nanoporetech -c bioconda medaka
fi

mamba activate Medaka

# STEP 2: Merging different fastq files together into a single fastq file

# Merging different basecalls together into a single fastq file
cat ${INPUTDIR}/*.fastq > ${OUTDIR}/merged_fastq/merged_fastq.fastq

# Step3: Running medaka consensus
# medaka tools list_models
BASECALLS="${OUTDIR}/merged_fastq/merged_fastq.fastq"
mkdir -p ${OUTDIR}/medaka_consensus_output

MODEL=$(medaka tools resolve_model --auto_model consensus ${BASECALLS})

medaka_consensus -i ${BASECALLS} -d ${DRAFT_REF} -o ${OUTDIR} -t ${THREADS} -m ${MODEL} -r N
