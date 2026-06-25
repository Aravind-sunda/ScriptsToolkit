#!/usr/bin/env bash

#SBATCH --job-name=cr_count
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --time=48:00:00
#SBATCH --array=0-1
#SBATCH --output=slurm_%x_%A_%a.log

set -euo pipefail

module load cellranger/9.0.1   # or set PATH to your install

REF="/path/to/refdata-gex-XXXX"
SAMPLES_CSV="/path/to/samples.csv"

# get line (skip header). array 0 -> line 2

LINE="$(awk -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print; exit}' "$SAMPLES_CSV")"
IFS=, read -r ID SAMPLE FASTQS <<< "$LINE"

cellranger count \
  --id="$ID" \
  --transcriptome="$REF" \
  --fastqs="$FASTQS" \
  --sample="$SAMPLE" \
  --create-bam=true \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --expect-cells=10000

# NOTE:
# SAMPLE has to be prefix of FASTQS files. If not, cellranger will not find the files and fail.
# Change array number in SBATCH --array=0-1 to match number of samples in count_samplesheet.csv (minus header).
# add the following 2 commands if you would like to use cellranger's auto cell annotation and have a tenx cloud token:
# --cell-annotation-model auto
# --tenx-cloud-token-path <PATH> 