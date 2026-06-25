#!/bin/bash

#SBATCH --job-name=annotate_macs3_peaks
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --output=slurm_%u_%x_%j.log


############################
# USER-EDITABLE VARIABLES
############################

# One or more directories that contain MACS3 narrowPeak files
INPUT_DIRS=(
  "/path/to/macs3_peaks/run1"
)

OUT_DIR="/path/to/annotated_peaks"
GENOME_OR_FASTA="/path/to/genome.fa"   # or "hg38"
GTF="/path/to/genes.gtf"
CPUS=$SLURM_CPUS_PER_TASK

############################
# LOAD MODULES AND RUN
############################

module load mamba
mamba activate
mamba activate homer

mkdir -p "${OUT_DIR}"

for in_dir in "${INPUT_DIRS[@]}"; do
  dir_name="$(basename "${in_dir%/}")"
  out_subdir="${OUT_DIR}/${dir_name}"
  mkdir -p "${out_subdir}"

  for peak in "${in_dir%/}"/*.narrowPeak; do
    base="$(basename "$peak")"
    prefix="${base%.narrowPeak}"

    annotatePeaks.pl \
      "$peak" \
      "$GENOME_OR_FASTA" \
      -gtf "$GTF" \
      -cpu "$CPUS" \
      > "${out_subdir}/${prefix}.annotatePeaks.txt"
  done
done
