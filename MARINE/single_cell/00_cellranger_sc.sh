#!/bin/bash
#SBATCH --partition=defq
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@institution.edu
#SBATCH --output=slurm_%u_%x_%j.out

# ══════════════════════════════════════════════════════════════════════════════
# 00_cellranger_sc.sh
# Per-sample SLURM job: run CellRanger count for one single-cell sample.
# Submitted by 00_marine_sc_pipeline.sh when START_FROM="fastq"; do not call directly.
#
# Required exports (injected by 00_marine_sc_pipeline.sh):
#   SAMPLENAME, FASTQ_DIR, CELLRANGER_OUTDIR,
#   TRANSCRIPTOME_REF, CELLRANGER_CPUS, CELLRANGER_MEM_GB
#
# ── FASTQ NAMING REQUIREMENT ──────────────────────────────────────────────────
# CellRanger discovers reads by filename convention, not by explicit R1/R2 paths.
# Files inside FASTQ_DIR must follow the 10x Genomics naming scheme:
#
#   {SAMPLENAME}_S{n}_L{lane}_R1_001.fastq.gz   ← barcodes + UMI  (28 bp)
#   {SAMPLENAME}_S{n}_L{lane}_R2_001.fastq.gz   ← cDNA read
#   {SAMPLENAME}_S{n}_L{lane}_I1_001.fastq.gz   ← index read (optional)
#
# If your FASTQs came from Illumina mkfastq or 10x mkfastq they are already
# named correctly. If they were downloaded (SRA/GEO) or have arbitrary names,
# rename them before running this pipeline. Example:
#   mv sample_R1.fastq.gz  SAMPLE1_S1_L001_R1_001.fastq.gz
#   mv sample_R2.fastq.gz  SAMPLE1_S1_L001_R2_001.fastq.gz
# ─────────────────────────────────────────────────────────────────────────────
#
# Output layout: ${CELLRANGER_OUTDIR}/${SAMPLENAME}/outs/
#   possorted_genome_bam.bam          → passed to 01_run_marine_sc.sh as USE_BAM
#   filtered_feature_bc_matrix/       → used as BARCODES source and MATRIXDIR
# ══════════════════════════════════════════════════════════════════════════════

CELLRANGER_CPUS=32      # keep in sync with #SBATCH --cpus-per-task
CELLRANGER_MEM_GB=300   # keep in sync with #SBATCH --mem (integer GB)

echo "===== CellRanger SC: ${SAMPLENAME} ====="
echo "[LOG] FASTQ dir:       ${FASTQ_DIR}"
echo "[LOG] Transcriptome:   ${TRANSCRIPTOME_REF}"
echo "[LOG] Output base:     ${CELLRANGER_OUTDIR}"
echo "[LOG] Cores:           ${CELLRANGER_CPUS}"
echo "[LOG] Memory (GB):     ${CELLRANGER_MEM_GB}"
echo "[LOG] Started:         $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# CellRanger will not overwrite an existing run directory — remove it first
CR_SAMPLE_DIR="${CELLRANGER_OUTDIR}/${SAMPLENAME}"
if [ -d "${CR_SAMPLE_DIR}" ]; then
    echo "[INFO] Removing existing CellRanger output: ${CR_SAMPLE_DIR}"
    rm -rf "${CR_SAMPLE_DIR}"
fi

mkdir -p "${CELLRANGER_OUTDIR}"

# CellRanger creates ${SAMPLENAME}/ inside the current working directory.
# cd here so the output lands in CELLRANGER_OUTDIR.
cd "${CELLRANGER_OUTDIR}"

module load cellranger/9.0.1

echo "[INFO] Running cellranger count..."
cellranger count \
    --id="${SAMPLENAME}" \
    --transcriptome="${TRANSCRIPTOME_REF}" \
    --fastqs="${FASTQ_DIR}" \
    --sample="${SAMPLENAME}" \
    --localcores="${CELLRANGER_CPUS}" \
    --localmem="${CELLRANGER_MEM_GB}" \
    --expect-cells=10000 \
    --create-bam=true \
    --cell-annotation-model=auto

echo ""
echo "[DONE] ${SAMPLENAME} completed at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[INFO] BAM:       ${CR_SAMPLE_DIR}/outs/possorted_genome_bam.bam"
echo "[INFO] Matrix:    ${CR_SAMPLE_DIR}/outs/filtered_feature_bc_matrix/"
echo "===== End: ${SAMPLENAME} ====="
