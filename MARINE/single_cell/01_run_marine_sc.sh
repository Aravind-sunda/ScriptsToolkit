#!/bin/bash
#SBATCH --partition=bigmemq
#SBATCH --cpus-per-task=32
#SBATCH --mem=900G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@institution.edu
#SBATCH --output=slurm_%u_%x_%j.out

# ══════════════════════════════════════════════════════════════════════════════
# 01_run_marine_sc.sh
# Per-sample SLURM job: run MARINE in single-cell mode.
# Submitted by 00_marine_sc_pipeline.sh via sbatch --export; do not call directly.
#
# Required exports (injected by 00_marine_sc_pipeline.sh):
#   SAMPLENAME, USE_BAM, BARCODES, SAMPLE_OUTDIR, MARINE_OUTDIR,
#   GENOME_FA, CONTAINER, MARINE_PY, ANNOTATION_BED,
#   BARCODE_TAG, STRANDEDNESS, CORES, MIN_READ_QUALITY, MIN_BASE_QUALITY
# ══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

echo "===== MARINE SC: ${SAMPLENAME} ====="
echo "[LOG] BAM:             ${USE_BAM}"
echo "[LOG] Barcodes:        ${BARCODES}"
echo "[LOG] Output:          ${SAMPLE_OUTDIR}"
echo "[LOG] Strandedness:    ${STRANDEDNESS}"
echo "[LOG] Started:         $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Create the parent directory; do NOT create SAMPLE_OUTDIR itself —
# MARINE creates it and will error if it already exists.
mkdir -p "$(dirname "${SAMPLE_OUTDIR}")"

# Remove any previous run so MARINE can write a clean output directory.
if [ -d "${SAMPLE_OUTDIR}" ]; then
    echo "[INFO] Removing existing output directory: ${SAMPLE_OUTDIR}"
    rm -rf "${SAMPLE_OUTDIR}"
fi

# ── BAM preparation: ensure MD tags are present ───────────────────────────────
# MARINE requires MD tags to call editing sites. CellRanger BAMs often lack them.
# Check the first 100 mapped reads; if none carry MD:Z: run samtools calmd.

CALMD_BAM="${MARINE_OUTDIR}/${SAMPLENAME}.md.bam"

echo "[INFO] Checking for MD tags in BAM..."
MD_TAG_COUNT=$(singularity exec \
    --bind /condo/brannanlab:/condo/brannanlab \
    "${CONTAINER}" \
    bash -c "samtools view -F 4 '${USE_BAM}' 2>/dev/null | head -100 | grep -c 'MD:Z:' || echo 0")

if [ "${MD_TAG_COUNT}" -gt 0 ]; then
    echo "[INFO] MD tags present — using original BAM."
    FINAL_BAM="${USE_BAM}"
else
    echo "[INFO] MD tags absent — running samtools calmd..."

    if [ -f "${CALMD_BAM}" ]; then
        echo "[INFO] calmd BAM already exists, skipping: ${CALMD_BAM}"
    else
        TMP_CALMD="${CALMD_BAM}.tmp"
        singularity exec \
            --bind /condo/brannanlab:/condo/brannanlab \
            "${CONTAINER}" \
            samtools calmd -b "${USE_BAM}" "${GENOME_FA}" \
            > "${TMP_CALMD}"
        mv "${TMP_CALMD}" "${CALMD_BAM}"
        echo "[INFO] calmd complete: ${CALMD_BAM}"
    fi

    FINAL_BAM="${CALMD_BAM}"
fi

# ── BAM indexing ──────────────────────────────────────────────────────────────
# Check both common BAI naming conventions.

BAI_SUFFIXED="${FINAL_BAM}.bai"
BAI_UNSUFFIXED="${FINAL_BAM%.bam}.bai"

if [ ! -f "${BAI_SUFFIXED}" ] && [ ! -f "${BAI_UNSUFFIXED}" ]; then
    echo "[INFO] BAI not found — indexing BAM with samtools..."
    singularity exec \
        --bind /condo/brannanlab:/condo/brannanlab \
        "${CONTAINER}" \
        samtools index -@ "${SLURM_CPUS_PER_TASK}" "${FINAL_BAM}"
    echo "[INFO] Indexing complete."
fi

# ── Run MARINE ────────────────────────────────────────────────────────────────

echo "[INFO] Running MARINE in single-cell mode..."
singularity exec \
    --bind /condo/brannanlab:/condo/brannanlab \
    "${CONTAINER}" \
    python "${MARINE_PY}" \
        --bam_filepath             "${FINAL_BAM}" \
        --output_folder            "${SAMPLE_OUTDIR}" \
        --barcode_whitelist_file   "${BARCODES}" \
        --annotation_bedfile_path  "${ANNOTATION_BED}" \
        --barcode_tag              "${BARCODE_TAG}" \
        --strandedness             "${STRANDEDNESS}" \
        --cores                    "${SLURM_CPUS_PER_TASK}" \
        --keep_intermediate_files \
        --min_read_quality         "${MIN_READ_QUALITY}" \
        --min_base_quality         "${MIN_BASE_QUALITY}"

echo ""
echo "[DONE] ${SAMPLENAME} completed at $(date '+%Y-%m-%d %H:%M:%S')"
echo "===== End: ${SAMPLENAME} ====="
