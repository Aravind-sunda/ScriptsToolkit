#!/bin/bash
# 03_normalize_marine_edits_umi_sc.sh
# Head-node script: normalize filtered MARINE SC edits by UMI counts from the
# CellRanger cell-by-gene matrix for all samples.
#
# Reads project variables exported by 00_marine_sc_pipeline.sh.
# Run manually on the head node after 02_filter_marine_edits_sc.sh completes.
#
# Usage: bash 03_normalize_marine_edits_umi_sc.sh

# set -euo pipefail

# If called stand-alone (not via 00_marine_sc_pipeline.sh), set variables here:
# HOMEDIR="/path/to/project"
# SAMPLESHEET="${HOMEDIR}/samplesheet.csv"
# FILTER_OUTDIR="${HOMEDIR}/02_filter"
# NORMALIZE_OUTDIR="${HOMEDIR}/03_normalize"
# ANNOTATION_BED="/path/to/gene_models.bed"
# EDIT_TYPE="C>T"         # strand conversion to retain: "C>T" or "A>G"
# CELLRANGER_OUTDIR="${HOMEDIR}/00_cellranger"   # needed when START_FROM="fastq"

module load mamba
mamba activate datascience

mkdir -p "${NORMALIZE_OUTDIR}"

echo "===== Normalize MARINE SC Edits (UMI) ====="
echo "[LOG] Started:          $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] FILTER_OUTDIR:    ${FILTER_OUTDIR}"
echo "[LOG] NORMALIZE_OUTDIR: ${NORMALIZE_OUTDIR}"
echo "[LOG] ANNOTATION_BED:   ${ANNOTATION_BED}"
echo "[LOG] EDIT_TYPE:        ${EDIT_TYPE}"
echo "============================================"

while IFS=',' read -r samplename bamfile matrixdir; do

    samplename=$(echo "${samplename}" | tr -d '\r')
    matrixdir=$(echo "${matrixdir}"   | tr -d '\r')

    # In fastq mode the samplesheet has no matrixdir column — derive from CellRanger output
    if [ -z "${matrixdir}" ]; then
        matrixdir="${CELLRANGER_OUTDIR}/${samplename}/outs/filtered_feature_bc_matrix"
    fi

    filtered_edits="${FILTER_OUTDIR}/${samplename}/filtered_edits.tsv"
    sample_outdir="${NORMALIZE_OUTDIR}/${samplename}"

    if [ ! -f "${filtered_edits}" ]; then
        echo "[WARN] Filtered edits not found for ${samplename} — skipping (${filtered_edits})"
        continue
    fi

    echo "[INFO] ${samplename} — filtered: ${filtered_edits} | matrixdir: ${matrixdir}"
    mkdir -p "${sample_outdir}"

    python helper_normalize_edits_umi_sc.py \
        --filtered_edits  "${filtered_edits}" \
        --counts_matrix   "${matrixdir}" \
        --bed             "${ANNOTATION_BED}" \
        --output_dir      "${sample_outdir}" \
        --edit-type       "${EDIT_TYPE}"

    echo "[DONE] ${samplename} → ${sample_outdir}"

done < <(tail -n +2 "${SAMPLESHEET}")

echo ""
echo "[DONE] Normalization complete at $(date '+%Y-%m-%d %H:%M:%S')"
