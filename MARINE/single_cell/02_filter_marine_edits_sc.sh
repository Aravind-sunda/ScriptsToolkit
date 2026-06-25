#!/bin/bash
# 02_filter_marine_edits_sc.sh
# Head-node script: filter MARINE SC edit calls for all samples against dbSNP
# and apply minimum coverage + maximum editing-fraction thresholds.
#
# Reads project variables exported by 00_marine_sc_pipeline.sh.
# Run manually on the head node after all 01_run_marine_sc.sh jobs complete.
#
# Usage: bash 02_filter_marine_edits_sc.sh

# set -euo pipefail

# If called stand-alone (not via 00_marine_sc_pipeline.sh), set variables here:
# HOMEDIR="/path/to/project"
# SAMPLESHEET="${HOMEDIR}/samplesheet.csv"
# MARINE_OUTDIR="${HOMEDIR}/01_marine"
# FILTER_OUTDIR="${HOMEDIR}/02_filter"
# DBSNP_BED="/path/to/dbsnp_combined.bed3"
# MIN_COUNT=3
# MAX_FRAC=0.1

module load mamba
mamba activate datascience

mkdir -p "${FILTER_OUTDIR}"

echo "===== Filter MARINE SC Edits ====="
echo "[LOG] Started:       $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] MARINE_OUTDIR: ${MARINE_OUTDIR}"
echo "[LOG] FILTER_OUTDIR: ${FILTER_OUTDIR}"
echo "[LOG] DBSNP_BED:     ${DBSNP_BED}"
echo "[LOG] MIN_COUNT:     ${MIN_COUNT}  MAX_FRAC: ${MAX_FRAC}"
echo "=================================="

while IFS=',' read -r samplename bamfile matrixdir; do

    samplename=$(echo "${samplename}" | tr -d '\r')

    sample_indir="${MARINE_OUTDIR}/${samplename}"
    sample_outdir="${FILTER_OUTDIR}/${samplename}"

    # use annotated output for single cell and strandedness 2
    marine_file="${sample_indir}/final_filtered_site_info_annotated.tsv"

    if [ ! -f "${marine_file}" ]; then
        echo "[WARN] MARINE output not found for ${samplename} — skipping (${sample_indir})"
        continue
    fi

    echo "[INFO] ${samplename} — input: ${marine_file}"
    mkdir -p "${sample_outdir}"

    python helper_filter_edits_sc.py \
        --marine-results "${marine_file}" \
        --dbsnp-bed      "${DBSNP_BED}" \
        --min-count      "${MIN_COUNT}" \
        --max-frac       "${MAX_FRAC}" \
        --output-dir     "${sample_outdir}"

    echo "[DONE] ${samplename} → ${sample_outdir}/filtered_edits.tsv"

done < <(tail -n +2 "${SAMPLESHEET}")

echo ""
echo "[DONE] Filtering complete at $(date '+%Y-%m-%d %H:%M:%S')"
