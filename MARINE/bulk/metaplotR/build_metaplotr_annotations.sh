#!/bin/bash
# build_metaplotr_annotations.sh
# Pre-build metaPlotR annotation files for common genomes (hg19, hg38, mm10).
# Run this once; the output files are then passed to 07_metaplotr.sh via ANNOT_DIR.
#
# Usage:
#   bash build_metaplotr_annotations.sh <sif_path> <output_dir> [genome1 genome2 ...]
#
# Example:
#   bash build_metaplotr_annotations.sh metaplotr.sif /data/metaplotr_annot hg38 mm10

set -euo pipefail

SIF="${1:?Usage: $0 <metaplotr.sif> <output_dir> [hg19 hg38 mm10 ...]}"
OUTDIR="${2:?Usage: $0 <metaplotr.sif> <output_dir> [hg19 hg38 mm10 ...]}"
shift 2
GENOMES=("${@:-hg19 hg38 mm10}")

declare -A REFGENE_URLS=(
    [hg19]="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
    [hg38]="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
    [mm10]="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz"
)

mkdir -p "$OUTDIR"

echo "====================================="
echo "[LOG] Building metaPlotR annotations"
echo "[LOG] SIF     : $SIF"
echo "[LOG] OUTDIR  : $OUTDIR"
echo "[LOG] Genomes : ${GENOMES[*]}"
echo "====================================="

for GENOME in "${GENOMES[@]}"; do

    if [[ -z "${REFGENE_URLS[$GENOME]+x}" ]]; then
        echo "[WARN] No known refGene URL for $GENOME — skipping. Add it to REFGENE_URLS in this script."
        continue
    fi

    GENOME_DIR="$OUTDIR/$GENOME"
    mkdir -p "$GENOME_DIR"

    REFGENE="$GENOME_DIR/refGene.txt"
    ANNOT_BED="$GENOME_DIR/annot.bed"
    SIZES="$GENOME_DIR/sizes.txt"

    echo ""
    echo "-------------------------------------"
    echo "[INFO] Processing $GENOME ..."
    echo "-------------------------------------"

    # Download refGene if not already present
    if [[ ! -f "$REFGENE" ]]; then
        echo "[INFO] Downloading refGene for $GENOME ..."
        wget -q -O "${REFGENE}.gz" "${REFGENE_URLS[$GENOME]}"
        gunzip "${REFGENE}.gz"
        echo "[DONE] Downloaded: $REFGENE"
    else
        echo "[INFO] refGene already exists, skipping download: $REFGENE"
    fi

    # Build annotation BED
    echo "[INFO] Building annotation BED for $GENOME ..."
    singularity exec "$SIF" \
        perl /opt/metaplotr/make_annot_bed.pl --genome "$REFGENE" \
        > "$ANNOT_BED"
    echo "[DONE] Annotation BED: $ANNOT_BED"

    # Build CDS/UTR sizes file
    echo "[INFO] Building sizes file for $GENOME ..."
    singularity exec "$SIF" \
        perl /opt/metaplotr/size_of_cds_utrs.pl --genome "$REFGENE" \
        > "$SIZES"
    echo "[DONE] Sizes file: $SIZES"

    echo "[DONE] $GENOME annotations complete → $GENOME_DIR"

done

echo ""
echo "====================================="
echo "[DONE] All annotations built at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[DONE] Pass ANNOT_DIR=$OUTDIR to 07_metaplotr.sh"
echo "====================================="
