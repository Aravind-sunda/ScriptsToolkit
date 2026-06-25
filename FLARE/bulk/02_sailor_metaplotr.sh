#!/bin/bash
#SBATCH --job-name=sailor_metaplotr_02
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=4:00:00
#SBATCH --output=slurm_%u_%x_%j.log

set -euo pipefail

# All variables exported from 00_saiilor_bulk_pipeline.sh:
#   OUTPUT_DIR, METAPLOT_DIR, GENOME, REFSEQ_DIR, METAPLOT_HELPER

# ── GENEPRED LOOKUP ───────────────────────────────────────────────────────────
declare -A GENEPRED_FILES=(
    [hg19]="$REFSEQ_DIR/hg19_ncbiRefSeqCurated.txt.gz"
    [hg38]="$REFSEQ_DIR/hg38_ncbiRefSeqCurated.txt.gz"
    [hg38_V44]="$REFSEQ_DIR/hg38_V44_E110_basic_knownGene_genePred.txt.gz"
    [mm10]="$REFSEQ_DIR/mm10_ncbiRefSeqCurated.txt.gz"
    [mm39]="$REFSEQ_DIR/mm39_ncbiRefSeqCurated.txt.gz"
)

GENEPRED="${GENEPRED_FILES[$GENOME]:-}"

if [[ -z "$GENEPRED" ]]; then
    echo "[ERROR] No genePred file mapped for GENOME='$GENOME'"
    echo "[ERROR] Supported: hg19, hg38, hg38_V44, mm10, mm39"
    exit 1
fi

if [[ ! -f "$GENEPRED" ]]; then
    echo "[ERROR] genePred file not found: $GENEPRED"
    exit 1
fi

mkdir -p "$METAPLOT_DIR"

echo "---------------------"
echo "[LOG] Starting SAILOR metaPlotR at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] OUTPUT_DIR   : $OUTPUT_DIR"
echo "[LOG] METAPLOT_DIR : $METAPLOT_DIR"
echo "[LOG] GENOME       : $GENOME"
echo "[LOG] GENEPRED     : $GENEPRED"
echo "---------------------"

module load mamba
mamba activate bioinformatics

# Combined BED pattern: {sample}.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed
# Column 4 = Bayesian confidence score (0–1)

for BED in "$OUTPUT_DIR"/*.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed; do

    [[ -f "$BED" ]] || { echo "[WARN] No combined BED files found in $OUTPUT_DIR"; break; }

    SAMPLE=$(basename "$BED" .combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed)
    SDIR="$METAPLOT_DIR/$SAMPLE"
    mkdir -p "$SDIR"

    echo "[INFO] Processing $SAMPLE ..."

    # ── Sort ──────────────────────────────────────────────────────────────────
    SORTED="$SDIR/${SAMPLE}.all.sorted.bed"
    sort -k1,1 -k2,2n "$BED" > "$SORTED"

    # ── Confidence filters ────────────────────────────────────────────────────
    CONF05="$SDIR/${SAMPLE}.conf0.5.sorted.bed"
    CONF09="$SDIR/${SAMPLE}.conf0.9.sorted.bed"
    awk '$4 >= 0.5' "$SORTED" > "$CONF05"
    awk '$4 >= 0.9' "$SORTED" > "$CONF09"

    echo "[INFO] $(wc -l < "$SORTED") total | $(wc -l < "$CONF05") conf>=0.5 | $(wc -l < "$CONF09") conf>=0.9"

    # ── metaPlotR distances ───────────────────────────────────────────────────
    # For strandedness 0 (unstranded): use --ignore-strand because SAILOR
    # assigns strand from chemistry (C>T/G>A), not transcript orientation.
    # For strandedness 1/2: strand is meaningful; do not ignore it.
    for TIER in "$SORTED" "$CONF05" "$CONF09"; do
        LABEL=$(basename "$TIER" .sorted.bed | sed "s/${SAMPLE}\.//")
        DIST_OUT="$SDIR/${SAMPLE}.${LABEL}.dist.measures.txt"

        echo "[INFO] metaPlotR distances: $LABEL ..."
        if [[ "$STRANDEDNESS" -eq 0 ]]; then
            python "$METAPLOT_HELPER" \
                --genePred    "$GENEPRED" \
                --bed         "$TIER" \
                --out         "$DIST_OUT" \
                --ignore-strand
        else
            python "$METAPLOT_HELPER" \
                --genePred    "$GENEPRED" \
                --bed         "$TIER" \
                --out         "$DIST_OUT"
        fi

        echo "[DONE] $DIST_OUT"
    done

done

echo "---------------------"
echo "[LOG] SAILOR metaPlotR completed at $(date '+%Y-%m-%d %H:%M:%S')"
echo "---------------------"
