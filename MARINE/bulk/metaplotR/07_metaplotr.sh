#!/bin/bash
#SBATCH --job-name=metaplotr_07
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=4:00:00
#SBATCH --output=slurm_%u_%x_%j.log

# SAMPLESHEET="/path/to/samplesheet.csv"  # columns: sample,r1,r2,libraryType (SE or PE)
# HOMEDIR="/path/to/project"

# Required — paths to pre-built annotation files (from build_metaplotr_annotations.sh)
# ANNOT_DIR="/path/to/metaplotr_annotations"
# GENOME="hg38"   # must match a subdirectory in ANNOT_DIR (hg19, hg38, mm10)

# Path to the metaplotr Singularity image
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
SIF="$SCRIPT_DIR/metaplotr.sif"

EDIT_TYPE="${EDIT_TYPE:-C>T}"    # must match the edit_type used in 06_marine_filter_normalize.sh
INPUTDIR="$HOMEDIR/06_filter_normalize"
OUTPUTDIR="$HOMEDIR/07_metaplotr"

# ANNOT_DIR, GENOME exported from 00_marine_bulk_pipeline.sh

mkdir -p $OUTPUTDIR

echo "---------------------"
echo "[LOG] Starting metaPlotR at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] Samplesheet : $SAMPLESHEET"
echo "[LOG] GENOME      : $GENOME"
echo "[LOG] ANNOT_DIR   : $ANNOT_DIR"
echo "[LOG] EDIT_TYPE   : $EDIT_TYPE"
echo "[LOG] SIF         : $SIF"
echo "---------------------"

ANNOT_BED="$ANNOT_DIR/$GENOME/annot.bed"
SIZES="$ANNOT_DIR/$GENOME/sizes.txt"

if [[ ! -f "$ANNOT_BED" ]]; then
    echo "[ERROR] Annotation BED not found: $ANNOT_BED"
    echo "[ERROR] Run build_metaplotr_annotations.sh first"
    exit 1
fi

if [[ ! -f "$SIZES" ]]; then
    echo "[ERROR] Sizes file not found: $SIZES"
    echo "[ERROR] Run build_metaplotr_annotations.sh first"
    exit 1
fi

EDIT_TAG="${EDIT_TYPE//>/_}"   # e.g. C>T → C_T

while IFS=',' read -r sample r1 r2 libraryType; do
    [[ "$sample" == "sample" ]] && continue

    echo "[INFO] Processing $sample ..."

    BED6="$INPUTDIR/$sample/bedgraphs/${sample}.${EDIT_TAG}.edit_fraction.bed"

    if [[ ! -f "$BED6" ]]; then
        echo "[WARN] BED6 not found for $sample, skipping: $BED6"
        continue
    fi

    SAMPLE_OUTDIR="$OUTPUTDIR/$sample"
    mkdir -p $SAMPLE_OUTDIR

    SORTED_BED="$SAMPLE_OUTDIR/${sample}.sorted.bed"
    ANNOTATED_BED="$SAMPLE_OUTDIR/${sample}.annotated.bed"
    DISTANCES="$SAMPLE_OUTDIR/${sample}.distances.txt"

    # ── Step 1: Sort BED6 ────────────────────────────────────────────────────
    echo "[INFO] Step 1: Sorting BED6 for $sample ..."
    sort -k1,1 -k2,2n "$BED6" > "$SORTED_BED"

    # ── Step 2: Annotate with gene features ──────────────────────────────────
    echo "[INFO] Step 2: Annotating BED with gene features ..."
    singularity exec --bind "$SAMPLE_OUTDIR,$ANNOT_DIR" "$SIF" \
        perl /opt/metaplotr/annotate_bed_file.pl \
            --bed "$SORTED_BED" \
            --bed2 "$ANNOT_BED" \
        > "$ANNOTATED_BED"

    # ── Step 3: Calculate relative and absolute distances ────────────────────
    echo "[INFO] Step 3: Calculating metagene distances ..."
    singularity exec --bind "$SAMPLE_OUTDIR,$ANNOT_DIR" "$SIF" \
        perl /opt/metaplotr/rel_and_abs_dist_calc.pl \
            --bed "$ANNOTATED_BED" \
            --regions "$SIZES" \
        > "$DISTANCES"

    # ── Step 4: Generate metagene plots ──────────────────────────────────────
    echo "[INFO] Step 4: Generating metagene plots ..."
    singularity exec --bind "$SAMPLE_OUTDIR" "$SIF" \
        Rscript /opt/metaplotr/visualize_metagene.R \
            "$DISTANCES" \
            "$SAMPLE_OUTDIR"

    echo "[DONE] metaPlotR complete for $sample → $SAMPLE_OUTDIR"

done < "$SAMPLESHEET"

echo "[DONE] metaPlotR completed at $(date '+%Y-%m-%d %H:%M:%S')"
