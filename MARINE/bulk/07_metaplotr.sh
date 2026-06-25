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

# EDIT_TYPE, GENOME, and REFSEQ_DIR exported from 00_marine_bulk_pipeline.sh

SCRIPT_DIR="$(dirname "$(realpath "$0")")"

INPUTDIR="$HOMEDIR/06_filter_normalize"
OUTPUTDIR="$HOMEDIR/07_metaplotr"

# ── genePred file lookup ───────────────────────────────────────────────────────
declare -A GENEPRED_FILES=(
	[hg19]="$REFSEQ_DIR/hg19_ncbiRefSeqCurated.txt.gz"
	[hg38]="$REFSEQ_DIR/hg38_ncbiRefSeqCurated.txt.gz"
	[hg38_V44]="$REFSEQ_DIR/hg38_V44_E110_basic_knownGene_genePred.txt.gz"
	[mm10]="$REFSEQ_DIR/mm10_ncbiRefSeqCurated.txt.gz"
	[mm39]="$REFSEQ_DIR/mm39_ncbiRefSeqCurated.txt.gz"
)

GENEPRED="${GENEPRED_FILES[$GENOME]:-}"

mkdir -p $OUTPUTDIR

echo "---------------------"
echo "[LOG] Starting metaPlotR distance calculation at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] Samplesheet : $SAMPLESHEET"
echo "[LOG] GENOME      : $GENOME"
echo "[LOG] GENEPRED    : $GENEPRED"
echo "[LOG] EDIT_TYPE   : $EDIT_TYPE"
echo "[LOG] INPUTDIR    : $INPUTDIR"
echo "[LOG] OUTPUTDIR   : $OUTPUTDIR"
echo "---------------------"

if [[ -z "$GENEPRED" ]]; then
	echo "[ERROR] No genePred file mapped for GENOME='$GENOME'"
	echo "[ERROR] Supported: hg19, hg38, hg38_V44, mm10, mm39"
	exit 1
fi

if [[ ! -f "$GENEPRED" ]]; then
	echo "[ERROR] genePred file not found: $GENEPRED"
	exit 1
fi

module load mamba
mamba activate bioinformatics

EDIT_TAG="${EDIT_TYPE//>/_}"   # e.g. C>T → C_T

while IFS=',' read -r sample r1 r2 libraryType; do
	[[ "$sample" == "sample" ]] && continue

	echo "[INFO] Processing $sample ..."

	BED6="$INPUTDIR/$sample/normalized/bedgraphs/${sample}.${EDIT_TAG}.edit_fraction.bed"

	if [[ ! -f "$BED6" ]]; then
		echo "[WARN] BED6 not found for $sample, skipping: $BED6"
		continue
	fi

	SAMPLE_OUTDIR="$OUTPUTDIR/$sample"
	mkdir -p $SAMPLE_OUTDIR

	SORTED_BED="$SAMPLE_OUTDIR/${sample}.${EDIT_TAG}.sorted.bed"
	DIST_OUT="$SAMPLE_OUTDIR/${sample}.${EDIT_TAG}.dist.measures.txt"

	# ── Sort BED6 ─────────────────────────────────────────────────────────────
	echo "[INFO] Sorting BED6 for $sample ..."
	sort -k1,1 -k2,2n "$BED6" > "$SORTED_BED"

	# ── Compute metagene distances ────────────────────────────────────────────
	echo "[INFO] Computing metagene distances for $sample ..."
	python "$SCRIPT_DIR/helper_calc_metaplot_dist.py" \
		--genePred "$GENEPRED" \
		--bed      "$SORTED_BED" \
		--out      "$DIST_OUT"

	echo "[DONE] Distance measures for $sample → $DIST_OUT"

done < "$SAMPLESHEET"

echo "[DONE] metaPlotR distances completed at $(date '+%Y-%m-%d %H:%M:%S')"
