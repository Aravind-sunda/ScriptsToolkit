#!/bin/bash
#SBATCH --job-name=marine_filter_normalize_06
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%u_%x_%j.log

# SAMPLESHEET="/path/to/samplesheet.csv"  # columns: sample,r1,r2,libraryType (SE or PE)
# HOMEDIR="/path/to/project"

MARINE_DIR="$HOMEDIR/04_marine"
STRANDEDNESS_DIR="$HOMEDIR/04_strandedness"
FC_MATRIX="$HOMEDIR/05_featurecounts/counts_matrix_combined.tsv"   # merged featureCounts output from 05_featurecounts.sh
OUTPUTDIR="$HOMEDIR/06_filter_normalize"

MAX_FRAC="0.10"   # max editing fraction cutoff for filter step
# EDIT_TYPE exported from 00_marine_bulk_pipeline.sh
SAVE_BEDGRAPH=true  # set to false to skip bedgraph output

# DBSNP_BED and GENE_BED exported from 00_marine_bulk_pipeline.sh

mkdir -p $OUTPUTDIR

echo "---------------------"
echo "[LOG] Starting MARINE filter + normalize at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] Samplesheet : $SAMPLESHEET"
echo "[LOG] MARINE_DIR  : $MARINE_DIR"
echo "[LOG] FC_MATRIX   : $FC_MATRIX"
echo "[LOG] DBSNP_BED   : $DBSNP_BED"
echo "[LOG] GENE_BED    : $GENE_BED"
echo "[LOG] EDIT_TYPE   : $EDIT_TYPE"
echo "[LOG] MAX_FRAC    : $MAX_FRAC"
echo "---------------------"

module load mamba
mamba activate bioinformatics  # environment with pandas, pybedtools, matplotlib

while IFS=',' read -r sample r1 r2 libraryType; do
	[[ "$sample" == "sample" ]] && continue

	echo "[INFO] Processing $sample ..."

	# strandedness from 04_infer_strandedness_marine.sh
	STRANDEDNESS=$(cat $STRANDEDNESS_DIR/${sample}_strand_code.txt)

	# For strandedness=0 always use the unannotated file so the filter script
	# runs its own LOJ annotation (no -s flag) — required for the strand ==
	# feature_strand filter in the normalize step.
	# For strandedness 1/2 prefer the MARINE-annotated file; fall back to unannotated.
	if [[ "$STRANDEDNESS" -eq 0 ]]; then
		MARINE_INPUT="$MARINE_DIR/${sample}/final_filtered_site_info.tsv"
	else
		MARINE_INPUT="$MARINE_DIR/${sample}/final_filtered_site_info_annotated.tsv"
		if [[ ! -f "$MARINE_INPUT" ]]; then
			MARINE_INPUT="$MARINE_DIR/${sample}/final_filtered_site_info.tsv"
		fi
	fi

	SAMPLE_OUTDIR="$OUTPUTDIR/$sample"
	mkdir -p $SAMPLE_OUTDIR/filtered
	mkdir -p $SAMPLE_OUTDIR/normalized

	echo "[INFO] MARINE input : $MARINE_INPUT"

	# ── Step 1: Filter ────────────────────────────────────────────────────────
	echo "[INFO] Step 1: Filtering edits for $sample (strandedness=$STRANDEDNESS)"

	FILTER_ARGS=(
		--marine-results "$MARINE_INPUT"
		--strandedness   "$STRANDEDNESS"
		--dbsnp-bed      "$DBSNP_BED"
		--max-frac       "$MAX_FRAC"
		--output-dir     "$SAMPLE_OUTDIR/filtered"
	)

	if [[ "$STRANDEDNESS" -eq 0 ]]; then
		FILTER_ARGS+=(--annotation-bed "$GENE_BED")
	fi

	python helper_filter_edits_bulk.py "${FILTER_ARGS[@]}"

	echo "[DONE] Filter complete for $sample : $SAMPLE_OUTDIR/filtered/filtered_edits.tsv"

	# ── Step 2: Normalize ─────────────────────────────────────────────────────
	echo "[INFO] Step 2: Normalizing edits for $sample (edit_type=$EDIT_TYPE)"

	NORMALIZE_ARGS=(
		-i   "$SAMPLE_OUTDIR/filtered/filtered_edits.tsv"
		-c   "$FC_MATRIX"
		-s   "$sample"
		-str "$STRANDEDNESS"
		-e   "$EDIT_TYPE"
		-d   "$SAMPLE_OUTDIR/normalized"
	)
	[[ "$SAVE_BEDGRAPH" == false ]] && NORMALIZE_ARGS+=(--no-bedgraph)

	python helper_normalize_edits_bulk.py "${NORMALIZE_ARGS[@]}"

	echo "[DONE] Normalize complete for $sample : $SAMPLE_OUTDIR/normalized/"

done < "$SAMPLESHEET"

echo "[DONE] MARINE filter + normalize completed at $(date '+%Y-%m-%d %H:%M:%S')"
