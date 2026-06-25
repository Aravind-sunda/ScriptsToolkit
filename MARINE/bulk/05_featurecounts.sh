#!/bin/bash
#SBATCH --job-name=bulk_featurecounts_05
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%u_%x_%j.log

module load mamba
mamba activate bioinformatics  # environment with subread/featureCounts v1.5.3

# SAMPLESHEET="/path/to/samplesheet.csv"  # columns: sample,r1,r2,libraryType (SE or PE)
# HOMEDIR="/path/to/project"

INPUTDIR="$HOMEDIR/03_star"
STRANDEDNESS_DIR="$HOMEDIR/04_strandedness"  # strand_code files written by 04_infer_strandedness_marine.sh
OUTPUTDIR="$HOMEDIR/05_featurecounts"

# GTF exported from 00_marine_bulk_pipeline.sh

mkdir -p $OUTPUTDIR

echo "---------------------"
echo "[LOG] Starting featureCounts at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] GTF    : $GTF"
echo "[LOG] Output : $OUTPUTDIR"
echo "---------------------"

while IFS=',' read -r sample r1 r2 libraryType; do
	[[ "$sample" == "sample" ]] && continue

	BAM="$INPUTDIR/${sample}.Aligned.sortedByCoord.out.bam"
	OUTFILE="$OUTPUTDIR/${sample}.featurecounts.txt"

	STRAND_CODE_FILE="$STRANDEDNESS_DIR/${sample}_strand_code.txt"
	if [[ ! -f "$STRAND_CODE_FILE" ]]; then
		echo "[ERROR] Strand code file not found for $sample: $STRAND_CODE_FILE"
		echo "[ERROR] Run 04_infer_strandedness_marine.sh before this script"
		exit 1
	fi
	STRANDEDNESS=$(cat $STRAND_CODE_FILE)

	echo "[INFO] featureCounts for $sample ($libraryType, strandedness=$STRANDEDNESS)"

	if [[ "$libraryType" == "SE" ]]; then

		featureCounts \
			-T $SLURM_CPUS_PER_TASK \
			-a $GTF \
			-t exon \
			-g gene_name \
			-s $STRANDEDNESS \
			-o $OUTFILE \
			$BAM

	elif [[ "$libraryType" == "PE" ]]; then

		featureCounts \
			-T $SLURM_CPUS_PER_TASK \
			-a $GTF \
			-t exon \
			-g gene_name \
			-s $STRANDEDNESS \
			-p --countReadPairs \
			-o $OUTFILE \
			$BAM

	fi

	echo "[DONE] featureCounts complete for $sample : $OUTFILE"

done < "$SAMPLESHEET"

echo "[DONE] featureCounts completed at $(date '+%Y-%m-%d %H:%M:%S')"

echo "---------------------"
echo "[LOG] Merging per-sample counts into matrix"
echo "---------------------"

python helper_merge_counts.py \
	--indir $OUTPUTDIR \
	--outfile $OUTPUTDIR/counts_matrix_combined.tsv
	# --drop-annotation

echo "[DONE] Count matrix written to $OUTPUTDIR/counts_matrix_combined.tsv"
