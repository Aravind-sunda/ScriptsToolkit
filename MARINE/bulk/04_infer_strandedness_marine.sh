#!/bin/bash
#SBATCH --job-name=bulk_strandedness_marine_04
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=72:00:00
#SBATCH --output=slurm_%u_%x_%j.log

# SAMPLESHEET="/path/to/samplesheet.csv"  # columns: sample,r1,r2,libraryType (SE or PE)
# HOMEDIR="/path/to/project"

INPUTDIR="$HOMEDIR/03_star"
STRANDEDNESS_DIR="$HOMEDIR/04_strandedness"
MARINE_DIR="$HOMEDIR/04_marine"

# GENE_BED, GTF, FASTA exported from 00_marine_bulk_pipeline.sh

mkdir -p $STRANDEDNESS_DIR
mkdir -p $MARINE_DIR

echo "---------------------"
echo "[LOG] Starting strandedness inference and MARINE at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] Samplesheet : $SAMPLESHEET"
echo "---------------------"

# ── INFER EXPERIMENT (RSeQC) ──────────────────────────────────────────────────

module load mamba
mamba activate bioinformatics  # environment with RSeQC (infer_experiment.py)

echo "[LOG] Step 1: infer_experiment.py"

while IFS=',' read -r sample r1 r2 libraryType; do
	[[ "$sample" == "sample" ]] && continue

	BAM="$INPUTDIR/${sample}.Aligned.sortedByCoord.out.bam"

	echo "[INFO] infer_experiment.py for $sample"
	infer_experiment.py \
		-i $BAM \
		-r $GENE_BED \
		> $STRANDEDNESS_DIR/${sample}_strandedness.txt

	# Anchor on "This is" line, then grab +2 (forward) and +3 (reverse) fractions
	# Robust to leading blank lines in infer_experiment.py output
	# Threshold 0.75: above → stranded; below for both → unstranded
	STRAND_CODE=$(awk '/This is/{found=NR} found && NR==found+2{fwd=$NF} found && NR==found+3{rev=$NF} END{
		if (fwd+0 >= 0.75) print 1
		else if (rev+0 >= 0.75) print 2
		else print 0
	}' $STRANDEDNESS_DIR/${sample}_strandedness.txt)

	echo $STRAND_CODE > $STRANDEDNESS_DIR/${sample}_strand_code.txt

	echo "[DONE] $sample strandedness code: $STRAND_CODE (0=unstranded, 1=forward, 2=reverse)"

done < "$SAMPLESHEET"

# ── MARINE (STAMP edit calling) ───────────────────────────────────────────────

mamba deactivate
mamba activate marine_environment  # environment with MARINE; update env name as needed

echo "[LOG] Step 2: MARINE"

while IFS=',' read -r sample r1 r2 libraryType; do
	[[ "$sample" == "sample" ]] && continue

	BAM="$INPUTDIR/${sample}.Aligned.sortedByCoord.out.bam"
	MARINE_SAMPLE_DIR="$MARINE_DIR/$sample"
	mkdir -p $MARINE_SAMPLE_DIR

	STRAND_CODE=$(cat $STRANDEDNESS_DIR/${sample}_strand_code.txt)

	PE_FLAG=""
	[[ "$libraryType" == "PE" ]] && PE_FLAG="--paired_end"

	echo "[INFO] Running MARINE for $sample (strandedness: $STRAND_CODE, libraryType: $libraryType)"

	/home/tmhaxs421/brannanlab/tmhaxs421/MARINE_NextFlow/trial_pipeline/containers/MARINE-main/marine.py \
		--bam $BAM \
		--output_folder $MARINE_SAMPLE_DIR \
		--annotation_bedfile_path $GENE_BED \
		--strandedness $STRAND_CODE \
		--cores $SLURM_CPUS_PER_TASK \
		--min_base_quality 30 \
		--min_read_quality 255 \
		--sailor CT \
		--bedgraphs CT \
		$PE_FLAG
		# add or remove flags as needed
	echo "[DONE] MARINE complete for $sample : $MARINE_SAMPLE_DIR"

done < "$SAMPLESHEET"

echo "[DONE] Strandedness and MARINE completed at $(date '+%Y-%m-%d %H:%M:%S')"
