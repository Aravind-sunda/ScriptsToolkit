#!/bin/bash
#SBATCH --job-name=bulk_star_03
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=72:00:00
#SBATCH --output=slurm_%u_%x_%j.log

module load star/2.7.10b
module load samtools

# SAMPLESHEET="/path/to/samplesheet.csv"  # columns: sample,r1,r2,libraryType (SE or PE)
# HOMEDIR="/path/to/project"

INPUTDIR="$HOMEDIR/02_fastp"
OUTPUTDIR="$HOMEDIR/03_star"

# GENOMEDIR exported from 00_marine_bulk_pipeline.sh

mkdir -p $OUTPUTDIR

echo "---------------------"
echo "[LOG] Starting STAR alignment at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] Genome : $GENOMEDIR"
echo "[LOG] Output : $OUTPUTDIR"
echo "---------------------"

while IFS=',' read -r sample r1 r2 libraryType; do
	[[ "$sample" == "sample" ]] && continue

	echo "[INFO] Aligning $sample ($libraryType)"

	if [[ "$libraryType" == "SE" ]]; then
		READ_FILES_IN="$INPUTDIR/${sample}.trimmed.fastq.gz"
	elif [[ "$libraryType" == "PE" ]]; then
		READ_FILES_IN="$INPUTDIR/${sample}_R1.trimmed.fastq.gz $INPUTDIR/${sample}_R2.trimmed.fastq.gz"
	fi

	STAR \
		--alignEndsType EndToEnd \
		--genomeDir $GENOMEDIR \
		--genomeLoad NoSharedMemory \
		--outBAMcompression 10 \
		--outFileNamePrefix $OUTPUTDIR/${sample}. \
		--outFilterMultimapNmax 10 \
		--outFilterMultimapScoreRange 1 \
		--outFilterScoreMin 10 \
		--outReadsUnmapped Fastx \
		--outSAMattributes All \
		--outSAMmode Full \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--readFilesCommand zcat \
		--readFilesIn $READ_FILES_IN \
		--runMode alignReads \
		--runThreadN $SLURM_CPUS_PER_TASK

	BAM="$OUTPUTDIR/${sample}.Aligned.sortedByCoord.out.bam"

	echo "[LOG] Indexing BAM for $sample"
	samtools index -@ $SLURM_CPUS_PER_TASK $BAM

	echo "[DONE] Alignment complete for $sample : $BAM"

done < "$SAMPLESHEET"

echo "[DONE] STAR alignment completed at $(date '+%Y-%m-%d %H:%M:%S')"
