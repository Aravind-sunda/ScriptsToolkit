#!/bin/bash
#SBATCH --partition=defq
#SBATCH --job-name=star
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=96:00:00
#SBATCH --output=slurm_%u_%x_%j.log

echo "Starting job at $(date '+%Y-%m-%d %H:%M:%S')"

WORKING_DIR=""
INPUTDIR_FASTQ=""

# REFERENCES
GENOME_STAR_INDEX=""
GTF=""

# SUFFIXES
FQ_SUFFIX=".fastq.gz"
BAM_suffix="Aligned.sortedByCoord.out.bam"

# MKDIR
mkdir -p $WORKING_DIR/star_output

# MODULES
module load star


for fq1 in $INPUTDIR_FASTQ/*_R1_*${FQ_SUFFIX}; do

    # get basename for sample
    sample=$(basename "$fq1" "$FQ_SUFFIX" | sed 's/_R1_.*//')

    # find mate (R2)
    fq2=${fq1/_R1_/_R2_}

    if [[ ! -f "$fq2" ]]; then
        echo "Warning: Mate file not found for $fq1"
        continue
    fi

    echo "Started to align paired-end reads for sample $sample"

    STAR --genomeDir "$GENOMEDIR" \
         --sjdbGTFfile "$GTF" \
         --runThreadN "$SLURM_CPUS_PER_TASK" \
         --readFilesIn "$fq1" "$fq2" \
         --readFilesCommand zcat \
         --sampleNamePrefix "$WORKING_DIR/star_output/${sample}." \
         --outSAMattributes All \
         --outSAMtype BAM SortedByCoordinate \
         --outReadsUnmapped Fastx \
         --outSAMmode Full \
         --quantMode GeneCounts TranscriptomeSAM \
         --outFilterMultimapNmax 10 \
         --outFilterMultimapScoreRange 1 \
         --outFilterScoreMin 10 \
         --alignEndsType EndToEnd 
    
done

# STAR OPTIONAL PARAMETERS:
# --outFilterScoreMinOverLread 0.1 # Used for shorter reads to set minimum score relative to read length
# --outFilterMatchNminOverLread 0.1 # Used for shorter reads to set minimum match score relative to read length
# --alignEndsType Local # used for soft clipping of adapters
# --alignEndsType EndToEnd # no soft clipping of adapters

# ========================================================================
# INDEX BAM FILES
module load mamba
mamba activate
mamba activate bioinformatics

for bam in $WORKING_DIR/star_output/*${BAM_suffix}; do
    sample=$(basename "$bam" ".$BAM_suffix")
    echo "Indexing BAM file for sample $sample"
    samtools index -@ "$SLURM_CPUS_PER_TASK" "$bam"
    echo "Done with $sample"
done


