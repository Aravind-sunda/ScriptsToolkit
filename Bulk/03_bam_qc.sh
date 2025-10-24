#!/bin/bash

#SBATCH --partition=defq
#SBATCH --job-name=stamp_bulk_qc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=96:00:00
#SBATCH --output=slurm_%u_%x_%j.log

echo "Starting job at $(date '+%Y-%m-%d %H:%M:%S')"

## INPUTS AND OUTPUTS
WORKING_DIR="/home/tmhaxs421/brannanlab/tmhaxs421/riboSTAMP_mouse" # Remove
INPUT_BAM_FILE_DIR="$WORKING_DIR/results/star_output" # Remove bam files in end and add /<path_to_bam_files>

OUTPUT_STRANDEDNESS_DIR="$WORKING_DIR/qc/strandedness"
OUTPUT_STATS_DIR="$WORKING_DIR/qc/bam_stats"
OUTPUT_QC_DIR="$WORKING_DIR/qc"


## REFERENCES
REFERENCE_PATH="" # adding the reference for alignment 
BED_PATH="" # for annotation and reseqc strandedness
GTF_PATH=""

PREFIX=".Aligned.sortedByCoord.out.bam" # remember to add the dot to the prefix 
# 1-Brain.Aligned.sortedByCoord.out.bam


MAPQ_THRESHOLD="5" # only Unique Aligners

# =========================================================================================================================================
# PART 1: QC AND STRANDEDNESS
# =========================================================================================================================================

module load mamba
mamba init bash
mamba activate
mamba activate bioinformatics

mkdir -p "$OUTPUT_STRANDEDNESS_DIR"
for bam_file in "$INPUT_BAM_FILE_DIR"/*"$PREFIX"; do # running RSEQC to determine strandedness
    sample=$(basename "$bam_file" "$PREFIX")
    echo "[INFO] Processing $sample (file: $(basename "$bam_file"))"

    echo "Running RSEQC for $sample"
    infer_experiment.py \
        -i "$bam_file" \
        -r "$BED_PATH" \
        -s 1000000 > "$OUTPUT_STRANDEDNESS_DIR/${sample}_strandedness.txt"

    echo "[$(date)] Done processing sample $sample"
done

mamba deactivate 
mamba activate gatk
for bam_file in "$INPUT_BAM_FILE_DIR"/*"$PREFIX"; do # Running GATK base quality score distribution to see how the 
    sample=$(basename "$bam_file" "$PREFIX")
    echo "[INFO] Processing $sample (file: $(basename "$bam_file"))"

    gatk QualityScoreDistribution \
        -I "$bam_file" \
        -O "$OUTPUT_STATS_DIR/$sample/${sample}_quality_score_distribution.tsv" \
        -CHART "$OUTPUT_STATS_DIR/$sample/${sample}_quality_score_distribution_chart.pdf" \
        --ALIGNED_READS_ONLY
done


#==========================================================================================================================================
# utlities
# # check strandedness of the data
# # Convert GTF to genePred
# /home/tmhaxs421/brannanlab/tmhaxs421/applications/ncbi_utilities/gtfToGenePred -genePredExt /home/tmhaxs421/brannanlab/VS_share/Genome_indexs/tdtmouse/gencode.vM37.tdTomato.annotation.gtf \
#     /home/tmhaxs421/brannanlab/VS_share/Genome_indexs/tdtmouse/annotation.genePred

# # Convert to BED12
# /home/tmhaxs421/brannanlab/tmhaxs421/applications/ncbi_utilities/genePredToBed /home/tmhaxs421/brannanlab/VS_share/Genome_indexs/tdtmouse/annotation.genePred \
#     /home/tmhaxs421/brannanlab/VS_share/Genome_indexs/tdtmouse/refgene.bed