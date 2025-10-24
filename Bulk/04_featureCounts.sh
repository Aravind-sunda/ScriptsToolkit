#!/bin/bash
#SBATCH --partition=defq
#SBATCH --job-name=featureCounts
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

# REFERENCES
GTF=""

# SUFFIXES
BAM_suffix="Aligned.sortedByCoord.out.bam"


module load mamba
mamba activate
mamba activate bioinformatics

# ========================================================================
# RUNNING FEATURECOUNTS
mkdir -p $WORKING_DIR/featureCounts

featureCounts -T $SLURM_CPUS_PER_TASK \
    -p --countReadPairs \
    -t exon -g gene_id  \
    --extraAttributes "gene_name" \
    -a $GTF \
    -s 0 --primary \
    -o $WORKING_DIR/featureCounts/featureCounts_trimmed.exon.counts.txt \
    $WORKING_DIR/star_output/*${BAM_suffix}
# ========================================================================


# featureCounts options
# -s  {0,1,2} #(0 = unstranded, 1 = stranded, 2 = reversely stranded)
# -p  --countReadPairs # paired end reads 
# --primary # only use primary alignments and multi-mapping reads will be considered automatically
# -t {gene,exon,CDS}
# -g {gene_id,gene_name,transcript_id}
# -o take reads which overlap with multiple features in the gtf file, use for unstranded data

# EXTRA PYTHON HELPER FUNCTION
# PYTHON FUNCTION TO PROCESS FEATURECOUNTS OUTPUT
