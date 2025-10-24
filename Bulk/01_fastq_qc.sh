#!/bin/bash

#SBATCH --job-name=bulk_fastq_qc
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_%u_%x_%j.log

echo "Starting job at $(date '+%Y-%m-%d %H:%M:%S')"

WORKINGDIR=""
INPUTDIR_FASTQ=""

## OUTPUT
OUTPUTDIR_QC="$WORKINGDIR/qc"


# =========================================================================================================================================
# QC AND STRANDEDNESS
# =========================================================================================================================================

module load mamba
mamba init bash
mamba activate
mamba activate bioinformatics

# Find all files that might be FASTQ
FASTQ_FILES=$(find "$INPUTDIR_FASTQ" -type f \
  \( -iname "*.fq" -o -iname "*.fastq" -o -iname "*.fq.gz" -o -iname "*.fastq.gz" \))

# Check what files were found
echo "$FASTQ_FILES"

mkdir -p "$OUTPUTDIR_QC/fastq_stats"
seqkit stats $FASTQ_FILES -T -a -j "$SLURM_CPUS_PER_TASK" --basename \
> "$OUTPUTDIR_QC/fastq_stats/fastq_stats.txt"

mkdir -p "$OUTPUTDIR_QC/fastqc"
fastqc -o "$OUTPUTDIR_QC/fastqc" -t "$SLURM_CPUS_PER_TASK" $FASTQ_FILES

mamba deactivate
mamba activate bioinformatics # running multi qc reports for all the stats and strandedness files
multiqc --force "$WORKINGDIR" -o "$OUTPUTDIR_QC/multiqc_reports" --filename "multiqc_report.html" --ignore ".*/" --fullnames