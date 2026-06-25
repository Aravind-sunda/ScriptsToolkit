#!/bin/bash

#SBATCH --job-name=<job_name>
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --output=slurm_%u_%x_%j.log

module load mamba
module load star/2.7.10b
sds
GTF="/home/tmhaxs421/brannanlab/tmhaxs421/BRCA/reference/gencode.v46.basic.annotation.gtf"
STAR_INDEX="/home/tmhaxs421/brannanlab/VS_share/Genome_indexs/STAR_indexes/HS38_Gencode/index"
THREADS=$SLURM_CPUS_PER_TASK


READ_LENGTH=50
LIBRARY_TYPE="paired" # paired or single
OUTPUT_DIR="/home/tmhaxs421/brannanlab/tmhaxs421/BRCA/results/normal_vs_TNBC_MDA-MB-231"
TEMP_DIR="${OUTPUT_DIR}/temp"
VRL=FALSE
if [ "$VRL" = TRUE ]; then
    VARIABLE_READ_LENGTH="--variable-read-length"
else
    VARIABLE_READ_LENGTH=""     
fi

mkdir -p $OUTPUT_DIR
mkdir -p $TEMP_DIR

# -------------------------------
# Helper Functions
# -------------------------------

module load bioinformatics
seqkit stats --threads 36 /home/tmhaxs421/brannanlab/Vrutant/bulk-RNAseq/TERT_IR/rawdata/TERTIR/*.fastq.gz
# -------------------------------

# add an indentifying line to see what is running
echo ""
/home/tmhaxs421/brannanlab/tmhaxs421/applications/rmats_turbo_v4_3_0/run_rmats \
--s1  <experiment_csv> \
--s2  <control_csv> \
--gtf $GTF \
--bi  $STAR_INDEX \
-t $LIBRARY_TYPE --readLength $READ_LENGTH --nthread $THREADS $VARIABLE_READ_LENGTH \
--od  $OUTPUT_DIR \
--tmp $TEMP_DIR \
--novelSS \
--individual-counts