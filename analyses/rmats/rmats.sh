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

GTF="/home/tmhaxs421/brannanlab/tmhaxs421/BRCA/reference/gencode.v46.basic.annotation.gtf"
STAR_INDEX="/home/tmhaxs421/brannanlab/VS_share/Genome_indexs/STAR_indexes/HS38_Gencode/index"
THREADS=$SLURM_CPUS_PER_TASK


READ_LENGTH=""
LIBRARY_TYPE="" #paired,single
VRL= #TRUE,FALSE
STRANDEDNESS=""  #fr-unstranded,fr-firststrand,fr-secondstrand
# will go the libtype argument in rmats command

# Strandedness equivalents
# -----------------------------------------------------------------------------------------------------------------
# rMATS --libType   MARINE / featureCounts -s   RSeQC infer_experiment.py (PE)   Meaning
# -----------------------------------------------------------------------------------------------------------------
# fr-unstranded     0                           ~50% / ~50%                     no strand info
# fr-firststrand    2 (reverse)                 1+-,1-+,2++,2-- dominates       read1 antisense to transcript (dUTP, TruSeq Stranded, NEB Directional)
# fr-secondstrand   1 (forward)                 1++,1--,2+-,2-+ dominates       read1 sense to transcript (ligation-based, SMARTer/Ovation)
# -----------------------------------------------------------------------------------------------------------------
# Single-end RSeQC patterns: firststrand = +-,-+ ; secondstrand = ++,--
# Salmon: ISR = firststrand, ISF = secondstrand, IU = unstranded

FASTQ_DIR="path/to/fastq/directory"
EXPERIMENT_CSV="path/to/experiment.csv"
CONTROL_CSV="path/to/control.csv"
OUTPUT_DIR="path/to/output/directory"
TEMP_DIR="${OUTPUT_DIR}/temp"

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
module load mamba
mamba activate
mamba activate bioinformatics
seqkit stats --threads $THREADS $FASTQ_DIR/*.fastq.gz
# -------------------------------

mamba activate
module load star/2.7.10b
module load samtools/1.22.1

echo "[INFO] $(date)"

for DIR in "$OUTPUT_DIR" "$TEMP_DIR"; do
    if [ -d "$DIR" ] && [ "$(ls -A "$DIR")" ]; then
        echo "[INFO] Cleaning up previous rMATS files in $DIR"
        rm -rf "$DIR"/*
    fi
done

echo "[INFO] Running rMATS"
/home/tmhaxs421/brannanlab/tmhaxs421/applications/rmats_turbo_v4_3_0/run_rmats \
--s1  $EXPERIMENT_CSV \
--s2  $CONTROL_CSV \
--gtf $GTF \
--bi  $STAR_INDEX \
-t $LIBRARY_TYPE --readLength $READ_LENGTH --nthread $THREADS $VARIABLE_READ_LENGTH \
--libType $STRANDEDNESS \
--od  $OUTPUT_DIR \
--tmp $TEMP_DIR \
--novelSS \
--individual-counts