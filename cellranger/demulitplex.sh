#!/bin/bash
#SBATCH --partition=defq
#SBATCH --job-name=demux_cc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=96:00:00
#SBATCH --output=slurm_%u_%x_%j.log

echo "Starting job at $(date '+%Y-%m-%d %H:%M:%S')"

module load cellranger/9.0.0
# Note: remove the lane column from the sample sheet if you want to demultiplex all lanes together


RUNFOLDER=""
SAMPLESHEET=""
OUTPUTDIR=""
# STATSDIR=""
# REPORTSDIR=""

mkdir -p "${OUTPUTDIR}"
cellranger mkfastq \
    --run "${RUNFOLDER}" \
    --csv "${SAMPLESHEET}" \
    --output-dir "${OUTPUTDIR}" \
    --barcode-mismatches 2 \
    --no-lane-splitting \
    --loading-threads 4 \
    --processing-threads 28 \
    --writing-threads 4 




# bcl2fastq \
#     --runfolder-dir "${RUNFOLDER}" \
#     --barcode-mismatches 1 \
#     --sample-sheet "${SAMPLESHEET}" \
#     --output-dir "${OUTPUTDIR}" \
#     --stats-dir "${STATSDIR}" \
#     --reports-dir "${REPORTSDIR}" \
#     --no-lane-splitting \
#     --loading-threads 4 \
#     --processing-threads 28 \
#     --writing-threads 4 