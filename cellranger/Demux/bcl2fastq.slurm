#!/bin/bash -l

#SBATCH --job-name=bcltofastq
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --output=slurm_%u_%x_%j.log


## Load modules, add version numbers
module load bcl2fastq/2.20

RUNFOLDER=""
SAMPLESHEET=""
OUTPUTDIR=""

# Cellranger mkfastq is deprecasted use bclconvert. Bclconvert is planned to replace bcl2fastq

bcl2fastq \
	--runfolder-dir "${RUNFOLDER}" \
    --sample-sheet "${SAMPLESHEET}" \
	--output-dir "${OUTPUTDIR}" \
	--barcode-mismatches 2 \
	--no-lane-splitting \
	--loading-threads 4 \
	--processing-threads 28 \
	--writing-threads 4 

