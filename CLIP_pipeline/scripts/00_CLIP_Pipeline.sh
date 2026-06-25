#!/bin/bash
#SBATCH --job-name=CLIP_Pipeline_00
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=72:00:00
#SBATCH --output=slurm_%u_%x_%j.log

set -euo pipefail

export HOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/delta_NRBD" # set homedir as the global variable. all scripts will use this variable
# export SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK # export the number of cpus per task as a global variable to be used in all scripts for parallel processing
    
if [ -z "$HOMEDIR" ]; then
    echo "HOMEDIR variable is not set. Please set it before running the script." # check if homedir exists or else fail script
    exit 1
fi

if [ ! -d "$HOMEDIR/data" ]; then
    echo "Data folder not found in HOMEDIR: $HOMEDIR/data" # check for data folder inside HOMEDIR
    exit 1
fi

echo "---------------------"
echo "[LOG] Starting job at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] Starting CLIP Pipeline"
echo "---------------------"
echo "[LOG] Step 1: Extract UMI"
echo "---------------------"
bash 01_extractUMI.sh

echo "---------------------"
echo "[LOG] Step 2: Trim Adapters and Sort Fastq"
echo "---------------------"
bash 02_trimAdapters_fastqSort.sh

echo "---------------------"
echo "[LOG] Step 3: Align to RepBase"
echo "---------------------"
bash 03_starAlignRepBase.sh

echo "---------------------"
echo "[LOG] Step 4: Align to Genome"
echo "---------------------"
bash 04_starAlignGenome.sh

echo "---------------------"
echo "[LOG] Step 5: Double Sort"
echo "---------------------"
bash 05_doubleSort.sh

echo "---------------------"
echo "[LOG] Step 6: Deduplicate and Sort"
echo "---------------------"
bash 06_dedup_sort.sh

echo "---------------------"
echo "[LOG] Step 7: Run Clipper"
echo "---------------------"
bash 07_clipper.sh

echo "---------------------"
echo "[LOG] Step 8: Make BigWig Files"
echo "---------------------"
bash 08_makeBigWigFiles.sh

echo "---------------------"
echo "[LOG] Step 9: Generate QC Read Count Summary"
echo "---------------------"
bash qc_read_counts.sh

echo "---------------------"
echo "[LOG] CLIP Pipeline Completed"
echo "[INFO] Run IDR Analysis Separately next"
echo "---------------------"