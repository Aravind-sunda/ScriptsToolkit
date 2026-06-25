#!/bin/bash
#SBATCH --job-name=bulk_qc_01
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%u_%x_%j.log

module load mamba
mamba activate bioinformatics  # environment with: fastqc, seqkit, multiqc

# SAMPLESHEET="/path/to/samplesheet.csv"  # columns: sample,r1,r2,libraryType (SE or PE)
# HOMEDIR="/path/to/project"

OUTPUTDIR="$HOMEDIR/01_qc"
FASTQC_DIR="$OUTPUTDIR/fastqc"
SEQKIT_DIR="$OUTPUTDIR/seqkit"

mkdir -p $FASTQC_DIR
mkdir -p $SEQKIT_DIR

echo "---------------------"
echo "[LOG] Starting bulk RNA-seq QC at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] Samplesheet : $SAMPLESHEET"
echo "[LOG] Output      : $OUTPUTDIR"
echo "---------------------"

# ── FASTQC ────────────────────────────────────────────────────────────────────

echo "[LOG] Step 1: FastQC"

while IFS=',' read -r sample r1 r2 libraryType; do
	[[ "$sample" == "sample" ]] && continue

	echo "[INFO] FastQC for $sample R1"
	fastqc \
		--threads $SLURM_CPUS_PER_TASK \
		--outdir $FASTQC_DIR \
		$r1

	if [[ "$libraryType" == "PE" ]]; then
		echo "[INFO] FastQC for $sample R2"
		fastqc \
			--threads $SLURM_CPUS_PER_TASK \
			--outdir $FASTQC_DIR \
			$r2
	fi

	echo "[DONE] FastQC complete for $sample"

done < "$SAMPLESHEET"

# ── SEQKIT STATS ──────────────────────────────────────────────────────────────

echo "[LOG] Step 2: seqkit stats"

while IFS=',' read -r sample r1 r2 libraryType; do
	[[ "$sample" == "sample" ]] && continue
	echo "$r1"
	[[ "$libraryType" == "PE" ]] && echo "$r2"
done < "$SAMPLESHEET" | xargs seqkit stats -a -j $SLURM_CPUS_PER_TASK -T \
	> $SEQKIT_DIR/seqkit_stats.tsv

echo "[DONE] seqkit stats written to $SEQKIT_DIR/seqkit_stats.tsv"

echo "[DONE] QC completed at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[INFO] MultiQC will run after fastp trimming (02_cutadapt.sh) to aggregate all reports"
