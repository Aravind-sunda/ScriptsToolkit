#!/bin/bash
#SBATCH --job-name=bulk_fastp_02
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
mamba activate bioinformatics  # environment with fastp

# SAMPLESHEET="/path/to/samplesheet.csv"  # columns: sample,r1,r2,libraryType (SE or PE)
# HOMEDIR="/path/to/project"

OUTPUTDIR="$HOMEDIR/02_fastp"
QC_DIR="$OUTPUTDIR/qc"

mkdir -p $OUTPUTDIR
mkdir -p $QC_DIR

echo "---------------------"
echo "[LOG] Starting fastp trimming at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] Samplesheet : $SAMPLESHEET"
echo "[LOG] Output      : $OUTPUTDIR"
echo "---------------------"

while IFS=',' read -r sample r1 r2 libraryType; do
	[[ "$sample" == "sample" ]] && continue

	echo "[INFO] fastp for $sample ($libraryType)"

	if [[ "$libraryType" == "SE" ]]; then

		fastp \
			--in1 $r1 \
			--out1 $OUTPUTDIR/${sample}.trimmed.fastq.gz \
			--thread $SLURM_CPUS_PER_TASK \
			--qualified_quality_phred 6 \
			--length_required 20 \
			--json $QC_DIR/${sample}.fastp.json \
			--html $QC_DIR/${sample}.fastp.html

		echo "[DONE] SE trimming complete for $sample"

	elif [[ "$libraryType" == "PE" ]]; then

		fastp \
			--in1 $r1 \
			--in2 $r2 \
			--out1 $OUTPUTDIR/${sample}_R1.trimmed.fastq.gz \
			--out2 $OUTPUTDIR/${sample}_R2.trimmed.fastq.gz \
			--detect_adapter_for_pe \
			--thread $SLURM_CPUS_PER_TASK \
			--qualified_quality_phred 6 \
			--length_required 20 \
			--json $QC_DIR/${sample}.fastp.json \
			--html $QC_DIR/${sample}.fastp.html

		echo "[DONE] PE trimming complete for $sample"

	fi

done < "$SAMPLESHEET"

echo "[LOG] Running MultiQC on all logs in $HOMEDIR"
multiqc $HOMEDIR \
	--outdir $HOMEDIR/multiqc \
	--filename multiqc_report \
	--force

echo "[DONE] fastp trimming completed at $(date '+%Y-%m-%d %H:%M:%S')"
