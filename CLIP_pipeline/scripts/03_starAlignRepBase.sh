#!/bin/bash
#SBATCH --job-name=starAlignRepBase_03
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=72:00:00
#SBATCH --output=slurm_%u_%x_%j.log

module load star

# HOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/ERBWTCLIP"

INPUTDIR="$HOMEDIR/02_cutadapt2"
OUTPUTDIR="$HOMEDIR/03_repbase_alignment"

GENOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/CLIP_pipeline/reference/Encode/repbase_index/index1"

mkdir -p $OUTPUTDIR


for fq in $INPUTDIR/*.umi.fqtrtr.sorted.fq; do
	# outfile=$(echo  $fq | cut -f1 -d '.' )
	outfile=$(basename "$fq" .umi.fqtrtr.sorted.fq)

	echo "Starting alignment for $fq."
	echo "Your output will be in $OUTPUTDIR/$outfile"

	STAR \
	--alignEndsType EndToEnd \
	--genomeDir $GENOMEDIR \
	--genomeLoad NoSharedMemory \
	--outBAMcompression 10 \
	--outFileNamePrefix $OUTPUTDIR/$outfile.repbaseMapped. \
	--outFilterMultimapNmax 30 \
	--outFilterMultimapScoreRange 1 \
	--outFilterScoreMin 10 \
	--outFilterType BySJout \
	--outReadsUnmapped Fastx \
	--outSAMattrRGline ID:foo \
	--outSAMattributes All \
	--outSAMmode Full \
	--outSAMtype BAM Unsorted \
	--outSAMunmapped Within \
	--outStd Log \
	--runMode alignReads \
	--runThreadN $SLURM_CPUS_PER_TASK \
	--readFilesIn $fq
	
	echo "Alignment is done"
done

