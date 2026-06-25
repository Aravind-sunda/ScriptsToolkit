#!/bin/bash
#SBATCH --job-name=makeBigWig_08
#SBATCH --nodes=1
#SBATCH --partition=bigmemq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=500G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=72:00:00
#SBATCH --output=slurm_%u_%x_%j.log

# module load mamba

module load apptainer
module load samtools

APPTAINER_CONTAINER="/home/tmhaxs421/brannanlab/tmhaxs421/applications/docker_container/makebigwigfiles_0.0.3.sif"
condo="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/"
home="/home/tmhaxs421"

# HOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/ERBWTCLIP"
INPUTDIR="$HOMEDIR/06_dedup_sort"
OUTPUTDIR="$HOMEDIR/08_bigwig"
ENCODE="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/CLIP_pipeline/reference/Encode"

mkdir -p $OUTPUTDIR

for bam in $INPUTDIR/*.genome-mappedSoSo.rmDupSo.bam; do
	# outfile=$(echo  $bam | cut -f1 -d '.' )
	outfile=$(basename "$bam" .bam)
	
	echo "Indexing bam $bam"
	samtools index -@ $SLURM_CPUS_PER_TASK $bam
	echo "Done indexing, now bigwigfile"

	apptainer run \
	--bind /condo/,$condo,$home,$HOMEDIR,$INPUTDIR,$OUTPUTDIR,$ENCODE \
	$APPTAINER_CONTAINER \
	makebigwigfiles \
		--bw_pos $OUTPUTDIR/$outfile.norm.pos.bw \
		--bw_neg $OUTPUTDIR/$outfile.norm.neg.bw \
		--bam $bam \
		--genome $ENCODE/hg38_chrom_sizes.txt \
		--direction f
	echo "Done Making BigWig files"

	
done


module load mamba
module load deeptools

INPUTDIR="$HOMEDIR/06_dedup_sort"
OUTPUTDIR="$HOMEDIR/08_bigwig_combined"
ENCODE="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/CLIP_pipeline/reference/Encode/"

mkdir -p $OUTPUTDIR

for bam in $INPUTDIR/*.genome-mappedSoSo.rmDupSo.bam; do
	outfile=$(basename "$bam" .bam)
	echo "Bedgraphing $bam"
	bamCoverage -b $bam -o $OUTPUTDIR/$outfile.bw -p $SLURM_CPUS_PER_TASK 
	echo "Done"
done


INPUTDIR="$HOMEDIR/06_dedup_sort"
OUTPUTDIR="$HOMEDIR/08_bigwig_combined_normalized"
ENCODE="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/CLIP_pipeline/reference/Encode/"

mkdir -p $OUTPUTDIR
for bam in $INPUTDIR/*.genome-mappedSoSo.rmDupSo.bam; do
	outfile=$(basename "$bam" .bam)
	echo "Bedgraphing $bam"
	bamCoverage -b $bam --normalizeUsing RPKM  -p $SLURM_CPUS_PER_TASK --effectiveGenomeSize 2805636231 -o $OUTPUTDIR/$outfile.normalized.bw
	echo "Done"
done
