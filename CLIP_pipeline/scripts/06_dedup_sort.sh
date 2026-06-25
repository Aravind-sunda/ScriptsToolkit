#!/bin/bash
#SBATCH --job-name=dedup_sort_06
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=72:00:00
#SBATCH --output=slurm_%u_%x_%j.log

module load mamba
mamba activate
mamba activate clipper3
module load samtools

# HOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/ERBWTCLIP"

INPUTDIR="$HOMEDIR/05_double_sort"
OUTPUTDIR="$HOMEDIR/06_dedup_sort"

mkdir -p $OUTPUTDIR
cd $INPUTDIR

for fq in $INPUTDIR/*.genome-mappedSoSo.bam; do
	# outfile=$(echo  $fq | cut -f1 -d '.' )
	outfile=$(basename "$fq" .genome-mappedSoSo.bam)

	# index sorted bam
	echo "Indexing bam $INPUTDIR/$outfile.genome-mappedSoSo.bam"
	samtools index -@ $SLURM_CPUS_PER_TASK $fq # This indexes the input bam file. the outut bam is automatically indexed with .bai extension

	echo "dedup file is $fq"
	umi_tools dedup \
		--random-seed 1 \
		-I $fq \
		--method unique \
		--output-stats $OUTPUTDIR/$outfile.genome-mappedSoSo \
		-S $OUTPUTDIR/$outfile.genome-mappedSoSo.rmDup.bam \
		--log $OUTPUTDIR/$outfile.genome-mappedSoSo.dedup.log
		
	echo "Check Log file for dedup stats: $OUTPUTDIR/$outfile.genome-mappedSoSo.dedup.log"
	echo "Done dedup, now sorting"

	samtools sort -@ $SLURM_CPUS_PER_TASK \
		-o $OUTPUTDIR/$outfile.genome-mappedSoSo.rmDupSo.bam \
		$OUTPUTDIR/$outfile.genome-mappedSoSo.rmDup.bam
	
	# removing the intermediate file to save space. the final output is the sorted bam file with duplicates removed
	echo "Removing intermediate file $outfile.genome-mappedSoSo.rmDup.bam"
	rm $OUTPUTDIR/$outfile.genome-mappedSoSo.rmDup.bam

	echo "Done"
done

