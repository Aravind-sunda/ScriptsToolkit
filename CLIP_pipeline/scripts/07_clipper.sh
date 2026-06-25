#!/bin/bash
#SBATCH --job-name=clipper_07
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=72:00:00
#SBATCH --output=slurm_%u_%x_%j.log

module load samtools
module load apptainer
echo "which modules are loaded?"
module list
module load samtools
# # Specify your Apptainer container
# APPTAINER_CONTAINER="/cm/shared/apptainer-images/clipper_5d865bb.sif"
# apptainer run --bind /condo/  $APPTAINER_CONTAINER /bin/bash

APPTAINER_CONTAINER="/cm/shared/apptainer-images/clipper_5d865bb.sif"
condo="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/"
home="/home/tmhaxs421"

# HOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/ERBWTCLIP"
INPUTDIR="$HOMEDIR/06_dedup_sort"
OUTPUTDIR="$HOMEDIR/07_clipper_output"

mkdir -p $OUTPUTDIR

for bam in $INPUTDIR/*.genome-mappedSoSo.rmDupSo.bam; do
	echo "your file is $bam"
    # outfile=$(echo $bam | cut -f1 -d '.' )
	outfile=$(basename "$bam" .bam)
	
	echo "Indexing bam $bam"
	samtools index -@ $SLURM_CPUS_PER_TASK $bam
	
	echo "Done indexing, now running clipper"
	apptainer run \
	--bind /condo/,$condo,$home,$HOMEDIR,$INPUTDIR,$OUTPUTDIR \
	$APPTAINER_CONTAINER \
	clipper \
		--species GRCh38_v29e \
		--bam $bam \
		--processors=$SLURM_CPUS_PER_TASK \
		--outfile $OUTPUTDIR/$outfile.peakClusters.bed
done