#!/bin/bash
#SBATCH --job-name=extractUMI_01
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
# module load eclip
mamba activate
mamba activate clipper3


# HOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/ERBWTCLIP"

INPUTDIR="$HOMEDIR/data"
OUTPUTDIR="$HOMEDIR/01_UMI_clip"

mkdir -p "$OUTPUTDIR"

for fq in "$INPUTDIR"/*.fastq.gz; 
do
	outfile=$(basename "$fq" .fastq.gz)
	
	echo "UMI cleanup started for $outfile"
	umi_tools extract \
		--random-seed 1 \
		--bc-pattern NNNNNNNNNN \
		--log "$OUTPUTDIR/$outfile.log.---.--.metrics" \
		--stdin "$fq" \
		--stdout "$OUTPUTDIR/$outfile.umi.fq"
	echo "Done"
done