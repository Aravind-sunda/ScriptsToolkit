#!/bin/bash
#SBATCH --job-name=doubleSort_05
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

# HOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/ERBWTCLIP"

INPUTDIR="$HOMEDIR/04_star_hg38_rep_sam"
OUTPUTDIR="$HOMEDIR/05_double_sort"

mkdir -p $OUTPUTDIR

# 294HATertNRBDdeltaIp1_S4_R1_001.repbaseMapped.genomeMapped.Aligned.out.bam

for bam in $INPUTDIR/*.repbaseMapped.genomeMapped.Aligned.out.bam; do
        # outfile=$(echo  $bam | cut -f1 -d '.' )
        outfile=$(basename "$bam" .repbaseMapped.genomeMapped.Aligned.out.bam)
        
        echo "Sorting bam $bam"
        samtools \
                sort \
                -@ 36 \
                -n \
                -o $OUTPUTDIR/$outfile.genome-mappedSo.bam \
                $bam

        samtools \
                sort \
                -@ 36 \
                -o $OUTPUTDIR/$outfile.genome-mappedSoSo.bam \
                $OUTPUTDIR/$outfile.genome-mappedSo.bam
        
        # removing one of the sort files to save space
        echo "Removing intermediate file $outfile.genome-mappedSo.bam"
        rm $OUTPUTDIR/$outfile.genome-mappedSo.bam

        echo "Done"
done

