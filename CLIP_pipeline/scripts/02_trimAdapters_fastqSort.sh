#!/bin/bash
#SBATCH --job-name=trimAdapters_fastqSort_02
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=72:00:00
#SBATCH --output=slurm_%u_%x_%j.log

# module load star
# module load samtools
# module load cutadapt
# module load fastq-tools

module load mamba
mamba activate clipper3

# cutadapt --version used here is 4.4 but eclip pipeline uses 1.14 which uses python 2.7. so create new environment 
# i think -f fastq is not available in 4.4 version of cutadapt

# HOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/ERBWTCLIP"
INPUTDIR="$HOMEDIR/01_UMI_clip"
OUTPUTDIR="$HOMEDIR"


mkdir -p $OUTPUTDIR/02_cutadapt1
mkdir -p $OUTPUTDIR/02_cutadapt2

cd $INPUTDIR

for fq in *.umi.fq;
do
        outfile=$(basename "$fq" .umi.fq)
	echo "Cutadapt for $fq which will be saved in as $OUTPUTDIR/02_cutadapt1/$outfile.umi.fqtr.fq"

	cutadapt --cores $SLURM_CPUS_PER_TASK \
                -O 1 \
                --match-read-wildcards \
                --times 1 \
                -e 0.1 \
                --quality-cutoff 6 \
                -m 18 \
		-o $OUTPUTDIR/02_cutadapt1/$outfile.umi.fqtr.fq \
		-a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                -a CTTCCGATCTACAAGTT \
                -a CTTCCGATCTTGGTCCT \
                -a AACTTGTAGATCGGA \
                -a AGGACCAAGATCGGA \
                -a ACTTGTAGATCGGAA \
                -a GGACCAAGATCGGAA \
                -a CTTGTAGATCGGAAG \
                -a GACCAAGATCGGAAG \
                -a TTGTAGATCGGAAGA \
                -a ACCAAGATCGGAAGA \
                -a TGTAGATCGGAAGAG \
                -a CCAAGATCGGAAGAG \
                -a GTAGATCGGAAGAGC \
                -a CAAGATCGGAAGAGC \
                -a TAGATCGGAAGAGCG \
                -a AAGATCGGAAGAGCG \
                -a AGATCGGAAGAGCGT \
                -a GATCGGAAGAGCGTC \
                -a ATCGGAAGAGCGTCG \
                -a TCGGAAGAGCGTCGT \
                -a CGGAAGAGCGTCGTG \
                -a GGAAGAGCGTCGTGT \
		$fq > $OUTPUTDIR/02_cutadapt1/$outfile.IP.umi.r1.fqTr.metrics

	echo "cutadapt1 done, cutadapt2 starting for $OUTPUTDIR/02_cutadapt2/$outfile.umi.fqtr.fq "

	cutadapt --cores $SLURM_CPUS_PER_TASK \
                -O 5 \
                --match-read-wildcards \
                --times 1 \
                -e 0.1 \
                --quality-cutoff 6 \
                -m 18 \
		-o $OUTPUTDIR/02_cutadapt2/$outfile.umi.fqtrtr.fq \
		-a AACTTGTAGATCGGA \
                -a AGGACCAAGATCGGA \
                -a ACTTGTAGATCGGAA \
                -a GGACCAAGATCGGAA \
                -a CTTGTAGATCGGAAG \
                -a GACCAAGATCGGAAG \
                -a TTGTAGATCGGAAGA \
                -a ACCAAGATCGGAAGA \
                -a TGTAGATCGGAAGAG \
                -a CCAAGATCGGAAGAG \
                -a GTAGATCGGAAGAGC \
                -a CAAGATCGGAAGAGC \
                -a TAGATCGGAAGAGCG \
                -a AAGATCGGAAGAGCG \
                -a AGATCGGAAGAGCGT \
                -a GATCGGAAGAGCGTC \
                -a ATCGGAAGAGCGTCG \
                -a TCGGAAGAGCGTCGT \
                -a CGGAAGAGCGTCGTG \
                -a GGAAGAGCGTCGTGT \
		$OUTPUTDIR/02_cutadapt1/$outfile.umi.fqtr.fq > $OUTPUTDIR/02_cutadapt2/$outfile.IP.umi.r1.fqTrTr.metrics
        
	echo "cutadapt2 done, starting sort"

	fastq-sort --id $OUTPUTDIR/02_cutadapt2/$outfile.umi.fqtrtr.fq > $OUTPUTDIR/02_cutadapt2/$outfile.umi.fqtrtr.sorted.fq

	echo "Fastq sort done for $outfile"
done
