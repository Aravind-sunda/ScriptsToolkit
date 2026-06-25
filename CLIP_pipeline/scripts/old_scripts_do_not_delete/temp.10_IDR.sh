#!/bin/bash

#SBATCH --partition=defq
#SBATCH --job-name=IDR_09
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%u_%x_%j.log


module load mamba
mamba init bash
module load shared
module load mergepeaks/0.1.0
module load cwltool/3.1.20250925164626 
mamba activate bioinformatics

# set -euo pipefail # this causes an internal error which does not change the scripts
# ============================================================================================================================================================
# DEFININING THE INPUT FILES
rep1_clip_bam=""
rep2_clip_bam=""
rep1_input_bam=""
rep2_input_bam=""

# Only peaks of IP files
rep1_peaks_bed=""
rep2_peaks_bed=""
# Format of peaks bed file
# chr1	15306	15338	ENSG00000227232.5_0_23	1.543195658818628e-08	-	15320	15324

chrom_sizes="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/TERTNRBDscripts/reference/Encode/hg38_chrom_sizes.txt"


HOMEDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_Tert_NRBD/IDR/Delta_293_Tert_NRBD"
WORKDIR="$HOMEDIR/IDR"
# ============================================================================================================================================================

# ============================================================================================================================================================
# STARTING THE SCRIPT
mkdir -p $WORKDIR
cd $WORKDIR
# ============================================================================================================================================================
module load samtools/1.22.1

samtools view -c -F 4 $rep1_clip_bam > rep1_clip.readnum
samtools view -c -F 4 $rep2_clip_bam > rep2_clip.readnum
samtools view -c -F 4 $rep1_input_bam > rep1_input.readnum
samtools view -c -F 4 $rep2_input_bam > rep2_input.readnum


/home/tmhaxs421/brannanlab/tmhaxs421/applications/eCLIP-master/bin/overlap_peakfi_with_bam.pl \
    $rep1_clip_bam \
    $rep1_input_bam \
    $rep1_peaks_bed \
    rep1_clip.readnum \
    rep1_input.readnum \
    rep1_normed_peaks.bed 2>/dev/null

echo "Normalize over input for rep1 done"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/eCLIP-master/bin/overlap_peakfi_with_bam.pl \
    $rep2_clip_bam \
    $rep2_input_bam \
    $rep2_peaks_bed \
    rep2_clip.readnum \
    rep2_input.readnum \
    rep2_normed_peaks.bed 2>/dev/null

echo "Normalize over input for rep2 done"

# ============================================================================================================================================================
/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl \
    rep1_normed_peaks.bed.full \
    rep1_normed_peaks.compressed.bed \
    rep1_normed_peaks.compressed.bed.full 

echo "Compress rep1 peaks done"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl \
    rep2_normed_peaks.bed.full \
    rep2_normed_peaks.compressed.bed \
    rep2_normed_peaks.compressed.bed.full

echo "Compress rep2 peaks done"

# ============================================================================================================================================================
echo "Computing Entropy for rep1 and rep2 peaks"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/make_informationcontent_from_peaks.pl \
    rep1_normed_peaks.compressed.bed.full \
    rep1_clip.readnum \
    rep1_input.readnum \
    rep1_normed_peaks.compressed.bed.entropy.full \
    rep1_normed_peaks.compressed.bed.entropy.excessreads 

echo " Computing Entropy for rep1 peaks done"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/make_informationcontent_from_peaks.pl \
    rep2_normed_peaks.compressed.bed.full \
    rep2_clip.readnum \
    rep2_input.readnum \
    rep2_normed_peaks.compressed.bed.entropy.full \
    rep2_normed_peaks.compressed.bed.entropy.excessreads 

echo "Computing Entropy for rep2 peaks done"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/full_to_bed.py \
    --input rep1_normed_peaks.compressed.bed.entropy.full \
    --output rep1_normed_peaks.compressed.bed.entropy.bed

echo "Convert rep1 peaks to bed done"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/full_to_bed.py \
    --input rep2_normed_peaks.compressed.bed.entropy.full \
    --output rep2_normed_peaks.compressed.bed.entropy.bed

echo "Convert rep2 peaks to bed done"

# ============================================================================================================================================================
# this idr environment must be installed carefully
# create the environment with the following command and update the script with the 2.0.3 application
# mamba create -n idr -c conda-forge -c bioconda idr=2.0.2
# mamba activate idr
# python3 /home/tmhaxs421/brannanlab/tmhaxs421/applications/idr-2.0.3/setup.py install
# removes shared modules since it has a more recent version of numpy it will conflict with the idr version of numpy

module purge
module load mamba
mamba init bash
mamba activate idr

echo "Running IDR analysis"


awk 'BEGIN{OFS="\t"}{print $0, 0, 0, 0}' rep1_normed_peaks.compressed.bed.entropy.bed > rep1_normed_peaks.compressed.bed.entropy.bed.9col
awk 'BEGIN{OFS="\t"}{print $0, 0, 0, 0}' rep2_normed_peaks.compressed.bed.entropy.bed > rep2_normed_peaks.compressed.bed.entropy.bed.9col

# python3 /home/tmhaxs421/brannanlab/tmhaxs421/applications/idr-2.0.2/idr/idr.py \
idr \
    --samples rep1_normed_peaks.compressed.bed.entropy.bed.9col rep2_normed_peaks.compressed.bed.entropy.bed.9col \
    --input-file-type bed \
    --rank 5 \
    --peak-merge-method max \
    --plot \
    -o 01v02.idr.out

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/parse_idr_peaks.pl \
    01v02.idr.out \
    rep1_normed_peaks.compressed.bed.entropy.full \
    rep2_normed_peaks.compressed.bed.entropy.full \
    01v02.idr.out.bed
    
echo "IDR analysis done"
# ============================================================================================================================================================
module purge
mamba deactivate
module load mamba
mamba init bash
module load shared
module load mergepeaks/0.1.0
module load cwltool/3.1.20250925164626 
mamba activate bioinformatics


/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/overlap_peakfi_with_bam.pl \
    $rep1_clip_bam \
    $rep1_input_bam \
    01v02.idr.out.bed \
    rep1_clip.readnum \
    rep1_input.readnum \
    01v02.IDR.out.0102merged.01.bed 2>/dev/null

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/overlap_peakfi_with_bam.pl \
    $rep2_clip_bam \
    $rep2_input_bam \
    01v02.idr.out.bed \
    rep2_clip.readnum \
    rep2_input.readnum \
    01v02.IDR.out.0102merged.02.bed 2>/dev/null

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/get_reproducing_peaks.pl \
    01v02.IDR.out.0102merged.01.bed.full \
    01v02.IDR.out.0102merged.02.bed.full \
    reproducible_peaks.01.bed.full \
    reproducible_peaks.02.bed.full \
    reproducible_peaks.bed \
    reproducible_peaks.custombed \
    rep1_normed_peaks.compressed.bed.entropy.full \
    rep2_normed_peaks.compressed.bed.entropy.full \
    01v02.idr.out


echo "Get reproducible peaks done"
echo "IDR analysis run done"