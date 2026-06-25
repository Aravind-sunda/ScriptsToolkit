if [ "$#" -ne 7 ]; then
  echo "Error: Usage: $0 REP1_CLIP_BAM REP2_CLIP_BAM SMI_input_bam REP1_PEAKS_BED REP2_PEAKS_BED CHROM_SIZES WORKDIR"
  exit 1
fi

rep1_clip_bam=$1
rep2_clip_bam=$2
SMI_input_bam=$3
rep1_peaks_bed=$4
rep2_peaks_bed=$5
chrom_sizes=$6
WORKDIR=$7

module load mamba
mamba init bash
module load shared
module load mergepeaks/0.1.0
module load cwltool/3.1.2 
mamba activate bioinformatics


# ============================================================================================================================================================
# DEFININING THE INPUT FILES
# rep1_clip_bam="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_HA_tert_fragments/dedup_sort/1to350Ip1_S1_R1_002.genome-mappedSoSo.rmDupSo.bam"
# rep2_clip_bam="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_HA_tert_fragments/dedup_sort/1to350Ip2_S2_R1_002.genome-mappedSoSo.rmDupSo.bam"

# SMI_input_bam=

# rep1_peaks_bed="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_HA_tert_fragments/Clipper_output/1to350Ip1_S1_R1_002.rerun.peaks.bed"
# rep2_peaks_bed="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_HA_tert_fragments/Clipper_output/1to350Ip2_S2_R1_002.rerun.peaks.bed"

# chrom_sizes="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/scripts/chrNameLength.txt"

# WORKDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/results/IDR/1to350Ip1_1to350Ip2"
# ============================================================================================================================================================
# STARTING THE SCRIPT
mkdir -p $WORKDIR
cd $WORKDIR
# ============================================================================================================================================================
module load samtools/1.16.1

samtools view -c -F 4 $rep1_clip_bam > rep1_clip.readnum
samtools view -c -F 4 $rep2_clip_bam > rep2_clip.readnum
samtools view -c -F 4 $SMI_input_bam > SMI_input.readnum

/home/tmhaxs421/brannanlab/tmhaxs421/applications/eCLIP-master/bin/overlap_peakfi_with_bam.pl \
    $rep1_clip_bam \
    $SMI_input_bam \
    $rep1_peaks_bed \
    rep1_clip.readnum \
    SMI_input.readnum \
    rep1_normed_peaks.bed 2>/dev/null

echo "Normalize over input fopr rep1 done"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/eCLIP-master/bin/overlap_peakfi_with_bam.pl \
    $rep2_clip_bam \
    $SMI_input_bam \
    $rep2_peaks_bed \
    rep2_clip.readnum \
    SMI_input.readnum \
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
    SMI_input.readnum \
    rep1_normed_peaks.compressed.bed.entropy.full \
    rep1_normed_peaks.compressed.bed.entropy.excessreads 

echo " Computing Entropy for rep1 peaks done"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/make_informationcontent_from_peaks.pl \
    rep2_normed_peaks.compressed.bed.full \
    rep2_clip.readnum \
    SMI_input.readnum \
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
module load cwltool/3.1.2 
mamba activate bioinformatics

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/overlap_peakfi_with_bam.pl \
    $rep1_clip_bam \
    $SMI_input_bam \
    01v02.idr.out.bed \
    rep1_clip.readnum \
    SMI_input.readnum \
    01v02.IDR.out.0102merged.01.bed 2>/dev/null

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/overlap_peakfi_with_bam.pl \
    $rep2_clip_bam \
    $SMI_input_bam \
    01v02.idr.out.bed \
    rep2_clip.readnum \
    SMI_input.readnum \
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
echo "IDR analysis run done, check for errors in the log file"