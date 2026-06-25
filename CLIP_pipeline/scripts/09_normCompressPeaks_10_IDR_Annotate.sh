#!/bin/bash

#SBATCH --partition=defq
#SBATCH --job-name=normCompressPeaks_IDR_Annotate
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

# chrom_sizes="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/TERTNRBDscripts/reference/Encode/hg38_chrom_sizes.txt"

HOMEDIR=""
WORKDIR="$HOMEDIR/09_normCompressPeaks"

# ============================================================================================================================================================
# You probably do not have to change the following variables, They are defined so that you know what you are using. The script locations will not change
# TODO: CHECK IF ALL THE SCRIPTS ARE COORECT THAT IS BEING USED
overlap_peakfi_with_bam_pl="/home/tmhaxs421/brannanlab/tmhaxs421/applications/eCLIP-master/bin/overlap_peakfi_with_bam.pl"
compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_pl="/home/tmhaxs421/brannanlab/tmhaxs421/applications/eCLIP-master/bin/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl"    
compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull_pl="/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl"
get_reproducing_peaks_pl="/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/get_reproducing_peaks.pl"
# ============================================================================================================================================================
# STARTING THE SCRIPT
mkdir -p $WORKDIR
cd $WORKDIR
# ============================================================================================================================================================
rep1_clip_basename=$(basename $rep1_clip_bam .genome-mappedSoSo.rmDupSo.bam)
rep2_clip_basename=$(basename $rep2_clip_bam .genome-mappedSoSo.rmDupSo.bam)
rep1_input_basename=$(basename $rep1_input_bam .genome-mappedSoSo.rmDupSo.bam)
rep2_input_basename=$(basename $rep2_input_bam .genome-mappedSoSo.rmDupSo.bam)

# get base name of rep1_clip_bam
# /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/CLIP_aravind/delta_NRBD/clipper_output/293HATertNRBDdeltaIn1_S1_R1_001.fq.repbase_mapped.Unmapped.out.fq.repbase.genome-mapped_hg38.Aligned.out.rerun.peaks.bed
samtools view -cF 4 $rep1_clip_bam > "$WORKDIR/$rep1_clip_basename.readnum"
samtools view -cF 4 $rep2_clip_bam > "$WORKDIR/$rep2_clip_basename.readnum"

samtools view -cF 4 $rep1_input_bam > "$WORKDIR/$rep1_input_basename.readnum"
samtools view -cF 4 $rep2_input_bam > "$WORKDIR/$rep2_input_basename.readnum"

# store readnumber file path in variables
rep1_clip_readnum="$WORKDIR/$rep1_clip_basename.readnum"
rep2_clip_readnum="$WORKDIR/$rep2_clip_basename.readnum"
rep1_input_readnum="$WORKDIR/$rep1_input_basename.readnum"
rep2_input_readnum="$WORKDIR/$rep2_input_basename.readnum"


perl $overlap_peakfi_with_bam_pl \
    $rep1_clip_bam \
    $rep1_input_bam \
    $rep1_peaks_bed \
    $rep1_clip_readnum \
    $rep1_input_readnum \
    $WORKDIR/$rep1_clip_basename.peakClusters.normed.bed 2>/dev/null


perl $compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_pl \
    $WORKDIR/$rep1_clip_basename.peakClusters.normed.bed \
    $WORKDIR/$rep1_clip_basename.peakClusters.normed.compressed.bed

# do the above for replicate 2 as well
perl $overlap_peakfi_with_bam_pl \
    $rep2_clip_bam \
    $rep2_input_bam \
    $rep2_peaks_bed \
    $rep2_clip_readnum \
    $rep2_input_readnum \
    $WORKDIR/$rep2_clip_basename.peakClusters.normed.bed 2>/dev/null

perl $compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_pl \
    $WORKDIR/$rep2_clip_basename.peakClusters.normed.bed \
    $WORKDIR/$rep2_clip_basename.peakClusters.normed.compressed.bed
# # ============================================================================================================================================================
# starting IDR analysis

# step 100 is the same as the first half of step 93. Step 101 is different from the second half of step 93. Step 93 second half does not output ful bed files
# the script for step 93 and step 101 is fifferent.

# takes the .full output from overlappeakfi_with_bam.pl and the output from compress bed and make the compressed.bed.full file
perl $compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull_pl \
    $WORKDIR/$rep1_clip_basename.peakClusters.normed.bed.full \
    $WORKDIR/$rep1_clip_basename.peakClusters.normed.compressed.bed \
    $WORKDIR/$rep1_clip_basename.peakClusters.normed.compressed.bed.full

echo "Compress rep1 peaks done"

perl $compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull_pl \
    $WORKDIR/$rep2_clip_basename.peakClusters.normed.bed.full \
    $WORKDIR/$rep2_clip_basename.peakClusters.normed.compressed.bed \
    $WORKDIR/$rep2_clip_basename.peakClusters.normed.compressed.bed.full

echo "Compress rep2 peaks done"

# ============================================================================================================================================================
# Annotating per-replicate normalized compressed peaks (Step 09a)
module purge
module load mamba
mamba activate
mamba activate annotator

SIF="/condo/brannanlab/tmhaxs421/applications/annotator/annotator.sif"
GTFDB="/condo/brannanlab/VS_share/Genome_indexs/GTFdb/Yeolab/gffutils_dbs/gencode.v40.annotation.gtf.db"
ANNOTDIR="$HOMEDIR/09a_annotated_normalized_clipper_peaks"
mkdir -p "$ANNOTDIR"

for basename in "$rep1_clip_basename" "$rep2_clip_basename"; do
    input_bed="$WORKDIR/${basename}.peakClusters.normed.compressed.bed"
    sorted_bed="$ANNOTDIR/${basename}.peakClusters.normed.compressed.sorted.bed"
    annotated_bed="$ANNOTDIR/${basename}.peakClusters.normed.compressed.sorted.annotated.bed"

    echo "Sorting $basename compressed peaks..."
    sort -k1,1 -k2,2n "$input_bed" -o "$sorted_bed"

    echo "Annotating $basename compressed peaks..."
    apptainer exec \
        --cleanenv \
        --env PYTHONNOUSERSITE=1 \
        --env GTFDB="$GTFDB" \
        --bind "$HOMEDIR",/condo/brannanlab \
        "$SIF" \
        annotator \
            --output "$annotated_bed" \
            --input  "$sorted_bed" \
            --gtfdb "$GTFDB"
done
echo "Per-replicate annotation done"

# ============================================================================================================================================================
# Restore modules for entropy/IDR steps
module purge
module load mamba
mamba init bash
module load shared
module load mergepeaks/0.1.0
module load cwltool/3.1.20250925164626
mamba activate bioinformatics

echo "Computing Entropy for rep1 and rep2 peaks"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/make_informationcontent_from_peaks.pl \
    $WORKDIR/$rep1_clip_basename.peakClusters.normed.compressed.bed.full \
    $WORKDIR/$rep1_clip_basename.readnum \
    $WORKDIR/$rep1_input_basename.readnum \
    $WORKDIR/$rep1_clip_basename.compressed.bed.entropy.full \
    $WORKDIR/$rep1_clip_basename.compressed.bed.entropy.excessreads 

echo " Computing Entropy for rep1 peaks done"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/make_informationcontent_from_peaks.pl \
    $WORKDIR/$rep2_clip_basename.peakClusters.normed.compressed.bed.full \
    $WORKDIR/$rep2_clip_basename.readnum \
    $WORKDIR/$rep2_input_basename.readnum \
    $WORKDIR/$rep2_clip_basename.compressed.bed.entropy.full \
    $WORKDIR/$rep2_clip_basename.compressed.bed.entropy.excessreads 

echo "Computing Entropy for rep2 peaks done"

echo "Converting compressed full to bed format"
/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/full_to_bed.py \
    --input $WORKDIR/$rep1_clip_basename.compressed.bed.entropy.full \
    --output $WORKDIR/$rep1_clip_basename.compressed.bed.entropy.bed
echo "Convert rep1 peaks to bed done"

/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/full_to_bed.py \
    --input $WORKDIR/$rep2_clip_basename.compressed.bed.entropy.full \
    --output $WORKDIR/$rep2_clip_basename.compressed.bed.entropy.bed

echo "Convert rep2 peaks to bed done"
# ============================================================================================================================================================
# running Part 10 of ths script to run the IDR analysis
# this idr environment must be installed carefully
# create the environment with the following command and update the script with the 2.0.3 application
# mamba create -n idr -c conda-forge -c bioconda idr=2.0.2
# mamba activate idr
# python3 /home/tmhaxs421/brannanlab/tmhaxs421/applications/idr-2.0.3/setup.py install
# removes shared modules since it has a more recent version of numpy it will conflict with the idr version of numpy
echo "Running IDR on peaks from rep1 and rep2"

module purge
module load mamba
mamba init bash
mamba activate idr

# adding 0's to 3
awk 'BEGIN{OFS="\t"}{print $0, 0, 0, 0}' $WORKDIR/$rep1_clip_basename.compressed.bed.entropy.bed > $WORKDIR/$rep1_clip_basename.compressed.bed.entropy.bed.9col
awk 'BEGIN{OFS="\t"}{print $0, 0, 0, 0}' $WORKDIR/$rep2_clip_basename.compressed.bed.entropy.bed > $WORKDIR/$rep2_clip_basename.compressed.bed.entropy.bed.9col

echo "Running IDR analysis"
WORKDIR2="$HOMEDIR/10_idr"
mkdir -p $WORKDIR2

idr \
    --samples "$WORKDIR/$rep1_clip_basename.compressed.bed.entropy.bed.9col" "$WORKDIR/$rep2_clip_basename.compressed.bed.entropy.bed.9col" \
    --input-file-type bed \
    --rank 5 \
    --peak-merge-method max \
    --plot \
    -o $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}.idr.out

# $HOMEDIR/idr/${rep1_clip_basename}_${rep2_clip_basename}.idr.out
/home/tmhaxs421/brannanlab/tmhaxs421/applications/merge_peaks-0.1.0/bin/perl/parse_idr_peaks.pl \
    $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}.idr.out \
    $WORKDIR/${rep1_clip_basename}.compressed.bed.entropy.full \
    $WORKDIR/${rep2_clip_basename}.compressed.bed.entropy.full \
    $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}.idr.out.bed
    
echo "IDR analysis done"
# ============================================================================================================================================================
module purge
module load mamba
mamba init bash
module load shared
module load mergepeaks/0.1.0
module load cwltool/3.1.20250925164626 
mamba activate bioinformatics

perl $overlap_peakfi_with_bam_pl \
    $rep1_clip_bam \
    $rep1_input_bam \
    $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}.idr.out.bed \
    $rep1_clip_readnum \
    $rep1_input_readnum \
    $WORKDIR2/${rep1_clip_basename}.idr.out.0102merged.01.bed 2>/dev/null

perl $overlap_peakfi_with_bam_pl \
    $rep2_clip_bam \
    $rep2_input_bam \
    $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}.idr.out.bed \
    $rep2_clip_readnum \
    $rep2_input_readnum \
    $WORKDIR2/${rep2_clip_basename}.idr.out.0102merged.02.bed 2>/dev/null

perl $get_reproducing_peaks_pl \
    $WORKDIR2/${rep1_clip_basename}.idr.out.0102merged.01.bed.full \
    $WORKDIR2/${rep2_clip_basename}.idr.out.0102merged.02.bed.full \
    $WORKDIR2/${rep1_clip_basename}_reproducible_peaks.01.bed.full \
    $WORKDIR2/${rep2_clip_basename}_reproducible_peaks.02.bed.full \
    $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.bed \
    $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.custombed \
    $WORKDIR/${rep1_clip_basename}.compressed.bed.entropy.full \
    $WORKDIR/${rep2_clip_basename}.compressed.bed.entropy.full \
    $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}.idr.out


echo "Get reproducible peaks done"
echo "IDR analysis run done"
#===========================================================================================================================================================
# Sort and Annotating the IDR peaks in bed files

# must sort the bed file before annotating
sort -k1,1 -k2,2n "${WORKDIR2}/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.bed" -o "${WORKDIR2}/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.sorted.bed"

module purge
module load mamba
mamba activate
mamba activate annotator

SIF="/condo/brannanlab/tmhaxs421/applications/annotator/annotator.sif"
GTFDB="/condo/brannanlab/VS_share/Genome_indexs/GTFdb/Yeolab/gffutils_dbs/gencode.v40.annotation.gtf.db"

apptainer exec \
    --cleanenv \
    --env PYTHONNOUSERSITE=1 \
    --env rep1_clip_basename="$rep1_clip_basename" \
    --env rep2_clip_basename="$rep2_clip_basename" \
    --env WORKDIR2="$WORKDIR2" \
    --env GTFDB="$GTFDB" \
    --bind "$WORKDIR2","$WORKDIR",/condo/brannanlab \
    "$SIF" \
    annotator \
        --output "$WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.sorted.annotated.bed" \
        --input  "$WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.sorted.bed" \
        --gtfdb "$GTFDB"



# use apptainer shell if you want to run on command line and apptainer exec if you want to submit as a job
# SIF="/condo/brannanlab/tmhaxs421/applications/annotator/annotator.sif"

# apptainer shell   --cleanenv   \
#     --env PYTHONNOUSERSITE=1 \
#     --env rep1_clip_basename="$rep1_clip_basename" \
#     --env rep2_clip_basename="$rep2_clip_basename"  \
#     --env WORKDIR2="$WORKDIR2" \
#     --bind "$WORKDIR2","$WORKDIR1",/condo/brannanlab   \
#     "$SIF"

# GTFDB="/condo/brannanlab/VS_share/Genome_indexs/GTFdb/Yeolab/gffutils_dbs/gencode.v40.annotation.gtf.db"

# annotator \
#     --output $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.sorted.annotated.bed \
#     --input $WORKDIR2/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.sorted.bed \
#     --gtfdb $GTFDB