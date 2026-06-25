#!/bin/bash

#SBATCH --partition=defq
#SBATCH --job-name=CLIP_IDR_delta
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
rep1_clip_bam="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_Tert_NRBD/dedup_sort/DeltaIp1_S5_R1_001.genome-mappedSoSo.rmDupSo.bam"
rep2_clip_bam="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_Tert_NRBD/dedup_sort/DeltaIp2.genome-mappedSoSo.rmDupSo.bam"
rep1_input_bam="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_Tert_NRBD/dedup_sort/DeltaIn1_S1_R1_001.genome-mappedSoSo.rmDupSo.bam"
rep2_input_bam="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_Tert_NRBD/dedup_sort/DeltaIn2_S2_R1_001.genome-mappedSoSo.rmDupSo.bam"

# Only peaks of IP files
rep1_peaks_bed="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_Tert_NRBD/Clipper_output/DeltaIp1_S5_R1_001.rerun.peaks.bed"
rep2_peaks_bed="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_Tert_NRBD/Clipper_output/DeltaIp2.rerun.peaks.bed"
# Format of peaks bed file
# chr1	15306	15338	ENSG00000227232.5_0_23	1.543195658818628e-08	-	15320	15324

chrom_sizes="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/scripts/chrNameLength.txt"


WORKDIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/293_Tert_NRBD/IDR/Delta_293_Tert_NRBD"
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



# # ============================================================================================================================================================
# # PART-2: STEPS 106 TO 112
# # ============================================================================================================================================================
# # We are going to use the same input files that we used for the first part of the analysis

# # rep1_clip_bam
# # rep2_clip_bam
# # SMI_input_bam
# # rep1_peaks_bed
# # rep2_peaks_bed

# echo "Running steps 106 to 112 for reproducibility analysis"

# mkdir -p $WORKDIR/replicate_reproducibility_mergeandsplit
# cd $WORKDIR/replicate_reproducibility_mergeandsplit

# # STEP 106: merging the replicates 
# samtools merge -@ 36 merged.bam \
#     $rep1_clip_bam \
#     $rep2_clip_bam

# # STEP 107: generating 2 random subsets of the merged bam file
# NLINES=$(samtools view -cF 4 merged.bam)
# HALFNLINES=$(($NLINES / 2))
# samtools view -@ 36 merged.bam | shuf | split -d -l ${HALFNLINES} - merged.bam 
# samtools view -@ 36 -H merged.bam | cat - merged.bam00 | samtools view -@ 36 -bS - > merged00.bam
# samtools view -@ 36 -H merged.bam | cat - merged.bam01 | samtools view -@ 36 -bS - > merged01.bam
# samtools sort -@ 36 merged00.bam \
# -o IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.bam 
# samtools sort -@ 36 merged01.bam \
# -o IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.bam

# # We will be merging the rep1_input_bam rep2_input_bam into a single SMInput file since we are generating the pesudo replicated from both the samples
# samtools merge -@ 36 merged_input.bam \
#     $rep1_input_bam \
#     $rep2_input_bam

# SMinput="$WORKDIR/replicate_reproducibility_mergeandsplit/merged_input.bam"

# # Step 108: Generating Pseudo peaks using CLIPPER for the merged bam file
# # running clipper from apptainer

# module purge
# module load apptainer


# SIF="/cm/shared/apptainer-images/clipper_5d865bb.sif"
# condo="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis"
# home="/home/tmhaxs421"
# INPUTDIR="$WORKDIR/replicate_reproducibility_mergeandsplit"

# cd $INPUTDIR

# apptainer run \
#     --bind /condo/,$home,$condo,$INPUTDIR \
#     $SIF \
#     clipper \
#         --species GRCh38_v29e \
#         --bam IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.bam \
#         --processors=36 \
#         --save-pickle \
#         --outfile IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.peakClusters.bed

# apptainer run \
#     --bind /condo/,$home,$condo,$INPUTDIR \
#     $SIF \
#     clipper \
#         --species GRCh38_v29e \
#         --bam $INPUTDIR/IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.bam \
#         --processors=36 \
#         --save-pickle \
#         --outfile $INPUTDIR/IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.peakClusters.bed


# echo "Clipper analysis done for merged and split bam files"

# # Step 109: Generating the pseudo peaks for the merged bam file
# module load mamba
# mamba init bash
# module load shared
# module load mergepeaks/0.1.0
# module load cwltool/3.1.20250925164626 
# mamba activate bioinformatics

# echo "Running steps 100 to 105 for reproducibility analysis for merge and split bam files \n"

# /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/scripts/steps100_105_for_reproducibility_single_SMI.sh \
#     IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.bam \
#     IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.bam \
#     $SMinput \
#     IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.peakClusters.bed \
#     IP.merged.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.peakClusters.bed \
#     $chrom_sizes \
#     $WORKDIR/replicate_reproducibility_mergeandsplit

# # Step 110: Split each bam file without merging the input bam files
# mkdir -p $WORKDIR/replicate_reproducibility_split_rep1
# cd $WORKDIR/replicate_reproducibility_split_rep1

# NLINES=$(samtools view -cF 4 $rep1_clip_bam)
# HALFNLINES=$(($NLINES / 2))

# samtools view $rep1_clip_bam | shuf | split -d -l ${HALFNLINES} - $WORKDIR/replicate_reproducibility_split_rep1/rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam

# samtools view -H $rep1_clip_bam | cat - rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam00 | samtools view -bS - > rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.bam 

# samtools view -H $rep1_clip_bam | cat - rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam01 | samtools view -bS - > rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.bam 

# samtools sort rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.bam \
#     -o rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.bam

# samtools sort rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.bam \
#     -o rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.bam

# SIF="/cm/shared/apptainer-images/clipper_5d865bb.sif"
# INPUTDIR="$WORKDIR/replicate_reproducibility_split_rep1"
# # condo="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis"
# home="/home/tmhaxs421"

# # Step 111: Generating the pseudo peaks for the split bam file
# apptainer run \
#     --bind /condo/,$home,$WORKDIR,$INPUTDIR \
#     $SIF \
#     clipper \
#     --species GRCh38_v29e \
#     --bam rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.bam \
#     --save-pickle \
#     --outfile rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.peakClusters.bed

# apptainer run \
#     --bind /condo/,$home,$WORKDIR,$INPUTDIR \
#     $SIF \
#     clipper \
#     --species GRCh38_v29e \
#     --bam rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.bam \
#     --save-pickle \
#     --outfile rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.peakClusters.bed

# echo "Clipper analysis done for split bam files"
# # Step 112: Generating the pseudo peaks for the split bam file


# /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/scripts/steps100_105_for_reproducibility_single_SMI.sh \
#     rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.bam \
#     rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.bam \
#     $rep1_input_bam \
#     rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.peakClusters.bed \
#     rep1.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.peakClusters.bed \
#     $chrom_sizes \
#     $WORKDIR/replicate_reproducibility_split_rep1

# # ============================================================================================================================================================
# # Step 113: Rep 2 analysis of steps 110 to 112
# mkdir -p $WORKDIR/replicate_reproducibility_split_rep2
# cd $WORKDIR/replicate_reproducibility_split_rep2

# NLINES=$(samtools view -cF 4 $rep2_clip_bam)
# HALFNLINES=$(($NLINES / 2))

# samtools view $rep2_clip_bam | shuf | split -d -l ${HALFNLINES} - $WORKDIR/replicate_reproducibility_split_rep2/rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam

# samtools view -H $rep2_clip_bam | cat - rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam00 | samtools view -bS - > rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.bam 

# samtools view -H $rep2_clip_bam | cat - rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam01 | samtools view -bS - > rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.bam 

# samtools sort rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.bam \
#     -o rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.bam

# samtools sort rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.bam \
#     -o rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.bam

# SIF="/cm/shared/apptainer-images/clipper_5d865bb.sif"
# INPUTDIR="$WORKDIR/replicate_reproducibility_split_rep2"
# # condo="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis"
# home="/home/tmhaxs421"

# module load apptainer
# apptainer run \
#     --bind /condo/,$home,$WORKDIR,$INPUTDIR \
#     $SIF \
#     clipper \
#     --species GRCh38_v29e \
#     --bam rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.bam \
#     --save-pickle \
#     --outfile rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.peakClusters.bed

# apptainer run \
#     --bind /condo/,$home,$WORKDIR,$INPUTDIR \
#     $SIF \
#     clipper \
#     --species GRCh38_v29e \
#     --bam rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.bam \
#     --save-pickle \
#     --outfile rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.peakClusters.bed

# echo "Clipper analysis done for split bam files"
# # Step 112: Generating the pseudo peaks for the split bam file

# /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/scripts/steps100_105_for_reproducibility_single_SMI.sh \
#     rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.bam \
#     rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.bam \
#     $rep2_input_bam \
#     rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split0.sorted.peakClusters.bed \
#     rep2.IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.split1.sorted.peakClusters.bed \
#     $chrom_sizes \
#     $WORKDIR/replicate_reproducibility_split_rep2