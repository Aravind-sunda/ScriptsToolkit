# step 1 Make consensus peak bed file from IDR output
bedtools merge -i <(
  cat \
    /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/IDR_filtered_bed3/1to350_filtered_IDR_peaks.bed \
    /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/IDR_filtered_bed3/293_HA_filtered_IDR_peaks.bed \
    /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/IDR_filtered_bed3/351to927_filtered_IDR_peaks.bed \
    /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/IDR_filtered_bed3/600to1132_filtered_IDR_peaks.bed \
    /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/IDR_filtered_bed3/928to1132_filtered_IDR_peaks.bed \
    /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/IDR_filtered_bed3/Delta_new_filtered_IDR_peaks.bed \
    /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/IDR_filtered_bed3/Flip_filtered_IDR_peaks.bed \
  | bedtools sort -i -
) > /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/IDR_filtered_bed3/consensus_peaks.bed3

# Step 2 use the normalized bigwigs to plot correlation heatmap. YOu can use un normalized bigwigs too if you are only usiung one experiment or single sequencing run since the depth should be the same
PEAKS="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/IDR_filtered_bed3/consensus_peaks.bed3"   # your merged BED3
OUT="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/results/peak_correlation/normalized_bigwigs/" 
mkdir -p ${OUT}
multiBigwigSummary BED-file \
  --BED "$PEAKS" \
  -b \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/1to350Ip1_S1_R1_002.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/1to350Ip2_S2_R1_002.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/293-HA-Tert-Ip-1.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/293-HA-Tert-Ip-2.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/351to927Ip1_S5_R1_001.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/351to927Ip2_S6_R1_001.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/600to1132Ip1_S3_R1_002.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/600to1132Ip2_S4_R1_002.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/928to1132Ip1_S7_R1_001.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/928to1132Ip2_S8_R1_001.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/293HATertNRBDdeltaIp1_S3_R1_001..genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/294HATertNRBDdeltaIp1_S4_R1_001..genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/reference/RIAN_seq_rloop/RIAN_BW/RIAN_REP1.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/reference/RIAN_seq_rloop/RIAN_BW/RIAN_REP2.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/FlipIp1.genome-mappedSoSo.rmDupSo_normalized.bw \
  /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/bigwigs_normalized/FlipIp2.genome-mappedSoSo.rmDupSo_normalized.bw \
  -o ${OUT}/correlation_RIAN_normalized.npz \
  --outRawCounts ${OUT}/correlation.rawCounts.normalized.tsv

  
plotCorrelation \
  -in ${OUT}/correlation_RIAN_normalized.npz \
  --corMethod spearman \
  --skipZeros \
  --whatToPlot heatmap \
  --plotNumbers \
  --labels 1to350_IP1 1to350_IP2 293HA_IP1 293HA_IP2 351to927_IP1 351to927_IP2 600to1132_IP1 600to1132_IP2 928to1132_IP1 928to1132_IP2 delta_IN1 delta_IP1 RIAN_REP1 RIAN_REP2 \
  --outFileCorMatrix ${OUT}/correlation_RIAN_normalized.spearman.matrix.tsv \
  --plotFile ${OUT}/correlation_RIAN_normalized.spearman.heatmap.pdf  
