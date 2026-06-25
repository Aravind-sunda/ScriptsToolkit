#!/usr/bin/env bash
set -euo pipefail

CSV_PATH="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/MassSpec_data/Individual_comparisons"
OUT_BASE="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/MassSpec_data/output"
SCRIPT="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/MassSpec_data/limma_de.R"

mkdir -p "$OUT_BASE"

shopt -s nullglob
for file in "$CSV_PATH"/*.csv; do
  base=$(basename "$file" .csv)
  outdir="$OUT_BASE/$base"
  mkdir -p "$outdir"

  echo "[`date`] Running: $base"
  Rscript "$SCRIPT" "$file" "$outdir" 
done &> "$OUT_BASE/log.txt"


# Rscript /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/MassSpec_data/limma_de.R \
#   /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/MassSpec_data/haec/HAEC_cells_ip.csv \
#   /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/analysis/MassSpec_data/haec/limma_output