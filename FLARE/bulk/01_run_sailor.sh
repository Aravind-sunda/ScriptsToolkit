#!/bin/bash
#SBATCH --job-name=sailor_run_01
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=96:00:00
#SBATCH --output=slurm_%u_%x_%j.log

set -euo pipefail

# All variables exported from 00_saiilor_bulk_pipeline.sh:
#   BAM_DIR, OUTPUT_DIR, CONFIG_JSON, SNAKEFILE, SINGULARITY_CACHE
#   LIBRARY, EDIT_TYPE, REVERSE_STRANDED, FASTA, KNOWN_SNPS

SCRIPT_DIR="$(dirname "$(realpath "$0")")"

module load snakemake/9.13.4

mkdir -p "$OUTPUT_DIR"
mkdir -p "$SINGULARITY_CACHE"
mkdir -p "$(dirname "$CONFIG_JSON")"

echo "---------------------"
echo "[LOG] Starting SAILOR pipeline at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] Snakefile  : $SNAKEFILE"
echo "[LOG] Config     : $CONFIG_JSON"
echo "[LOG] BAM_DIR    : $BAM_DIR"
echo "[LOG] Output     : $OUTPUT_DIR"
echo "[LOG] CPUs       : ${SLURM_CPUS_PER_TASK:-36}"
echo "---------------------"

# ── GENERATE CONFIG ───────────────────────────────────────────────────────────
if [ ! -f "$CONFIG_JSON" ]; then
    echo "[LOG] Config JSON not found — generating it now"
    python3 "$SCRIPT_DIR/helper_make_sailor_json.py" \
        --samples_path     "$BAM_DIR" \
        --output_dir       "$OUTPUT_DIR" \
        --output_json      "$CONFIG_JSON" \
        --library          "$LIBRARY" \
        --edit_type        "$EDIT_TYPE" \
        --reverse_stranded "$REVERSE_STRANDED" \
        --reference_fasta  "$FASTA" \
        --known_snps       "$KNOWN_SNPS" \
        --mm_tolerance     "$MM_TOLERANCE"
else
    echo "[LOG] Config JSON already exists — skipping generation: $CONFIG_JSON"
fi

# ── RUN SAILOR ────────────────────────────────────────────────────────────────
snakemake \
    --snakefile      "$SNAKEFILE" \
    --configfile     "$CONFIG_JSON" \
    --cores          "${SLURM_CPUS_PER_TASK:-36}" \
    --use-singularity \
    --singularity-args "--bind /home/tmhaxs421/brannanlab" \
    --singularity-prefix "$SINGULARITY_CACHE" \
    --latency-wait   120 \
    --printshellcmds \
    --verbose \
    --rerun-incomplete

echo "---------------------"
echo "[LOG] SAILOR pipeline completed at $(date '+%Y-%m-%d %H:%M:%S')"
echo "---------------------"
