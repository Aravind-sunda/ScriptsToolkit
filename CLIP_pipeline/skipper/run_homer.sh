#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# run_homer.sh — Standalone HOMER motif analysis from SKIPPER output
#
# Replicates the run_homer Snakemake rule so you can re-run HOMER with
# different flags without touching the rest of the pipeline.
#
# Run from your SKIPPER working directory (where output/ lives).
#
# Usage:
#   ./run_homer.sh --experiment-label LABEL --genome /path/to/genome.fa [OPTIONS]
#
# Required:
#   --experiment-label    Experiment label matching SKIPPER output filenames
#   --genome              Path to genome FASTA (or .2bit) file
#
# Optional path overrides (default to standard SKIPPER output locations):
#   --finemapped          Path to finemapped_windows.bed.gz
#   --background          Path to sampled_fixed_windows.bed.gz
#   --outdir              Output directory
#   --preparsed-dir       HOMER preparsed dir (shared cache across runs)
# =============================================================================

# =============================================================================
# HOMER FLAGS — edit these to customize the motif search
# =============================================================================
HOMER_SIZE="given"       # -size: use actual region sizes (or a fixed bp, e.g. 200)
HOMER_RNA=false           # -rna: search RNA motif database (false to disable)
HOMER_NOFACTS=true       # -nofacts: skip JASPAR/TRANSFAC background (false to disable)
HOMER_S=20               # -S: number of motifs to find
HOMER_LEN="5,6,7,8,9,10,12"   # -len: motif lengths to search (comma-separated)
HOMER_NLEN=1             # -nlen: lower-order oligo normalization length
HOMER_EXTRA_FLAGS="-mset vertebrates"     # any additional flags, e.g. "-p 8 -mset vertebrates"
# =============================================================================

usage() {
    sed -n '/^# Usage/,/^# =====/{ /^# =====/d; s/^# \?//; p }' "$0"
    exit 1
}

EXPERIMENT_LABEL=""
GENOME=""
FINEMAPPED=""
BACKGROUND=""
OUTDIR=""
PREPARSED_DIR="output/homer/preparsed"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --experiment-label) EXPERIMENT_LABEL="$2"; shift 2 ;;
        --genome)           GENOME="$2";            shift 2 ;;
        --finemapped)       FINEMAPPED="$2";        shift 2 ;;
        --background)       BACKGROUND="$2";        shift 2 ;;
        --outdir)           OUTDIR="$2";            shift 2 ;;
        --preparsed-dir)    PREPARSED_DIR="$2";     shift 2 ;;
        -h|--help)          usage ;;
        *) echo "ERROR: Unknown argument: $1"; usage ;;
    esac
done

[[ -z "$EXPERIMENT_LABEL" ]] && { echo "ERROR: --experiment-label is required"; exit 1; }
[[ -z "$GENOME" ]]           && { echo "ERROR: --genome is required"; exit 1; }

# Default paths match the SKIPPER output directory structure
FINEMAPPED="${FINEMAPPED:-output/secondary_results/finemapping/mapped_sites/${EXPERIMENT_LABEL}.finemapped_windows.bed.gz}"
BACKGROUND="${BACKGROUND:-output/homer/region_matched_background/fixed/${EXPERIMENT_LABEL}.sampled_fixed_windows.bed.gz}"
OUTDIR="${OUTDIR:-output/homer/finemapped_results/${EXPERIMENT_LABEL}}"

[[ ! -f "$FINEMAPPED" ]] && { echo "ERROR: finemapped windows not found: $FINEMAPPED"; exit 1; }
[[ ! -f "$BACKGROUND" ]] && { echo "ERROR: background windows not found: $BACKGROUND"; exit 1; }
[[ ! -f "$GENOME" ]]     && { echo "ERROR: genome file not found: $GENOME"; exit 1; }

mkdir -p "$OUTDIR" "$PREPARSED_DIR"

# Build optional boolean flags
BOOL_FLAGS=""
$HOMER_RNA     && BOOL_FLAGS="$BOOL_FLAGS -rna"
$HOMER_NOFACTS && BOOL_FLAGS="$BOOL_FLAGS -nofacts"

echo "[$(date)] Starting HOMER for: $EXPERIMENT_LABEL"
echo "  Foreground  : $FINEMAPPED"
echo "  Background  : $BACKGROUND"
echo "  Genome      : $GENOME"
echo "  Output      : $OUTDIR"
echo "  Preparsed   : $PREPARSED_DIR"
echo "  HOMER flags : -size $HOMER_SIZE -S $HOMER_S -len $HOMER_LEN -nlen $HOMER_NLEN$BOOL_FLAGS $HOMER_EXTRA_FLAGS"

# Use temp files instead of process substitution so paths are stable
TMPDIR_HOMER=$(mktemp -d)
trap 'rm -rf "$TMPDIR_HOMER"' EXIT

FG_TMP="$TMPDIR_HOMER/fg.bed"
BG_TMP="$TMPDIR_HOMER/bg.bed"

echo "[$(date)] Preparing foreground input..."
zcat "$FINEMAPPED" \
    | awk -v OFS="\t" '{print $4 ":" $9, $1, $2+1, $3, $6}' \
    > "$FG_TMP"

echo "[$(date)] Preparing background input..."
zcat "$BACKGROUND" \
    | awk -v OFS="\t" '{print $4, $1, $2+1, $3, $6}' \
    > "$BG_TMP"

echo "[$(date)] Running findMotifsGenome.pl..."
# shellcheck disable=SC2086
findMotifsGenome.pl "$FG_TMP" "$GENOME" "$OUTDIR" \
    -preparsedDir "$PREPARSED_DIR" \
    -size "$HOMER_SIZE" \
    -S "$HOMER_S" \
    -len "$HOMER_LEN" \
    -nlen "$HOMER_NLEN" \
    $BOOL_FLAGS \
    $HOMER_EXTRA_FLAGS \
    -bg "$BG_TMP"

echo "[$(date)] Done. Results in: $OUTDIR"
