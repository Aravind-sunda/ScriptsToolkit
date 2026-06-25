#!/bin/bash
#SBATCH --job-name=sailor_bulk_00
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

# ── PROJECT-LEVEL VARIABLES — update these for every project ──────────────────
# Point HOMEDIR at the same project directory used for MARINE.
# SAILOR reads BAMs from $HOMEDIR/03_star/ and writes to $HOMEDIR/04_sailor/.

export HOMEDIR="/path/to/project"

# ── REQUIRED: these change every project — no defaults are applied ────────────
export LIBRARY="single"            # single | paired
export EDIT_TYPE="CT"              # CT (C>T / A-to-I) | AG (A>G) | etc.
export STRANDEDNESS=2              # from 04_strandedness: 0=unstranded, 1=forward, 2=reverse
export FASTA="/path/to/genome.fa"
export KNOWN_SNPS="/path/to/dbsnp_combined.bed3"

# ── GENOME BUILD (metaPlotR genePred lookup) ──────────────────────────────────
export GENOME="hg38_V44"   # hg19 | hg38 | hg38_V44 | mm10 | mm39
export REFSEQ_DIR="/home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/refseq"

# ── FIXED PATHS (rarely change between projects) ──────────────────────────────
export SNAKEFILE="/home/tmhaxs421/brannanlab/tmhaxs421/applications/FLARE-master/workflow_sailor/Snakefile"
export SINGULARITY_CACHE="/home/tmhaxs421/brannanlab/tmhaxs421/applications/FLARE-master/.singularity_cache"
export METAPLOT_HELPER="/home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/helper_calc_metaplot_dist.py"

# ── DERIVED PATHS ─────────────────────────────────────────────────────────────
export BAM_DIR="$HOMEDIR/03_star"
export OUTPUT_DIR="$HOMEDIR/04_sailor"
export METAPLOT_DIR="$HOMEDIR/05_sailor_metaplotr"
export CONFIG_JSON="$HOMEDIR/scripts/sailor_config.json"

# ── DERIVE SAILOR PARAMS FROM STRANDEDNESS ────────────────────────────────────
# Strandedness 0 (unstranded):
#   reverse_stranded=false — no way to preserve both strands; half the reads are dropped
#   mm_tolerance=2         — relaxed to compensate for unstranded noise
# Strandedness 1 (forward stranded):
#   reverse_stranded=false, mm_tolerance=1
# Strandedness 2 (reverse stranded):
#   reverse_stranded=true,  mm_tolerance=1

case "$STRANDEDNESS" in
    0)
        export REVERSE_STRANDED="false"
        export MM_TOLERANCE=2
        ;;
    1)
        export REVERSE_STRANDED="false"
        export MM_TOLERANCE=1
        ;;
    2)
        export REVERSE_STRANDED="true"
        export MM_TOLERANCE=1
        ;;
    *)
        echo "[ERROR] STRANDEDNESS must be 0, 1, or 2 (got: '$STRANDEDNESS')"
        exit 1
        ;;
esac

# ── VALIDATION ────────────────────────────────────────────────────────────────
_check_var()  { [[ -n "${!1:-}" ]]  || { echo "[ERROR] $1 is not set"; exit 1; }; }
_check_dir()  { [[ -d "$1" ]]       || { echo "[ERROR] Directory not found: $1"; exit 1; }; }
_check_file() { [[ -f "$1" ]]       || { echo "[ERROR] File not found: $1"; exit 1; }; }

_check_var HOMEDIR;             _check_dir  "$HOMEDIR"
_check_var LIBRARY
_check_var EDIT_TYPE
_check_var FASTA;               _check_file "$FASTA"
_check_var KNOWN_SNPS;          _check_file "$KNOWN_SNPS"
_check_dir  "$BAM_DIR"
_check_dir  "$REFSEQ_DIR"
_check_file "$SNAKEFILE"
_check_file "$METAPLOT_HELPER"

# ── LOG CONFIG ────────────────────────────────────────────────────────────────
echo "---------------------"
echo "[LOG] Starting SAILOR Bulk Pipeline at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] HOMEDIR          : $HOMEDIR"
echo "[LOG] LIBRARY          : $LIBRARY"
echo "[LOG] EDIT_TYPE        : $EDIT_TYPE"
echo "[LOG] FASTA            : $FASTA"
echo "[LOG] KNOWN_SNPS       : $KNOWN_SNPS"
echo "[LOG] STRANDEDNESS     : $STRANDEDNESS"
echo "[LOG] REVERSE_STRANDED : $REVERSE_STRANDED"
echo "[LOG] MM_TOLERANCE     : $MM_TOLERANCE"
echo "[LOG] GENOME           : $GENOME"
echo "[LOG] REFSEQ_DIR       : $REFSEQ_DIR"
echo "[LOG] BAM_DIR          : $BAM_DIR"
echo "[LOG] OUTPUT_DIR       : $OUTPUT_DIR"
echo "[LOG] METAPLOT_DIR     : $METAPLOT_DIR"
echo "[LOG] CONFIG_JSON      : $CONFIG_JSON"
echo "---------------------"

SCRIPT_DIR="$(dirname "$(realpath "$0")")"

# ── STEP 1: Run SAILOR snakemake ──────────────────────────────────────────────
echo "---------------------"
echo "[LOG] Step 1: Run SAILOR snakemake pipeline"
echo "---------------------"
bash "$SCRIPT_DIR/01_run_sailor.sh"

# ── STEP 2: metaPlotR metagene distances ──────────────────────────────────────
echo "---------------------"
echo "[LOG] Step 2: metaPlotR metagene distance calculation"
echo "---------------------"
bash "$SCRIPT_DIR/02_sailor_metaplotr.sh"

echo "---------------------"
echo "[LOG] SAILOR Bulk Pipeline completed at $(date '+%Y-%m-%d %H:%M:%S')"
echo "---------------------"
