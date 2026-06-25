#!/bin/bash
#SBATCH --job-name=marine_bulk_pipeline_00
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

export HOMEDIR="/path/to/project"
export SAMPLESHEET="$HOMEDIR/samplesheet.csv"   # columns: sample,r1,r2,libraryType (SE or PE)

# reference files
export GENOMEDIR="/path/to/STAR/index"      # STAR genome index (built with STAR v2.5.2b)
export FASTA="/path/to/genome.fa"           # genome FASTA (for MARINE/SAILOR)
export GTF="/path/to/annotation.gtf"        # GTF annotation (featureCounts + MARINE)
export GENE_BED="/path/to/gene_models.bed"  # BED6 for infer_experiment.py (RSeQC) and edit-site annotation
export DBSNP_BED="/path/to/dbsnp_combined.bed3"  # dbSNP BED for edit filtering (step 06)
export GENOME="hg38_V44"    # genome build: hg19, hg38, hg38_V44, mm10, mm39 — used by step 07 to select genePred file
export REFSEQ_DIR="/home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/refseq"  # directory containing genePred .txt.gz files (step 07)
export EDIT_TYPE="C>T"  # strand_conversion type to filter/normalize (step 06) and metagene distances (step 07)

# ── VALIDATION ────────────────────────────────────────────────────────────────

if [ -z "$HOMEDIR" ]; then
    echo "[ERROR] HOMEDIR is not set"
    exit 1
fi

if [ ! -d "$HOMEDIR" ]; then
    echo "[ERROR] HOMEDIR does not exist: $HOMEDIR"
    exit 1
fi

if [ ! -f "$SAMPLESHEET" ]; then
    echo "[ERROR] Samplesheet not found: $SAMPLESHEET"
    exit 1
fi

if [ ! -d "$GENOMEDIR" ]; then
    echo "[ERROR] STAR genome index not found: $GENOMEDIR"
    exit 1
fi

if [ ! -f "$FASTA" ]; then
    echo "[ERROR] Genome FASTA not found: $FASTA"
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "[ERROR] GTF not found: $GTF"
    exit 1
fi

if [ ! -f "$GENE_BED" ]; then
    echo "[ERROR] Gene BED file not found: $GENE_BED"
    exit 1
fi

if [ ! -f "$DBSNP_BED" ]; then
    echo "[ERROR] dbSNP BED file not found: $DBSNP_BED"
    exit 1
fi

if [ ! -d "$REFSEQ_DIR" ]; then
    echo "[ERROR] REFSEQ_DIR not found: $REFSEQ_DIR"
    exit 1
fi


# ── PIPELINE ──────────────────────────────────────────────────────────────────

echo "---------------------"
echo "[LOG] Starting MARINE Bulk RNA-seq Pipeline at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] HOMEDIR    : $HOMEDIR"
echo "[LOG] SAMPLESHEET: $SAMPLESHEET"
echo "[LOG] GENOMEDIR  : $GENOMEDIR"
echo "[LOG] FASTA      : $FASTA"
echo "[LOG] GTF        : $GTF"
echo "[LOG] GENE_BED   : $GENE_BED"
echo "[LOG] DBSNP_BED  : $DBSNP_BED"
echo "[LOG] GENOME     : $GENOME"
echo "[LOG] REFSEQ_DIR : $REFSEQ_DIR"
echo "[LOG] EDIT_TYPE  : $EDIT_TYPE"
echo "---------------------"

echo "---------------------"
echo "[LOG] Step 1: QC (FastQC, seqkit stats, MultiQC)"
echo "---------------------"
bash 01_qc.sh

echo "---------------------"
echo "[LOG] Step 2: Adapter Trimming (fastp)"
echo "---------------------"
bash 02_cutadapt.sh

echo "---------------------"
echo "[LOG] Step 3: Alignment (STAR)"
echo "---------------------"
bash 03_star_align.sh

echo "---------------------"
echo "[LOG] Step 4: Infer Strandedness and MARINE"
echo "---------------------"
bash 04_infer_strandedness_marine.sh

echo "---------------------"
echo "[LOG] Step 5: featureCounts"
echo "---------------------"
bash 05_featurecounts.sh

echo "---------------------"
echo "[LOG] Step 6: Filter and Normalize MARINE edits"
echo "---------------------"
bash 06_marine_filter_normalize.sh

echo "---------------------"
echo "[LOG] Step 7: metaPlotR metagene plots"
echo "---------------------"
bash 07_metaplotr.sh

echo "---------------------"
echo "[LOG] MARINE Bulk RNA-seq Pipeline Completed at $(date '+%Y-%m-%d %H:%M:%S')"
echo "---------------------"
