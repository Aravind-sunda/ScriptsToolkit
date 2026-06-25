#!/bin/bash

# 00_marine_sc_pipeline.sh
# Master driver for the MARINE single-cell RNA-editing pipeline.
#
# Supports two starting points controlled by START_FROM:
#
#   START_FROM="fastq"  — runs CellRanger first, then chains MARINE per sample.
#                         Samplesheet columns: sample,fastq_dir
#
#   START_FROM="bam"    — skips CellRanger; BAMs are already available.
#                         Samplesheet columns: sample,bam,matrixdir
#
# ── FASTQ NAMING REQUIREMENT (START_FROM="fastq" only) ───────────────────────
# CellRanger discovers reads by filename, not by explicit R1/R2 paths.
# Files in each sample's fastq_dir MUST follow the 10x Genomics naming scheme:
#
#   {sample}_S{n}_L{lane}_R1_001.fastq.gz   ← barcodes + UMI (28 bp)
#   {sample}_S{n}_L{lane}_R2_001.fastq.gz   ← cDNA read
#   {sample}_S{n}_L{lane}_I1_001.fastq.gz   ← index read (optional)
#
# FASTQs from Illumina mkfastq / 10x mkfastq are already correctly named.
# Downloaded FASTQs (SRA/GEO) must be renamed before running this pipeline:
#   mv sample_R1.fastq.gz  SAMPLE1_S1_L001_R1_001.fastq.gz
#   mv sample_R2.fastq.gz  SAMPLE1_S1_L001_R2_001.fastq.gz
# ─────────────────────────────────────────────────────────────────────────────
#
# After MARINE jobs finish, run on the head node:
#   bash 02_filter_marine_edits_sc.sh
#   bash 03_normalize_marine_edits_umi_sc.sh
#
# Usage: bash 00_marine_sc_pipeline.sh

# ── STARTING POINT ─────────────────────────────────────────────────────────────

START_FROM="bam"   # "bam"  → samplesheet: sample,bam,matrixdir
                   # "fastq"→ samplesheet: sample,fastq_dir

# ── PROJECT-LEVEL VARIABLES — update these for every project ──────────────────

export HOMEDIR="/path/to/project"
export SAMPLESHEET="${HOMEDIR}/samplesheet.csv"

# Reference files
export GENOME_FA="/path/to/genome.fa"              # genome FASTA — required by samtools calmd
export ANNOTATION_BED="/path/to/gene_models.bed"  # BED6 for MARINE edit-site annotation
export DBSNP_BED="/path/to/dbsnp_combined.bed3"   # dbSNP BED for filtering (step 02)

# MARINE container + parameters
export CONTAINER="/path/to/marine_container.sif"  # Singularity SIF with MARINE + samtools
export MARINE_PY="/opt/MARINE/marine.py"           # path inside the container
export BARCODE_TAG="CB"                            # CellRanger barcode tag (default: CB)
export STRANDEDNESS=2                              # 0=unstranded, 1=forward, 2=reverse
export MIN_READ_QUALITY=255                        # CellRanger-mapped reads
export MIN_BASE_QUALITY=30

# Filter parameters (used by 02_filter_marine_edits_sc.sh)
export MIN_COUNT=3
export MAX_FRAC=0.1

# Normalize parameters (used by 03_normalize_marine_edits_umi_sc.sh)
export EDIT_TYPE="C>T"    # strand conversion to retain: "C>T" (A-to-I, reverse strand) or "A>G" (sense strand)

# Output directories (shared with steps 02 and 03)
export CELLRANGER_OUTDIR="${HOMEDIR}/00_cellranger"   # only used when START_FROM="fastq"
export MARINE_OUTDIR="${HOMEDIR}/01_marine"
export FILTER_OUTDIR="${HOMEDIR}/02_filter"
export NORMALIZE_OUTDIR="${HOMEDIR}/03_normalize"

# LOG_DIR="${HOMEDIR}/logs/sc_pipeline"
MAIL_USER="your_email@institution.edu"

# ── CellRanger RESOURCES (START_FROM="fastq" only) ─────────────────────────────
# Partition, CPUs, memory, and walltime are set as #SBATCH directives in 00_cellranger_sc.sh.

TRANSCRIPTOME_REF="/path/to/cellranger_reference"   # pre-built CellRanger reference

# ── MARINE RESOURCES ──────────────────────────────────────────────────────────
# Partition, CPUs, memory, and walltime are set as #SBATCH directives in 01_run_marine_sc.sh.

# ── VALIDATION ────────────────────────────────────────────────────────────────

[ ! -f "${SAMPLESHEET}" ]    && echo "[ERROR] Samplesheet not found: ${SAMPLESHEET}"          && exit 1
[ ! -f "${GENOME_FA}" ]      && echo "[ERROR] GENOME_FA not found: ${GENOME_FA}"               && exit 1
[ ! -f "${ANNOTATION_BED}" ] && echo "[ERROR] ANNOTATION_BED not found: ${ANNOTATION_BED}"    && exit 1
[ ! -f "${DBSNP_BED}" ]      && echo "[ERROR] DBSNP_BED not found: ${DBSNP_BED}"              && exit 1
[ ! -f "${CONTAINER}" ]      && echo "[ERROR] Container not found: ${CONTAINER}"               && exit 1

if [[ "${START_FROM}" == "fastq" ]]; then
    [ ! -d "${TRANSCRIPTOME_REF}" ] && echo "[ERROR] TRANSCRIPTOME_REF not found: ${TRANSCRIPTOME_REF}" && exit 1
    mkdir -p "${CELLRANGER_OUTDIR}"
elif [[ "${START_FROM}" != "bam" ]]; then
    echo "[ERROR] START_FROM must be 'fastq' or 'bam', got: '${START_FROM}'"
    exit 1
fi

mkdir -p "${MARINE_OUTDIR}" "${FILTER_OUTDIR}" "${NORMALIZE_OUTDIR}" 

# "${LOG_DIR}"

# ── HEADER LOG ────────────────────────────────────────────────────────────────

echo "===== MARINE Single-Cell Pipeline ====="
echo "[LOG] Started:        $(date '+%Y-%m-%d %H:%M:%S')"
echo "[LOG] START_FROM:     ${START_FROM}"
echo "[LOG] HOMEDIR:        ${HOMEDIR}"
echo "[LOG] SAMPLESHEET:    ${SAMPLESHEET}"
echo "[LOG] GENOME_FA:      ${GENOME_FA}"
echo "[LOG] ANNOTATION_BED: ${ANNOTATION_BED}"
echo "[LOG] DBSNP_BED:      ${DBSNP_BED}"
echo "[LOG] CONTAINER:      ${CONTAINER}"
echo "[LOG] STRANDEDNESS:   ${STRANDEDNESS}"
echo "[LOG] MIN_COUNT:      ${MIN_COUNT} | MAX_FRAC: ${MAX_FRAC}"
if [[ "${START_FROM}" == "fastq" ]]; then
    echo "[LOG] TRANSCRIPTOME:  ${TRANSCRIPTOME_REF}"
    echo "[LOG] CELLRANGER_OUT: ${CELLRANGER_OUTDIR}"
fi
echo "======================================="

# Use process substitution (not a pipe) so ALL_JIDS persists in this shell
declare -a ALL_JIDS=()

# ── FASTQ MODE: CellRanger → MARINE (chained per sample) ─────────────────────

if [[ "${START_FROM}" == "fastq" ]]; then

    echo ""
    echo "[LOG] Step 0: Submitting CellRanger jobs (one per sample)"

    while IFS=',' read -r samplename fastq_dir; do

        samplename=$(echo "${samplename}" | tr -d '\r')
        fastq_dir=$(echo "${fastq_dir}"   | tr -d '\r')

        # Submit CellRanger job
        CR_JID=$(sbatch \
            --job-name="sc_cr_${samplename}" \
            --mail-user="${MAIL_USER}" \
            --export=ALL,SAMPLENAME="${samplename}",FASTQ_DIR="${fastq_dir}",CELLRANGER_OUTDIR="${CELLRANGER_OUTDIR}",TRANSCRIPTOME_REF="${TRANSCRIPTOME_REF}" \
            "00_cellranger_sc.sh" \
            | awk '{print $NF}')

        # CellRanger always writes to this predictable location — safe to define
        # before the job runs since the paths are deterministic
        derived_bam="${CELLRANGER_OUTDIR}/${samplename}/outs/possorted_genome_bam.bam"
        derived_matrixdir="${CELLRANGER_OUTDIR}/${samplename}/outs/filtered_feature_bc_matrix"
        barcodes="${derived_matrixdir}/barcodes.tsv.gz"
        sample_outdir="${MARINE_OUTDIR}/${samplename}"

        # Submit MARINE — runs only after CellRanger succeeds
        JID=$(sbatch \
            --job-name="sc_marine_${samplename}" \
            --mail-user="${MAIL_USER}" \
            --dependency="afterok:${CR_JID}" \
            --export=ALL,SAMPLENAME="${samplename}",USE_BAM="${derived_bam}",BARCODES="${barcodes}",SAMPLE_OUTDIR="${sample_outdir}",MARINE_OUTDIR="${MARINE_OUTDIR}",GENOME_FA="${GENOME_FA}",CONTAINER="${CONTAINER}",MARINE_PY="${MARINE_PY}",ANNOTATION_BED="${ANNOTATION_BED}",BARCODE_TAG="${BARCODE_TAG}",STRANDEDNESS="${STRANDEDNESS}",MIN_READ_QUALITY="${MIN_READ_QUALITY}",MIN_BASE_QUALITY="${MIN_BASE_QUALITY}" \
            "01_run_marine_sc.sh" \
            | awk '{print $NF}')

        echo "  ${samplename}: CellRanger → ${CR_JID}, MARINE → ${JID} (after CellRanger)"
        ALL_JIDS+=("${JID}")

    done < <(tail -n +2 "${SAMPLESHEET}")

# ── BAM MODE: MARINE only ─────────────────────────────────────────────────────

elif [[ "${START_FROM}" == "bam" ]]; then

    echo ""
    echo "[LOG] Step 1: Submitting MARINE jobs (one per sample)"

    while IFS=',' read -r samplename bamfile matrixdir; do

        samplename=$(echo "${samplename}" | tr -d '\r')
        bamfile=$(echo "${bamfile}"       | tr -d '\r')
        matrixdir=$(echo "${matrixdir}"   | tr -d '\r')

        barcodes="${matrixdir}/barcodes.tsv.gz"
        sample_outdir="${MARINE_OUTDIR}/${samplename}"

        JID=$(sbatch \
            --job-name="sc_marine_${samplename}" \
            --mail-user="${MAIL_USER}" \
            --export=ALL,SAMPLENAME="${samplename}",USE_BAM="${bamfile}",BARCODES="${barcodes}",SAMPLE_OUTDIR="${sample_outdir}",MARINE_OUTDIR="${MARINE_OUTDIR}",GENOME_FA="${GENOME_FA}",CONTAINER="${CONTAINER}",MARINE_PY="${MARINE_PY}",ANNOTATION_BED="${ANNOTATION_BED}",BARCODE_TAG="${BARCODE_TAG}",STRANDEDNESS="${STRANDEDNESS}",MIN_READ_QUALITY="${MIN_READ_QUALITY}",MIN_BASE_QUALITY="${MIN_BASE_QUALITY}" \
            "01_run_marine_sc.sh" \
            | awk '{print $NF}')

        echo "  ${samplename} → job ${JID}"
        ALL_JIDS+=("${JID}")

    done < <(tail -n +2 "${SAMPLESHEET}")

fi

# ── SUMMARY ────────────────────────────────────────────────────────────────────

echo ""
echo "[LOG] All jobs submitted at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[INFO] Monitor with:  squeue -u ${USER}"
echo "[INFO] Logs in:       ${LOG_DIR}/"
echo "[INFO] MARINE job IDs: ${ALL_JIDS[*]}"
echo ""
echo "===== NEXT STEPS ====="
echo "[INFO] Once all MARINE jobs finish, comment out PART 1 above and uncomment PART 2 below,"
echo "[INFO] then re-run this script to filter and normalize all samples."
echo "======================"

# ══════════════════════════════════════════════════════════════════════════════
# ── PART 2: FILTER + NORMALIZE (head node) ────────────────────────────────────
#
# HOW TO USE:
#   1. Comment out the PART 1 block above (from "declare -a ALL_JIDS" to the
#      "NEXT STEPS" echo block, inclusive).
#   2. Uncomment the block below (remove the leading #).
#   3. Re-run:  bash 00_marine_sc_pipeline.sh
#
# All project variables (SAMPLESHEET, MARINE_OUTDIR, FILTER_OUTDIR, etc.)
# are already defined above — no need to rewrite them.
# ══════════════════════════════════════════════════════════════════════════════

# echo "===== MARINE SC Pipeline — Part 2: Filter + Normalize ====="
# echo "[LOG] Started:          $(date '+%Y-%m-%d %H:%M:%S')"
# echo "[LOG] SAMPLESHEET:      ${SAMPLESHEET}"
# echo "[LOG] MARINE_OUTDIR:    ${MARINE_OUTDIR}"
# echo "[LOG] FILTER_OUTDIR:    ${FILTER_OUTDIR}"
# echo "[LOG] NORMALIZE_OUTDIR: ${NORMALIZE_OUTDIR}"
# echo "[LOG] DBSNP_BED:        ${DBSNP_BED}"
# echo "[LOG] ANNOTATION_BED:   ${ANNOTATION_BED}"
# echo "[LOG] MIN_COUNT:        ${MIN_COUNT} | MAX_FRAC: ${MAX_FRAC}"
# echo "[LOG] EDIT_TYPE:        ${EDIT_TYPE}"
# echo "============================================================="
#
# echo ""
# echo "[LOG] Step 2: Filtering MARINE edits..."
# bash "02_filter_marine_edits_sc.sh"
#
# echo ""
# echo "[LOG] Step 3: Normalizing MARINE edits..."
# bash "03_normalize_marine_edits_umi_sc.sh"
#
# echo ""
# echo "[DONE] Part 2 complete at $(date '+%Y-%m-%d %H:%M:%S')"
