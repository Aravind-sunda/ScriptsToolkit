#!/bin/bash
#SBATCH --job-name=qc_read_counts
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=2:00:00
#SBATCH --output=slurm_%u_%x_%j.log

# Usage: bash qc_read_counts.sh [HOMEDIR]
# If HOMEDIR is not passed as $1, it falls back to the $HOMEDIR environment variable
# (which is exported in 00_CLIP_Pipeline.sh).
#
# Output files written to $HOMEDIR/:
#   qc_read_counts_summary.tsv     -- read counts at each pipeline stage per sample
#   qc_unmapped_dedup_removed.tsv  -- unmapped reads (genome) + reads removed by dedup

HOMEDIR="${1:-$HOMEDIR}"

if [ -z "$HOMEDIR" ]; then
    echo "ERROR: HOMEDIR not set. Run as: bash qc_read_counts.sh /path/to/HOMEDIR"
    exit 1
fi

if [ ! -d "$HOMEDIR" ]; then
    echo "ERROR: HOMEDIR directory not found: $HOMEDIR"
    exit 1
fi

OUTPUT1="$HOMEDIR/qc_read_counts_summary.tsv"
OUTPUT2="$HOMEDIR/qc_unmapped_dedup_removed.tsv"

echo "---------------------"
echo "[LOG] Starting QC Read Count Summary"
echo "[LOG] HOMEDIR: $HOMEDIR"
echo "---------------------"

# ── File 1 header ────────────────────────────────────────────────────────────
echo -e "Sample\tRaw_reads\tAfter_trim\tPct_surviving_trim\tAfter_repbase_filter\tRepbase_dropped\tPct_repbase_dropped\tAfter_genome_unique_map\tPct_genome_mapped\tAfter_dedup\tPct_dedup_remaining\tPct_usable" \
    > "$OUTPUT1"

# ── File 2 header ────────────────────────────────────────────────────────────
echo -e "Sample\tReads_removed_dedup\tPct_removed_dedup\tGenome_unmapped_mismatch\tGenome_unmapped_tooshort\tGenome_unmapped_other\tTotal_genome_unmapped\tPct_genome_unmapped" \
    > "$OUTPUT2"

# ── Main loop over samples (ordered alphabetically by UMI metrics filename) ──
for umi_metrics in $(ls -1 "$HOMEDIR/01_UMI_clip"/*.log.---.--.metrics | sort); do

    sample=$(basename "$umi_metrics" .log.---.--.metrics)
    echo "[LOG] Processing: $sample"

    # ------------------------------------------------------------------
    # STEP 1 – Raw reads
    # Source: umi_tools extract log (01_UMI_clip/*.log.---.--.metrics)
    # Line format:  "2026-... INFO Input Reads: 42455957"
    # ------------------------------------------------------------------
    raw_reads=$(grep "Input Reads:" "$umi_metrics" | awk '{print $NF}')

    # ------------------------------------------------------------------
    # STEP 2 – Reads after adapter trimming
    # ------------------------------------------------------------------

    # Method A (ACTIVE): count lines in the final sorted trimmed fastq / 4.
    # Uses wc -l (fast, reads only newlines) then integer division by 4.
    # Works for all runs regardless of whether per-sample cutadapt logs exist.
    trim_fq="$HOMEDIR/02_cutadapt2/${sample}.umi.fqtrtr.sorted.fq"
    after_trim=$(( $(wc -l < "$trim_fq") / 4 ))

    # Method B (COMMENTED): parse the per-sample cutadapt2 metrics file.
    # Requires 02_trimAdapters_fastqSort.sh to have the per-sample redirect fix.
    # Line format:  "Reads written (passing filters):    32,136,312 (98.9%)"
    # trim2_metrics="$HOMEDIR/02_cutadapt2/${sample}.IP.umi.r1.fqTrTr.metrics"
    # after_trim=$(grep "Reads written" "$trim2_metrics" | awk '{gsub(/,/,""); print $(NF-1)}')

    # ------------------------------------------------------------------
    # STEP 3 – Reads surviving repbase filter
    # Source: genome STAR log (04_star_hg38_rep_sam/*.genomeMapped.Log.final.out)
    # "Number of input reads" in the genome log = reads that did NOT map to
    # repbase and therefore proceed to genome alignment.
    # ------------------------------------------------------------------
    genome_log="$HOMEDIR/04_star_hg38_rep_sam/${sample}.repbaseMapped.genomeMapped.Log.final.out"
    after_repbase=$(grep "Number of input reads" "$genome_log" | awk '{print $NF}')
    repbase_dropped=$((after_trim - after_repbase))

    # ------------------------------------------------------------------
    # STEP 4 – Reads uniquely mapped to genome
    # Source: same genome STAR log
    # ------------------------------------------------------------------
    genome_unique=$(grep "Uniquely mapped reads number" "$genome_log" | awk '{print $NF}')

    # Genome unmapped breakdown (used in File 2)
    genome_unmapped_mismatch=$(grep "Number of reads unmapped: too many mismatches" "$genome_log" | awk '{print $NF}')
    genome_unmapped_tooshort=$(grep "Number of reads unmapped: too short"           "$genome_log" | awk '{print $NF}')
    genome_unmapped_other=$(   grep "Number of reads unmapped: other"               "$genome_log" | awk '{print $NF}')
    total_genome_unmapped=$((genome_unmapped_mismatch + genome_unmapped_tooshort + genome_unmapped_other))

    # ------------------------------------------------------------------
    # STEP 5/6 – Reads after deduplication
    # ------------------------------------------------------------------

    # Method A (ACTIVE): sum total_counts_post column from per_umi.tsv.
    # This file is always produced by umi_tools dedup --output-stats.
    # Column layout: UMI | median_pre | times_pre | total_pre | median_post | times_post | total_post
    per_umi_tsv="$HOMEDIR/06_dedup_sort/${sample}.genome-mappedSoSo_per_umi.tsv"
    after_dedup=$(awk 'NR>1 {sum+=$7} END{print sum}' "$per_umi_tsv")

    # Method B (COMMENTED): parse the umi_tools dedup log file.
    # Requires 06_dedup_sort.sh to include:  --log $OUTPUTDIR/$outfile.genome-mappedSoSo.dedup.log
    # Line format:  "... INFO Reads in: 14272817, reads out: 14272817"
    # dedup_log="$HOMEDIR/06_dedup_sort/${sample}.genome-mappedSoSo.dedup.log"
    # after_dedup=$(grep "reads out:" "$dedup_log" | awk -F'reads out: ' '{print $2}')

    # ------------------------------------------------------------------
    # Percentages (all calculated with awk to support decimal output)
    # ------------------------------------------------------------------
    pct_surviving_trim=$( awk "BEGIN{printf \"%.2f\", ($after_trim/$raw_reads)*100}")
    pct_repbase_dropped=$(awk "BEGIN{printf \"%.2f\", ($repbase_dropped/$after_trim)*100}")
    pct_genome_mapped=$(  awk "BEGIN{printf \"%.2f\", ($genome_unique/$after_repbase)*100}")
    pct_dedup_remaining=$(awk "BEGIN{printf \"%.2f\", ($after_dedup/$genome_unique)*100}")
    pct_usable=$(         awk "BEGIN{printf \"%.2f\", ($after_dedup/$raw_reads)*100}")

    reads_removed_dedup=$((genome_unique - after_dedup))
    pct_removed_dedup=$(  awk "BEGIN{printf \"%.2f\", ($reads_removed_dedup/$genome_unique)*100}")
    pct_genome_unmapped=$(awk "BEGIN{printf \"%.2f\", ($total_genome_unmapped/$after_repbase)*100}")

    # ------------------------------------------------------------------
    # Write rows
    # ------------------------------------------------------------------
    printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%d\t%s\t%s\n" \
        "$sample"           \
        "$raw_reads"        \
        "$after_trim"       \
        "$pct_surviving_trim"   \
        "$after_repbase"    \
        "$repbase_dropped"  \
        "$pct_repbase_dropped"  \
        "$genome_unique"    \
        "$pct_genome_mapped"    \
        "$after_dedup"      \
        "$pct_dedup_remaining"  \
        "$pct_usable"       \
        >> "$OUTPUT1"

    printf "%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\n" \
        "$sample"                   \
        "$reads_removed_dedup"      \
        "$pct_removed_dedup"        \
        "$genome_unmapped_mismatch" \
        "$genome_unmapped_tooshort" \
        "$genome_unmapped_other"    \
        "$total_genome_unmapped"    \
        "$pct_genome_unmapped"      \
        >> "$OUTPUT2"

done

echo "---------------------"
echo "[DONE] Read counts summary : $OUTPUT1"
echo "[DONE] Unmapped/dedup report: $OUTPUT2"
echo "---------------------"
