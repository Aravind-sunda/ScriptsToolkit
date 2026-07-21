#!/bin/bash

#SBATCH --partition=defq
#SBATCH --job-name=homer_findMotifs_11
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%u_%x_%j.log

# ============================================================================================================================================================
# HOMER motif analysis on reproducible IDR peaks from step 10
#
# Mirrors the Skipper approach:
#   - findMotifsGenome.pl with -rna -size given -len 5,6,7,8,9 -nlen 1 -S 20
#   - Background: either HOMER's built-in random genomic (default) or a custom matched BED
#   - Input: reproducible peaks BED from step 10 IDR
# ============================================================================================================================================================
module load mamba
mamba activate
mamba activate homer

export PATH="$PATH:/home/tmhaxs421/brannanlab/tmhaxs421/applications/HOMER/bin"

HOMER="/home/tmhaxs421/brannanlab/tmhaxs421/applications/HOMER/bin/findMotifsGenome.pl"
GENOME="hg38"   # HOMER has hg38 installed at HOMER/data/genomes/hg38
PREPARSED_DIR="/condo/brannanlab/tmhaxs421/homer_preparsed"

# ---- SET THESE to match your step 10 output ----
HOMEDIR=""
rep1_clip_basename=""
rep2_clip_basename=""

IDR_DIR="$HOMEDIR/10_idr"
OUTDIR="$HOMEDIR/11_homer"
mkdir -p "$OUTDIR" "$PREPARSED_DIR"

# Reproducible peaks BED from step 10 (chr, start, end, name, score, strand)
PEAKS_BED="$IDR_DIR/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.bed"

# ============================================================================================================================================================
# Convert reproducible peaks to HOMER's expected BED5 format:
#   col1=name  col2=chr  col3=start(1-based)  col4=end  col5=strand
# CLIPper BED is 0-based; HOMER expects 1-based start → add 1 to col2
# BED format is 0-based half-open
# HOMER's peak format is 1-based closed
# ============================================================================================================================================================
HOMER_INPUT="$OUTDIR/homer_fg.bed"

awk 'BEGIN{OFS="\t"} {
    name = ($4 != "" && $4 != ".") ? $4 : $1":"$2"-"$3
    strand = ($6 != "" && $6 != ".") ? $6 : "+"
    print name, $1, $2+1, $3, strand
}' "$PEAKS_BED" > "$HOMER_INPUT"

echo "Foreground peaks: $(wc -l < "$HOMER_INPUT") regions"

# ============================================================================================================================================================
# Run findMotifsGenome.pl
#
# Key flags matching Skipper:
#   -rna         : treat sequences as RNA (T→U), use RNA motif databases
#   -size given  : use the actual peak coordinates (not a fixed window)
#   -len 5,6,7,8,9 : discover motifs of lengths 5–9 nt
#   -nlen 1      : normalize by 1-nt background frequencies
#   -S 20        : report top 20 motifs
#   -p           : threads
#   -preparsedDir: cache parsed genome sequences (speeds up reruns)
#   -bg          : (optional) custom background BED — see BACKGROUND section below
# ============================================================================================================================================================

"$HOMER" \
    "$HOMER_INPUT" \
    "$GENOME" \
    "$OUTDIR/motifs" \
    -preparsedDir "$PREPARSED_DIR" \
    -rna \
    -size given \
    -len 5,6,7,8,9 \
    -nlen 1 \
    -S 20 \
    -p "$SLURM_CPUS_PER_TASK"

echo "HOMER motif analysis done. Results in $OUTDIR/motifs/"

# ============================================================================================================================================================
# OPTION A: Re-run HOMER with a matched background (non-reproducible CLIPper peaks)
#
# Background = all peaks detected by CLIPper in either replicate (step 09 compressed beds)
#              MINUS the reproducible IDR peaks used as foreground.
# These are regions the RBP signal was detected in but did not pass reproducibility
# filtering — they live in the same transcriptomic space as true binding sites,
# which avoids the GC/feature-composition bias of HOMER's random genomic background.
#
# Output goes to $OUTDIR/motifs_matched_bg/ for direct comparison with $OUTDIR/motifs/
# ============================================================================================================================================================


module load mamba
mamba activate
mamba activate bioinformatics   # needs bedtools

WORKDIR_09="$HOMEDIR/09_normCompressPeaks"
rep1_compressed="$WORKDIR_09/$rep1_clip_basename.peakClusters.normed.compressed.bed"
rep2_compressed="$WORKDIR_09/$rep2_clip_basename.peakClusters.normed.compressed.bed"

BG_CANDIDATES="$OUTDIR/background_candidates.bed"
HOMER_BG="$OUTDIR/homer_bg.bed"

# Step 1: merge peaks from both reps (strand-aware), subtract reproducible peaks.
#
# Pipeline:
#   cat both compressed beds → sort → bedtools merge -s (strand-aware merge)
#     → output: chr, start, end, strand (4 cols)
#   awk: pad to 6 cols so bedtools subtract -s can find strand in col6
#     → chr, start, end, name, len, strand
#   bedtools subtract -s: remove any region overlapping a reproducible peak
#     → background candidates, 6 cols, strand in col6
cat "$rep1_compressed" "$rep2_compressed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -s -i stdin -c 6 -o distinct \
    | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $1":"$2"-"$3, $3-$2, $4}' \
    | bedtools subtract -s -a stdin -b "$PEAKS_BED" \
    > "$BG_CANDIDATES"

n_bg=$(wc -l < "$BG_CANDIDATES")
echo "Background candidates: $n_bg regions"

if [ "$n_bg" -eq 0 ]; then
    echo "ERROR: no background candidates found — check paths to compressed beds"
    exit 1
fi

mamba deactivate 
mamba activate homer

# Step 2: convert to HOMER BED5 (name, chr, start+1, end, strand)
# col1=chr col2=start(0-based) col3=end col4=name col5=len col6=strand
awk 'BEGIN{OFS="\t"}{print $4, $1, $2+1, $3, $6}' "$BG_CANDIDATES" > "$HOMER_BG"

echo "Running HOMER with matched background ($n_bg regions)..."

"$HOMER" \
    "$HOMER_INPUT" \
    "$GENOME" \
    "$OUTDIR/motifs_matched_bg" \
    -preparsedDir "$PREPARSED_DIR" \
    -rna \
    -size given \
    -len 5,6,7,8,9 \
    -nlen 1 \
    -S 20 \
    -p "$SLURM_CPUS_PER_TASK" \
    -bg "$HOMER_BG"

echo "HOMER matched-background run done. Results in $OUTDIR/motifs_matched_bg/"

# ============================================================================================================================================================
# OPTIONAL: Run per-annotation-category (like Skipper's region-stratified approach)
# Uncomment and set ANNOTATED_BED to the step 10 annotated output if you want
# separate motif runs per genomic feature (CDS, UTR3, intron, etc.)
#
# ANNOTATED_BED="$IDR_DIR/${rep1_clip_basename}_${rep2_clip_basename}_reproducible_peaks.sorted.annotated.bed"
#
# for region in UTR3 CDS intron UTR5 ncRNA; do
#     region_bed="$OUTDIR/homer_fg_${region}.bed"
#     awk -v r="$region" 'BEGIN{OFS="\t"} $0 ~ r {
#         name = ($4 != "" && $4 != ".") ? $4 : $1":"$2"-"$3
#         strand = ($6 != "" && $6 != ".") ? $6 : "+"
#         print name, $1, $2+1, $3, strand
#     }' "$ANNOTATED_BED" > "$region_bed"
#
#     n=$(wc -l < "$region_bed")
#     if [ "$n" -lt 10 ]; then
#         echo "Skipping $region — only $n peaks"
#         continue
#     fi
#
#     echo "Running HOMER on $n $region peaks..."
#     "$HOMER" \
#         "$region_bed" \
#         "$GENOME" \
#         "$OUTDIR/motifs_${region}" \
#         -preparsedDir "$PREPARSED_DIR" \
#         -rna \
#         -size given \
#         -len 5,6,7,8,9 \
#         -nlen 1 \
#         -S 20 \
#         -p "$SLURM_CPUS_PER_TASK"
# done
# ============================================================================================================================================================
