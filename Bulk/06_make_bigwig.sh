#!/bin/bash
#SBATCH --partition=defq
#SBATCH --job-name=make_bigwig
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=96:00:00
#SBATCH --output=slurm_%u_%x_%j.log

echo "Starting job at $(date '+%Y-%m-%d %H:%M:%S')"

WORKING_DIR=""
INPUT_BAM_DIR="$WORKING_DIR/star_output"

# CHROM SIZES: path to chromosome sizes file (tab-separated: chrom\tsize)
# Generate with: samtools view -H sample.bam | awk '/^@SQ/{split($2,a,":");split($3,b,":");print a[2]"\t"b[2]}'
# Or from genome fasta: samtools faidx genome.fa && cut -f1,2 genome.fa.fai > chrom.sizes
CHROM_SIZES=""

# STRANDEDNESS: set to "unstranded", "forward", or "reverse"
# Determine from 03_bam_qc.sh infer_experiment.py output
# forward  = dUTP second-strand (most TruSeq Stranded kits)
# reverse  = first-strand (uncommon)
# unstranded = no strand info
STRANDEDNESS="reverse"

# SUFFIX
BAM_SUFFIX="Aligned.sortedByCoord.out.bam"

# OUTPUT
OUTDIR="$WORKING_DIR/bigwig"
mkdir -p "$OUTDIR"

# MODULES
module load mamba
mamba activate
mamba activate bioinformatics

# ========================================================================
# HELPER: compute CPM scale factor from mapped reads in a BAM
# scale = 1,000,000 / total_mapped_reads
# bedtools genomecov -scale multiplies each bin count by this value,
# producing CPM-normalised coverage comparable across all samples.
# ========================================================================
get_scale_factor() {
    local bam="$1"
    local mapped
    mapped=$(samtools flagstat "$bam" | awk '/mapped \(/ {print $1; exit}')
    if [[ -z "$mapped" || "$mapped" -eq 0 ]]; then
        echo "ERROR: could not get mapped read count for $bam" >&2
        echo "1"
    else
        awk -v m="$mapped" 'BEGIN {printf "%.10f", 1000000/m}'
    fi
}

# ========================================================================
# HELPER: bedgraph -> sorted -> bigWig, then clean up bedgraph
# ========================================================================
make_bigwig() {
    local prefix="$1"   # full output path prefix (no extension)
    local bam="$2"
    local scale="$3"
    local strand="$4"   # "+" | "-" | "" (empty = combined/unstranded)

    local strand_arg=""
    [[ -n "$strand" ]] && strand_arg="-strand $strand"

    bedtools genomecov \
        -ibam "$bam" \
        -bg \
        -scale "$scale" \
        $strand_arg \
        | LC_ALL=C sort --parallel="$SLURM_CPUS_PER_TASK" -k1,1 -k2,2n \
        > "${prefix}.bedgraph"

    bedGraphToBigWig \
        "${prefix}.bedgraph" \
        "$CHROM_SIZES" \
        "${prefix}.bigWig"

    rm "${prefix}.bedgraph"
}

# ========================================================================
# MAKE BIGWIGS
#
# All samples:      combined (no strand filter) bigwig
# Stranded samples: additionally fwd (+ strand) and rev (- strand) bigwigs
#
# Which strand-split track is "sense" depends on your library prep:
#   forward-stranded (e.g. ScriptSeq):        fwd = sense
#   reverse-stranded (e.g. dUTP/TruSeq Stranded): rev = sense
# ========================================================================

for bam in "$INPUT_BAM_DIR"/*"$BAM_SUFFIX"; do

    sample=$(basename "$bam" ".$BAM_SUFFIX")
    echo "[INFO] Processing $sample"

    scale=$(get_scale_factor "$bam")
    echo "[INFO] $sample: scale factor = $scale (CPM)"

    # Combined track for all samples
    make_bigwig "$OUTDIR/${sample}.combined.CPM" "$bam" "$scale" ""

    # Strand-split tracks for stranded libraries
    if [[ "$STRANDEDNESS" != "unstranded" ]]; then
        make_bigwig "$OUTDIR/${sample}.fwd.CPM" "$bam" "$scale" "+"
        make_bigwig "$OUTDIR/${sample}.rev.CPM" "$bam" "$scale" "-"
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done: $sample"
done

echo "All bigwigs complete."
