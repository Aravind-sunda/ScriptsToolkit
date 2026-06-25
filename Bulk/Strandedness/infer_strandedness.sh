#!/bin/bash
#SBATCH --partition=defq
#SBATCH --job-name=infer_strandedness
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=96:00:00
#SBATCH --output=slurm_%u_%x_%j.log

# ── Variables ──────────────────────────────────────────────────────────────────
BED=path/to/reference.bed
# check /home/tmhaxs421/brannanlab/tmhaxs421/reference/rseqc

BAM_FILES=(
    /path/to/sample1.bam
    /path/to/sample2.bam
    /path/to/sample3.bam
)

OUTDIR=path/to/output/directory
OUTFILE=${OUTDIR}/strandedness_results.txt
# ──────────────────────────────────────────────────────────────────────────────

# Clear output file at the start of each run
> "$OUTFILE"

module load mamba
mamba activate
mamba activate bioinformatics

for BAM in "${BAM_FILES[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)
    echo "====== $SAMPLE ======" >> "$OUTFILE"
    infer_experiment.py -i "$BAM" -r "$BED" >> "$OUTFILE" 2>&1
    echo "" >> "$OUTFILE"
done

echo "Done. Results written to $OUTFILE"
