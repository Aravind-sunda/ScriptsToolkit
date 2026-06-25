#!/bin/bash

# Submit one sbatch job per sample row in the samplesheet CSV.
# Usage: bash submit_09_10_jobs.sh [samplesheet.csv]
#
# CSV columns (with header):
#   rep1_input_bam, rep2_input_bam, rep1_clip_bam, rep2_clip_bam,
#   rep1_peaks_bed, rep2_peaks_bed, outputdir, homedir
#
# All BAM/BED paths must be full absolute paths.
# homedir: the per-sample working root where 09_normCompressPeaks/, 09a_..., 10_idr/ are created.

SAMPLESHEET="/home/tmhaxs421/brannanlab/tmhaxs421/CLIP/TERT_Frags/scripts/samplesheet.csv"
SCRIPT="/home/tmhaxs421/brannanlab/tmhaxs421/CLIP/TERT_Frags/scripts/09_10_normCompressPeaks_IDR_Annotate.sh"
HOMEDIR="/home/tmhaxs421/brannanlab/tmhaxs421/CLIP/TERT_Frags"


if [[ ! -f "$SAMPLESHEET" ]]; then
    echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2
    exit 1
fi

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: analysis script not found: $SCRIPT" >&2
    exit 1
fi

submitted=0
skipped=0

while IFS=',' read -r sample rep1_input_bam rep2_input_bam rep1_clip_bam rep2_clip_bam \
                         rep1_peaks_bed rep2_peaks_bed; do

    # Skip header row
    [[ "$sample" == "sample" ]] && continue

    # Skip blank lines
    [[ -z "${rep1_input_bam// }" ]] && continue

    # Validate required fields
    missing=()
    [[ -z "$sample"         ]] && missing+=("sample")
    [[ -z "$rep1_input_bam" ]] && missing+=("rep1_input_bam")
    [[ -z "$rep2_input_bam" ]] && missing+=("rep2_input_bam")
    [[ -z "$rep1_clip_bam"  ]] && missing+=("rep1_clip_bam")
    [[ -z "$rep2_clip_bam"  ]] && missing+=("rep2_clip_bam")
    [[ -z "$rep1_peaks_bed" ]] && missing+=("rep1_peaks_bed")
    [[ -z "$rep2_peaks_bed" ]] && missing+=("rep2_peaks_bed")

    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "WARNING: skipping incomplete row (missing: ${missing[*]})" \
             "— clip bam: ${rep1_clip_bam:-<empty>}" >&2
        skipped=$((skipped + 1))
        continue
    fi


    job_name="${sample}_idr"
    
    # adding folder paths to the file names in the samplesheet
    # appending paths to each variable since they will only contain the file names and not the full path. This is because the samplesheet is easier to read with just the file names and not the full paths
    rep1_input_bam_full="${HOMEDIR}/06_dedup_sort/${rep1_input_bam}"
    rep2_input_bam_full="${HOMEDIR}/06_dedup_sort/${rep2_input_bam}"
    rep1_clip_bam_full="${HOMEDIR}/06_dedup_sort/${rep1_clip_bam}"
    rep2_clip_bam_full="${HOMEDIR}/06_dedup_sort/${rep2_clip_bam}"
    rep1_peaks_bed_full="${HOMEDIR}/07_clipper_output/${rep1_peaks_bed}"
    rep2_peaks_bed_full="${HOMEDIR}/07_clipper_output/${rep2_peaks_bed}"

    # so that folders are created correctly in 09_10_normCompressPeaks_IDR_Annotate.sh
    homedir_full="${HOMEDIR}/${sample}"

    # check if all files are present
    for file in "$rep1_input_bam_full" "$rep2_input_bam_full" "$rep1_clip_bam_full" "$rep2_clip_bam_full" "$rep1_peaks_bed_full" "$rep2_peaks_bed_full"; do
        if [[ ! -f "$file" ]]; then
            echo "ERROR: file not found: $file — skipping $sample" >&2
            skipped=$((skipped + 1))
            continue 2   # continue the outer while loop
        fi
    done

    echo "Submitting: $job_name"
    sbatch \
        --job-name="$job_name" \
        --output="slurm_%u_%x_%j.log" \
        "$SCRIPT" \
        --rep1-input-bam "$rep1_input_bam_full" \
        --rep2-input-bam "$rep2_input_bam_full" \
        --rep1-clip-bam  "$rep1_clip_bam_full"  \
        --rep2-clip-bam  "$rep2_clip_bam_full"  \
        --rep1-peaks-bed "$rep1_peaks_bed_full" \
        --rep2-peaks-bed "$rep2_peaks_bed_full" \
        --homedir        "$homedir_full"

    submitted=$((submitted + 1))

done < "$SAMPLESHEET"

echo ""
echo "Submitted: $submitted  |  Skipped (incomplete): $skipped"
