#!/bin/bash

#SBATCH --partition=defq
#SBATCH --job-name=annotate_clipper_output
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%u_%x_%j.log



module purge
module load mamba
mamba activate
mamba activate annotator

SIF="/condo/brannanlab/tmhaxs421/applications/annotator/annotator.sif"
GTFDB="/condo/brannanlab/VS_share/Genome_indexs/GTFdb/Yeolab/gffutils_dbs/gencode.v40.annotation.gtf.db"
INPUT_DIRS=(
    "/path/to/first_input_directory"
    "/path/to/second_input_directory"
    "/path/to/third_input_directory"
)

# 3. Loop through each directory in the array
for current_dir in "${INPUT_DIRS[@]}"; do
    echo "Processing directory: $current_dir"

    # use the filename of the input bed file
    for bed_file in ${current_dir}/09_normCompressPeaks/*.compressed.bed; do
    echo "Processing file: $bed_file"
    
    input_bed="$bed_file"
    sorted_bed="${current_dir}/09a_annotated_normalized_clipper_peaks/$(basename "$input_bed" .compressed.bed).compressed.sorted.bed"
    annotated_bed="${current_dir}/09a_annotated_normalized_clipper_peaks/$(basename "$sorted_bed" .compressed.sorted.bed).compressed.sorted.annotated.bed"

    mkdir -p "${current_dir}/09a_annotated_normalized_clipper_peaks"

    # Optional but recommended: Check if the input file actually exists before running
    if [ ! -f "$input_bed" ]; then
        echo "Warning: $input_bed not found in $current_dir. Skipping."
        continue
    fi

    echo "Sorting peaks..."
    sort -k1,1 -k2,2n "$input_bed" -o "$sorted_bed"

    echo "Running annotator..."
    # Note: The --bind flag has been updated to include current_dir and OUTPUT_DIR
    apptainer exec \
        --cleanenv \
        --env PYTHONNOUSERSITE=1 \
        --env GTFDB="$GTFDB" \
        --bind "$current_dir",/condo/brannanlab \
        "$SIF" \
        annotator \
            --output "$annotated_bed" \
            --input  "$sorted_bed" \
            --gtfdb "$GTFDB"
    echo "Finished processing $current_dir"
    echo "------------------------------------------------"

    # add another additional script to create helper plots to show the total number o
    done
done

echo "All done! All annotated files are located in $OUTPUT_DIR"