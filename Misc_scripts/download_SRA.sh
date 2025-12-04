# Load mamba and activate environment with SRA Toolkit
# module load sratoolkit

module load mamba
mamba activate bioinformatics   # or: mamba activate sra-tools

# Input arguments (set these properly)
outputDir="/home/tmhaxs421/brannanlab/tmhaxs421/CHIP_Kim/data/GEO_ChIP_Tert"        # output directory
samplesheet="/home/tmhaxs421/brannanlab/tmhaxs421/CHIP_Kim/data/GEO_ChIP_Tert/SRA_ChIP_TERT_Samples.csv"  # CSV file with SRR numbers and sample names with no header

mkdir -p $outputDir/data

# # Example: keep only Run and Sample_Name
# cut -d "," -f1,18 $runInfo | tail -n +2 > $outputDir/SRR_SampleName.tsv   # adjust "5" to the right column for Sample_Name

# # Extract SRR IDs
# cut -d "," -f 1 "$runInfo" | tail -n +2 > "$outputDir/SRR.numbers"

# Download + convert using fasterq-dump and 
# store the log in a log file in output dir

cat "$samplesheet" | parallel -j 36 --colsep ',' \
  "prefetch {1} --output-directory $outputDir/sra && \
   fasterq-dump --split-3 --threads 4 -O $outputDir/data $outputDir/sra/{1}" >> "$outputDir/download.log" 2>&1

pigz -p 36 "$outputDir"/data/*.fastq

# renaming the files to include sample names

cat "$samplesheet" | while IFS=, read -r srr sample_name; do

  if [ -f "$outputDir/data/${srr}_1.fastq.gz" ] && [ -f "$outputDir/data/${srr}_2.fastq.gz" ]; then
    mv "$outputDir/data/${srr}_1.fastq.gz" "$outputDir/data/${sample_name}_R1.fastq.gz"
    mv "$outputDir/data/${srr}_2.fastq.gz" "$outputDir/data/${sample_name}_R2.fastq.gz"

  elif [ -f "$outputDir/data/${srr}.fastq.gz" ]; then
    mv "$outputDir/data/${srr}.fastq.gz" "$outputDir/data/${sample_name}.fastq.gz"

  else
    echo "FASTQ files for ${srr} not found!"

  fi
done < "$samplesheet"


# # checking strandedness for 2 samples
# infer_experiment.py -r {params.bed} -i {input.bam} -s 1000000 > {output.txt} 2> {log}