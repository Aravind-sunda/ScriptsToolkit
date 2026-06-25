#!/bin/bash

#SBATCH --job-name=mkref
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --output=slurm_%u_%x_%j.log
# =============================================================
# IMPORTANT!!!
# RUN THIS SCRIPT ONLY ONCE SO THAT YOU DO NOT CREATE DUPLICATE ENTRIES IN THE FASTA AND GTF FILES
# LINK: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr
# Remember custom fasta headers should not have underscore, since seurat does not like underscore in the featurenames(genes) and will change them to hyphens
# =============================================================
# VARIABLES
WORKDIR="/home/tmhaxs421/brannanlab/tmhaxs421/tmhvvs2/hESC_10x/reference"
fasta_custom="/home/tmhaxs421/brannanlab/tmhaxs421/tmhvvs2/hESC_10x/reference/custom.fa"
fasta="/home/tmhaxs421/brannanlab/tmhaxs421/tmhvvs2/hESC_10x/reference/genome.fa"
gtf=/home/tmhaxs421/brannanlab/tmhaxs421/tmhvvs2/hESC_10x/reference/genes.gtf
# =============================================================
cd "$WORKDIR"
# =============================================================
# PREPROCESSING
# STEP 1: Add the new fasta to the old reference fasta
cat "$fasta_custom" >> "$fasta"

# STEP 2: Count the number of bases in the custom fasta file (multiple entries)
# cat "$fasta_custom" | grep -v "^>" | tr -d "\n" | wc -c
awk '
  /^>/ {
    if (NR>1) print hdr "\t" len
    hdr=$0
    len=0
    next
  }
  { gsub(/[ \t\r\n]/,""); len += length($0) }
  END { if (NR) print hdr "\t" len }
' "$fasta_custom"

# STEP 3: Create a custom GTF file that includes the new genes
# echo -e '<Name>\tunknown\texon\t1\t<NumberOfBases>\t.\t+\t.\tgene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";' > custom.gtf
# GFP      unknown exon 1  <len_of_GFP>      . + . gene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";
awk 'BEGIN{OFS="\t"}
  /^>/{
    if(seen) print name,"unknown","exon",1,len,".","+",".","gene_id \""name"\"; transcript_id \""name"\"; gene_name \""name"\"; gene_biotype \"protein_coding\";"
    name=substr($0,2); sub(/ .*/,"",name)   # take header, drop leading ">", keep first token
    len=0; seen=1
    next
  }
  { gsub(/[ \t\r\n]/,""); len += length($0) }
  END{
    if(seen) print name,"unknown","exon",1,len,".","+",".","gene_id \""name"\"; transcript_id \""name"\"; gene_name \""name"\"; gene_biotype \"protein_coding\";"
  }' "$fasta_custom" > custom.gtf


# STEP 4: Add the custom gtf file to the old reference gtf file
cat custom.gtf >> "$gtf"
# =============================================================
# CHECKS
grep ">" "$fasta" | tail -n 20 # gets all the headers in the fasta file
tail -n 20 "$gtf" # gets the last 20 lines of the gtf file to check if the new entries are there
# =============================================================

module load cellranger/9.0.1
genome="GRCh38"
version="2024-A"

cellranger mkref --ref-version="$version" \
    --genome="$genome" \
    --fasta="$fasta" \
    --genes="$gtf"