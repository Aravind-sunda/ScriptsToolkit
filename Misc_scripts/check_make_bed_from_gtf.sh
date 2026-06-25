#!/bin/bash

#SBATCH --job-name=make_bed_from_gtf
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --output=slurm_%u_%x_%j.log

# this code tries 2 different methods and checks for consistency between the 2 methods. the first method uses the bedops gtf2bed convert tool, and the second method uses awk to parse the gtf file. both methods should produce the same bed file after adding mRUBY3 to both files manually.
# usign the second method moving forward since it is faster


# module load mamba 
# mamba activate 
# mamba activate bioinformatics

# # using the bedops gtf2bed convert tool

# gunzip -k genes.gtf.gz

# convert2bed -i gtf -m 36G  < genes.gtf  > genes_gtf2bed.bed

# awk -F'\t' 'BEGIN{OFS="\t"}
# $8=="gene"{
#   name=$4                                      # fallback to gene_id
#   if (match($0, /gene_name "[^"]*"/))
#       name=substr($0, RSTART+11, RLENGTH-12)
#   print $1,$2,$3,name,".",$6
# }' genes_gtf2bed.bed > genes_gtf2bed.bed6


# awk -F'\t' '$3=="gene" {
#   match($9, /gene_name "([^"]+)"/, a);
#   match($9, /gene_type "([^"]+)"/, b);
#   print $1 "\t" ($4-1) "\t" $5 "\t" a[1] "\t" b[1] "\t" $7
# }' genes.gtf > genes_awk.bed6

# # add mruby3 to the bed file at the end of the line
# # mRUBY3	1	711	mRUBY3	protein_coding	+
# # add mRUBY3 to the end of the bed file

# echo -e "mRUBY3\t1\t711\tmRUBY3\t.\t+" >> genes_gtf2bed.bed6
# echo -e "mRUBY3\t1\t711\tmRUBY3\tprotein_coding\t+" >> genes_awk.bed6


# rm -rf genes.gtf

# # sanity check to make sure that 

# bedtools intersect -v -a genes_awk.bed6 -b genes_gtf2bed.bed6 > genes_awk_no_overlap.bed

# This test passed, so the 2 files are the same after adding mRUBY3 to both files. manually
# code to  run