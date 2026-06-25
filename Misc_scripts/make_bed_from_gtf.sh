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

module load mamba 
mamba activate 
mamba activate bioinformatics

gunzip -k genes.gtf.gz

awk -F'\t' '$3=="gene" {
  match($9, /gene_name "([^"]+)"/, a);
  match($9, /gene_type "([^"]+)"/, b);
  print $1 "\t" ($4-1) "\t" $5 "\t" a[1] "\t" b[1] "\t" $7
}' genes.gtf > genes_awk.bed6

# add any custom sequences to bed file at the end of the line
echo -e "mRUBY3\t1\t711\tmRUBY3\tprotein_coding\t+" >> genes_awk.bed6

rm -rf genes.gtf