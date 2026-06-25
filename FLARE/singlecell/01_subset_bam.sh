#!/bin/bash

#SBATCH --partition=defq
#SBATCH --job-name=subset-bam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%u_%x_%j.log


module load mamba
mamba init bash
conda activate subset-bam
module load subset-bam/1.1.0
module load samtools

#OLD CODE-----------------------------------------------------
# subset-bam --bam /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/bam/possorted_genome_bam.bam \
#     --cell-barcodes /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/barcode/barcodes.tsv \
#     --out-bam /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/subset-bam/possorted_genome_bam_subset.bam \
#     --cores 36 \
#     --log-level debug

# RUNNING CODE-----------------------------------------------------
# The following script only takes 30 minutes and 5gb of ram to run
# I think using more cores will speed up the processing of the files
# Maybe can delete the subset bam files for each cell once running 

mkdir -p /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/subset-bam/

subset_bam_per_cb.sh /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/bam/possorted_genome_bam.bam \
/home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/barcode/barcodes.tsv \
/home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/subset-bam/t1_shCTRL

echo "Done-completed processing and splitting the bam file by cell barcode"

# The following can be used in the snakemake file:
# subset_bam_per_cb.sh /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/bam/possorted_genome_bam.bam \
# /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/barcode/barcodes.tsv \
# /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/subset-bam/{sample}
# In the above command the {sample} will be replaced by the sample name in the snakemake file followed by the cell barcode.bam




