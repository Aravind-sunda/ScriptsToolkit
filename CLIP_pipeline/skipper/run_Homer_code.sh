#!/bin/bash

#SBATCH --job-name=custom_Homer
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --output=slurm_%u_%x_%j.log

# --------------------------------------------------------------------------------
# NOTES
# remember to change the rna flag inside the script of if you are changing the type of motif you want to search for (rna vs dna)
# also remember to remove vertebrates flag if you are running in rna mode
# you can add more flags to the HOMER_EXTRA_FLAGS variable in the run_homer.sh script
# --------------------------------------------------------------------------------
# HELPERS
# 1. Filter unnecessary rows from finemapped windows file (e.g. remove rows with TERT in gene_name column)
# delta_2to4_finemapped_windows = pd.read_csv("/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/skipper/results/output/secondary_results/finemapping/mapped_sites/Delta_2to4.finemapped_windows.annotated.tsv", sep="\t")
# delta_2to4_finemapped_windows_no_TERT = delta_2to4_finemapped_windows[~delta_2to4_finemapped_windows["gene_name"].str.contains("TERT")].copy()
# # convert delta_2to4_finemapped_windows_no_TERT to a bed format with 6 columns 
# delta_2to4_finemapped_windows_no_TERT_bed = delta_2to4_finemapped_windows_no_TERT[["chrom","start",	"end",	"name","score","strand","thickStart","thickEnd","itemRgb"]].copy()
# # save to new folder
# outputdir = "/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/skipper/analysis/delta_2to4_analysis"
# delta_2to4_finemapped_windows_no_TERT_bed.to_csv(f"{outputdir}/delta_2to4_finemapped_windows_no_TERT.bed.gz", index=False, sep="\t", compression="gzip")
# # print number of rows removed
# num_removed = delta_2to4_finemapped_windows.shape[0] - delta_2to4_finemapped_windows_no_TERT.shape[0]
# print(f"Removed {num_removed} rows with TERT in gene_name column")
# --------------------------------------------------------------------------------

module load mamba
mamba activate 
mamba activate homer

GENOME="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/CLIP_pipeline/reference/Encode/hg38_star/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
OUTDIR="/path/to/outdir/"
# clear output directory to avoid conflicts with previous runs 

rm -rf $OUTDIR

/home/tmhaxs421/brannanlab/tmhaxs421/CLIP/TERT_DELTA/scripts_results_skipper/run_homer.sh --experiment-label <experiment_label> \
    --genome $GENOME \
    --finemapped <path/to/skipper/results/output/secondary_results/finemapping/mapped_sites/> \
    --background <path/to/results/output/homer/region_matched_background/fixed/> \
    --outdir $OUTDIR \
    --preparsed-dir $OUTDIR/preparsed