#!/bin/bash

#SBATCH --partition=bigmemq
#SBATCH --job-name=flare-snakemake-false_CT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=1300G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%u_%x_%j.log

module load python3
module load python-site-packages
module load mamba
module load snakemake
# module load flare/2024.03

snakemake --snakefile /home/tmhaxs421/brannanlab/VS_share/FLARE/workflow_sailor/Snakefile \
    --configfile /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/json/t1_shCTRL_sailor_input_false_CT.json \
    --verbose --use-singularity \
    --cores 90 \
    --singularity-args '--bind /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR --bind /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/dbsnp142-mm10 --bind /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/json --bind /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/sailor_op_false_CT --bind /home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/subset-bam --bind /home/tmhaxs421/brannanlab/VS_share/FLARE/workflow_sailor --bind /home/tmhaxs421/brannanlab/10x_genomics/Mouse_genome_10x/mRuby3_mm10/mRUBY3_mm10/fasta' \
    --latency-wait 600 \
    --retries 2