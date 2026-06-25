#!/bin/bash

#SBATCH --job-name=downloadSRA
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --output=slurm_%u_%x_%j.log



/home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/SRAToolkit/downloadSRA.py \
--out-dir /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/skipper/data/RBFOX2_EVN_293_fasterq_usename \
--acc-list /home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/SRAToolkit/eg_SraAccList_py.tsv \
--threads 36


# NOTE: THis following run uses fastqdump in the internal code so that we can get orignal read headers from sequencing so that we can find out what the different files in a single run of SRR represents
# example use case https://www.ncbi.nlm.nih.gov/sra/SRX1563182[accn]
# orignal read headers are not saved in the GEO, so they are not present in this


# /home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/SRAToolkit/downloadSRA_fastqdump.py \
# --out-dir /home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/skipper/data/RBFOX2_EVN_293_fastqDump \
# --acc-list /home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/SRAToolkit/eg_SraAccList_py.tsv \
# --threads 36