#!/bin/bash

#SBATCH --job-name=<job_name>
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --output=slurm_%u_%x_%j.log


# 1. Load the R module (change this to match your HPC's exact R module name)
module load  R/4.5.2  # version used in my Rstudio is 4.5.3. so there might be conflicts sometimes


# 2. Define the exact path to your custom R library
export R_LIBS_USER="/home/tmhaxs421/R/x86_64-redhat-linux-gnu-library/4.5"

# 3. Run your R script
Rscript my_script.R