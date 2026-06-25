#!/bin/bash

#SBATCH --job-name=check_file_size
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_%u_%x_%j.log


# du -sh  * 2>/home/tmhaxs421/brannanlab/tmhaxs421/du_noaccess.log

path="/home/tmhaxs421/brannanlab/*"

ls -d $path | xargs -P 6 -I {} du -sh "{}" 2>/home/tmhaxs421/brannanlab/tmhaxs421/du_noaccess.log