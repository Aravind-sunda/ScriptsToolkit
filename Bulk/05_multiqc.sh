#!/bin/bash
#SBATCH --partition=defq
#SBATCH --job-name=multiqc_bulk
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --time=96:00:00
#SBATCH --output=slurm_%u_%x_%j.log

mamba deactivate
mamba activate bioinformatics # running multi qc reports for all the stats and strandedness files
multiqc --force "$WORKING_DIR" -o "$WORKING_DIR/multiqc_reports" --filename "multiqc_report.html" --ignore ".*/"