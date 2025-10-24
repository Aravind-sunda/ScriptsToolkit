mamba deactivate
mamba activate bioinformatics # running multi qc reports for all the stats and strandedness files
multiqc --force "$WORKING_DIR" -o "$WORKING_DIR/multiqc_reports" --filename "multiqc_report.html" --ignore ".*/"