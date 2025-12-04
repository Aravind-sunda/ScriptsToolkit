module load mamba

# Making the dbSNP for mm10
# for b in $(ls *.bed); do echo $b; tail -n+2 $b | cut -f1,2,3  | sort >> ../mm10_dbsnp_combined.bed3; done


# Explanation of the arguments:
# 1. samples_path
# 2. sample_name-> for the json file output
# 3. json_path
# 4. sailor output_path
# 5. genome_fasta
# 6. dbsnp_bed
# 7. edit to calculate
# 8. strandedness

python make_sailor_json.py \
/home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/subset-bam \
t1_shCTRL \
/home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/json \
/home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/sailor_op \
/home/tmhaxs421/brannanlab/10x_genomics/Mouse_genome_10x/mRuby3_mm10/mRUBY3_mm10/fasta/genome.fa \
/home/tmhaxs421/brannanlab/tmhaxs421/TCA/Short_Read_SAILOR/SAILOR/dbsnp142-mm10/mm10_dbsnp_combined.bed3 \
GA \
true