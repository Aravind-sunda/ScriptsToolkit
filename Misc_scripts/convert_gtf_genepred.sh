# conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed
GTF=""
OUTPUT=""       # base name for output files (e.g. myAnnotation)
GENOMEFASTA=""

GTFTOGENEPRED=/home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/NCBI_utilities/gtfToGenePred
IGVTOOLS=/home/tmhaxs421/brannanlab/tmhaxs421/applications/IGV_2.19.8/igvtools

# Convert GTF → genePred
"$GTFTOGENEPRED" -genePredExt -geneNameAsName2 "$GTF" "${OUTPUT}.genePred"

# Prepend bin column
awk 'BEGIN{OFS="\t"} {print "0", $0}' "${OUTPUT}.genePred" > "${OUTPUT}.refGene.txt"

# Sort → produces the IGV-ready file
"$IGVTOOLS" sort "${OUTPUT}.refGene.txt" "${OUTPUT}.sorted.refGene.txt"

# Remove intermediate files
rm -f "${OUTPUT}.genePred" "${OUTPUT}.refGene.txt"

# Index the sorted file → produces the .idx IGV needs
"$IGVTOOLS" index "${OUTPUT}.sorted.refGene.txt"

# Index the genome FASTA for IGV
module load samtools
samtools faidx "$GENOMEFASTA"