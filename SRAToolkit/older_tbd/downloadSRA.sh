# TODO: This following script will work well for single end reads but will not work for paired end reads. So add the functionality for downloading paired end reads.
# Download the SRA accession numbers to download the required files


module load sratoolkit/3.2.0


OUT_DIR="/home/tmhaxs421/brannanlab/VS_share/VS_AS_folder/Clip_analysis/skipper/data/RBFOX2_Nature_methods"  # change
ACC_LIST="/home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/SRAToolkit/eg_SraAccList.tsv"        # change (TSV!)

SRA_DIR="${OUT_DIR}/sra"
FASTQ_DIR="${OUT_DIR}/fastq"
TMP_DIR="${OUT_DIR}/tmp"
THREADS=32

mkdir -p "${SRA_DIR}" "${FASTQ_DIR}" "${TMP_DIR}"

# 1) prefetch needs only accessions (col1) - The file will be deleted after prefetch ---------------------------------------------------------------------
ACC_ONLY="${TMP_DIR}/acc_only.txt"
cut -f1 "${ACC_LIST}" | sed '/^$/d' > "${ACC_ONLY}"

cd "${SRA_DIR}"
prefetch --option-file "${ACC_ONLY}"

rm -f "${ACC_ONLY}"


while IFS=$'\t' read -r acc sample; do

  fasterq-dump "${acc}" \
    --outdir "${FASTQ_DIR}" \
    --temp "${TMP_DIR}" \
    -e "${THREADS}" \
    -f \
    --split-3

  fq="${FASTQ_DIR}/${acc}.fastq"

  mv "${fq}" "${FASTQ_DIR}/${sample}.fastq"

  # compress renamed fastqs
  pigz -p "${THREADS}" "${FASTQ_DIR}/${sample}"*.fastq

done < "${ACC_LIST}"