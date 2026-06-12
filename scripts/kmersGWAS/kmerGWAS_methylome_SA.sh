KMC='~/kmersGWAS/external_programs/kmc_v3'
BIN='~/kmersGWAS/bin/'
DIR='~/output_kmersGWAS/'
FASTQ='~/path_to_fastq.txt'

index=$LSB_JOBINDEX
sample=$(sed -n ${index}p samples.txt)

input_file="${DIR}samples/$sample/input_files_$sample.txt"

mkdir ${DIR}samples/$sample

grep -F ${sample}. $FASTQ > $input_file

$KMC -t2 -k31 -ci2 @$input_file ${DIR}samples/$sample/kmc_canon \
  ${DIR}samples/$sample > ${DIR}samples/$sample/kmc_canon.log.out

$KMC -t2 -k31 -ci0 -b @$input_file \
  ${DIR}samples/$sample/kmc_non_canon \
    ${DIR}samples/$sample > ${DIR}samples/$sample/kmc_non_canon.log.out

${BIN}kmers_add_strand_information -c ${DIR}samples/$sample/kmc_canon \
  -n ${DIR}samples/$sample/kmc_non_canon -k 31 \
  -o ${DIR}samples/$sample/kmers_with_strand > ${DIR}samples/$sample/add_strand.log.out

rm ${DIR}samples/$sample/*.kmc*
