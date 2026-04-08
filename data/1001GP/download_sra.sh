i=$1

if [[ ! -e ${i}.sra ]]; then
  first_6_chars=$(echo $i | cut -c1-6)
  accession="${i%.*}"
  
  # Download SRA file
  echo "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${first_6_chars}/${accession}/${accession}.sra"
  wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${first_6_chars}/${accession}/${accession}.sra
  
  # Extract fastq file from SRA
  echo "fastq-dump --split-3 ${i}.sra"
  fastq-dump --split-3 ${i}.sra
  
  # Compress fastq file(s) (1 or 2 files for SE or PE libraries, respectively)
  gzip ${i}*.fastq
  
  # Remove SRA file
  rm ${i}.sra
else
  echo "${i}.sra already exists"
fi
