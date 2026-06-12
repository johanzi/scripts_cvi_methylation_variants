High impact mutations drive DNA methylation variation after colonization
of a novel habitat
================
Johan Zicola
2026-06-12 12:13:22

- [Overview](#overview)
- [Softwares required](#softwares-required)
- [WGBS library preparation](#wgbs-library-preparation)
- [Sequencing](#sequencing)
- [Reanalysis of the 1001 GP data](#reanalysis-of-the-1001-gp-data)
  - [Get information samples](#get-information-samples)
  - [Download data](#download-data)
  - [Mapping with Bismark](#mapping-with-bismark)
    - [Get reference genome](#get-reference-genome)
    - [Mapping single-end data](#mapping-single-end-data)
    - [Mapping paired-end data](#mapping-paired-end-data)
- [Analysis of the accessions from Cape Verde and
  Morocco](#analysis-of-the-accessions-from-cape-verde-and-morocco)
  - [Download fastq files](#download-fastq-files)
  - [Run Bismark](#run-bismark)
  - [Methylation call with methylKit](#methylation-call-with-methylkit)
  - [Annotation Araport11](#annotation-araport11)
  - [R libraries and functions](#r-libraries-and-functions)
  - [Create methylKit objects](#create-methylkit-objects)
  - [Create methylRawListDB objects](#create-methylrawlistdb-objects)
  - [Load methylRawListDB objects](#load-methylrawlistdb-objects)
  - [Filter methylRawList raw](#filter-methylrawlist-raw)
  - [Load filtered methylRawListDB
    objects](#load-filtered-methylrawlistdb-objects)
  - [Create methlBaseDB objects](#create-methlbasedb-objects)
  - [Load methylBaseDB objects](#load-methylbasedb-objects)
  - [Subset genomic regions](#subset-genomic-regions)
  - [Load methylRaw subset data](#load-methylraw-subset-data)
  - [Methylation levels](#methylation-levels)
    - [Whole genome](#whole-genome)
    - [Genes](#genes)
    - [Plot gene body methylation all world (figure
      1b)](#plot-gene-body-methylation-all-world-figure-1b)
    - [Map of the world (Figure 1a)](#map-of-the-world-figure-1a)
    - [Map of Santo Antao (Figure 1a)](#map-of-santo-antao-figure-1a)
    - [All TEs](#all-tes)
    - [Long TEs](#long-tes)
  - [TE methylation compared to world-wide accessions (figure
    3)](#te-methylation-compared-to-world-wide-accessions-figure-3)
    - [mCG all TEs](#mcg-all-tes)
    - [mCHG all TEs](#mchg-all-tes)
    - [mCHH all TEs](#mchh-all-tes)
    - [mCG long TEs](#mcg-long-tes)
    - [mCHG long TEs](#mchg-long-tes)
    - [mCHH long TEs](#mchh-long-tes)
- [GWAS analysis](#gwas-analysis)
  - [Download raw fastq reads whole-genome
    sequencing](#download-raw-fastq-reads-whole-genome-sequencing)
  - [SNP calling](#snp-calling)
  - [SNP annotation in CVI](#snp-annotation-in-cvi)
  - [Prepare VCF file for the 83 CVI
    accession](#prepare-vcf-file-for-the-83-cvi-accession)
  - [Prepare phenotype](#prepare-phenotype)
    - [Generate phenotype files](#generate-phenotype-files)
  - [Run Gemma](#run-gemma)
    - [GWAS whole genome](#gwas-whole-genome)
    - [GWAS genes](#gwas-genes)
    - [GWAS all TEs](#gwas-all-tes)
    - [GWAS long TEs](#gwas-long-tes)
  - [Variants at SUVH4, AGO9, DRM1, and
    MET1](#variants-at-suvh4-ago9-drm1-and-met1)
    - [Subsample VCF file for the 190 Santo Antao
      accessions](#subsample-vcf-file-for-the-190-santo-antao-accessions)
    - [Run SnpEff](#run-snpeff)
    - [Summary variants](#summary-variants)
- [Allele status in CPV](#allele-status-in-cpv)
  - [Allele distribution by
    population](#allele-distribution-by-population)
  - [Map allele distribution by
    population](#map-allele-distribution-by-population)
  - [Coordinate file](#coordinate-file)
  - [VIM2 distribution](#vim2-distribution)
  - [CMT2 distribution](#cmt2-distribution)
  - [FBX5 distribution](#fbx5-distribution)
  - [Plot diagram VIM2](#plot-diagram-vim2)
  - [Plot diagram CMT2](#plot-diagram-cmt2)
  - [Plot diagram FBX5](#plot-diagram-fbx5)
  - [Plot diagram all alleles](#plot-diagram-all-alleles)
- [Plot gbM by VIM2/4 allele](#plot-gbm-by-vim24-allele)
- [Plot TE methylation by FBX5
  allele](#plot-te-methylation-by-fbx5-allele)
  - [mCG in long TEs](#mcg-in-long-tes)
  - [mCHG in long TEs](#mchg-in-long-tes)
    - [Statistics](#statistics)
  - [mCHH in long TEs](#mchh-in-long-tes)
    - [Statistics](#statistics-1)
    - [Test difference mCHH for FBX5 in CMT2ref/stop
      background](#test-difference-mchh-for-fbx5-in-cmt2refstop-background)
- [GWAS gbM with VIM2/4 insertion as
  covariate](#gwas-gbm-with-vim24-insertion-as-covariate)
- [GWAS mCG genome-wide with VIM2del as
  covariate](#gwas-mcg-genome-wide-with-vim2del-as-covariate)
- [GWAS mCHH long TEs with CMT2 as
  covariate](#gwas-mchh-long-tes-with-cmt2-as-covariate)
- [GWAS mCHH genome-wide CMT2 as
  covariate](#gwas-mchh-genome-wide-cmt2-as-covariate)
- [GWAS mCHG genome-wide CMT2 as
  covariate](#gwas-mchg-genome-wide-cmt2-as-covariate)
- [GWAS mCG long TEs with FBX5 as
  covariate](#gwas-mcg-long-tes-with-fbx5-as-covariate)
- [kmersGWAS](#kmersgwas)
  - [Running KMC](#running-kmc)
  - [Create k-mers list to be used in GWAS from all
    individuals](#create-k-mers-list-to-be-used-in-gwas-from-all-individuals)
  - [Filter k-mers from separate lists to one list with all k-mers to
    use](#filter-k-mers-from-separate-lists-to-one-list-with-all-k-mers-to-use)
  - [Run the GWAS](#run-the-gwas)
  - [Prepare output from GWAS for
    Bowtie2](#prepare-output-from-gwas-for-bowtie2)
  - [Map kmers to TAIR10](#map-kmers-to-tair10)
  - [Manhattan plots](#manhattan-plots)
    - [Figure 2b](#figure-2b)
    - [Figure S1a](#figure-s1a)
    - [Figure S1b](#figure-s1b)
  - [Kmer mapping in the four CVI
    assemblies](#kmer-mapping-in-the-four-cvi-assemblies)
    - [Identify coordinates VIM2/4 and VIM3 in the four
      assemblies](#identify-coordinates-vim24-and-vim3-in-the-four-assemblies)
    - [Blast kmer on the five
      assemblies](#blast-kmer-on-the-five-assemblies)
- [DMR analysis for VIM2](#dmr-analysis-for-vim2)
  - [Pooling of the data](#pooling-of-the-data)
  - [Run Bismark](#run-bismark-1)
  - [Create methylKit objects](#create-methylkit-objects-1)
    - [Create methylRawListDB
      objects](#create-methylrawlistdb-objects-1)
    - [Load methylRawListDB objects](#load-methylrawlistdb-objects-1)
    - [Filter methylRawList raw](#filter-methylrawlist-raw-1)
    - [Load filtered methylRawListDB
      objects](#load-filtered-methylrawlistdb-objects-1)
  - [Gene methylation](#gene-methylation)
  - [Create methlyBaseDB objects](#create-methlybasedb-objects)
  - [Load methylBaseDB objects](#load-methylbasedb-objects-1)
  - [Create tiles objects](#create-tiles-objects)
  - [Load tiles objects](#load-tiles-objects)
  - [Create methylDiffDB objects](#create-methyldiffdb-objects)
  - [Load methylDiffDB objects](#load-methyldiffdb-objects)
  - [Create DMRs](#create-dmrs)
  - [Distribution of mCG DMRs](#distribution-of-mcg-dmrs)
  - [Merge mCG DMRs](#merge-mcg-dmrs)
  - [Overlap with genes](#overlap-with-genes)
  - [GO analysis](#go-analysis)
  - [Overlap of mCG DMRs with promoter
    regions](#overlap-of-mcg-dmrs-with-promoter-regions)
- [Analysis DNA methylation in vim
  mutants](#analysis-dna-methylation-in-vim-mutants)
  - [Run Bismark](#run-bismark-2)
  - [Create methylKit objects](#create-methylkit-objects-2)
  - [Create methylRawListDB objects](#create-methylrawlistdb-objects-2)
  - [Load methylRawListDB objects](#load-methylrawlistdb-objects-2)
  - [Filter methylRawList raw](#filter-methylrawlist-raw-2)
  - [Load filtered methylRawListDB
    objects](#load-filtered-methylrawlistdb-objects-2)
  - [Subset genomic regions](#subset-genomic-regions-1)
  - [Subset data](#subset-data)
  - [Load methylRaw subset data](#load-methylraw-subset-data-1)
  - [Genome-wide methylation](#genome-wide-methylation)
  - [Genes](#genes-1)
  - [All TEs](#all-tes-1)
  - [Long TEs](#long-tes-1)
  - [Whole genome methylation compared to
    SA](#whole-genome-methylation-compared-to-sa)
  - [Gene body methylation compared to
    SA](#gene-body-methylation-compared-to-sa)
  - [TEs](#tes)
- [Analysis DNA methylation in SAIL and SALK fbx5 and cmt2
  mutants](#analysis-dna-methylation-in-sail-and-salk-fbx5-and-cmt2-mutants)
  - [Samples](#samples)
  - [Run Bismark](#run-bismark-3)
  - [Create methylKit objects](#create-methylkit-objects-3)
  - [Create methylRawListDB objects](#create-methylrawlistdb-objects-3)
  - [Load methylRawListDB objects](#load-methylrawlistdb-objects-3)
  - [Filter methylRawList raw](#filter-methylrawlist-raw-3)
  - [Load filtered methylRawListDB
    objects](#load-filtered-methylrawlistdb-objects-3)
  - [Subset genomic regions](#subset-genomic-regions-2)
  - [Subset data](#subset-data-1)
  - [Load methylRaw subset data](#load-methylraw-subset-data-2)
  - [mCG in long TEs for FBX5
    mutants](#mcg-in-long-tes-for-fbx5-mutants)
  - [mCHG in long TEs for FBX5
    mutants](#mchg-in-long-tes-for-fbx5-mutants)
  - [mCHH in long TEs for FBX5
    mutants](#mchh-in-long-tes-for-fbx5-mutants)
  - [mCHH in whole genome for FBX5
    mutant](#mchh-in-whole-genome-for-fbx5-mutant)
  - [Statistical tests](#statistical-tests)
    - [Col-0 background](#col-0-background)
    - [Col-3 background](#col-3-background)
  - [mCHH in long TEs for cmt2-5
    mutants](#mchh-in-long-tes-for-cmt2-5-mutants)
    - [Statistical test](#statistical-test)
  - [Compare cmt2 and SA with
    CMT2stop](#compare-cmt2-and-sa-with-cmt2stop)
    - [Statistical difference](#statistical-difference)
- [Analysis DNA methylation fbx5 CRISPR
  lines](#analysis-dna-methylation-fbx5-crispr-lines)
  - [Samples](#samples-1)
  - [Create methylRawListDB objects](#create-methylrawlistdb-objects-4)
  - [Load methylRawListDB objects](#load-methylrawlistdb-objects-4)
  - [Filter methylRawList raw](#filter-methylrawlist-raw-4)
  - [Load filtered methylRawListDB
    objects](#load-filtered-methylrawlistdb-objects-4)
  - [Subset genomic regions](#subset-genomic-regions-3)
  - [Subset data](#subset-data-2)
  - [Load methylRaw objects per
    regions](#load-methylraw-objects-per-regions)
  - [Create methlBaseDB objects](#create-methlbasedb-objects-1)
  - [Load methylBaseDB objects](#load-methylbasedb-objects-2)
  - [Long TE methylation](#long-te-methylation)
    - [mCG](#mcg-6)
    - [mCHG](#mchg-6)
    - [mCHH](#mchh-6)
- [RNA-seq library preparation](#rna-seq-library-preparation)
  - [Read trimming](#read-trimming)
  - [Mapping](#mapping)
    - [Get reference genome](#get-reference-genome-1)
    - [Get gene annotation Araport11](#get-gene-annotation-araport11)
    - [Mapping](#mapping-1)
  - [Read counting](#read-counting)
    - [Install HTseq](#install-htseq)
    - [Gene annotation for htseq](#gene-annotation-for-htseq)
    - [Keep only CDS and exons](#keep-only-cds-and-exons)
    - [Remove miRNAs genes](#remove-mirnas-genes)
    - [Counting](#counting)
    - [Count merging](#count-merging)
  - [Analysis in R](#analysis-in-r)
    - [R Libraries](#r-libraries)
    - [Load cts and coldata](#load-cts-and-coldata)
    - [PCA](#pca)
  - [Replicates analysis (20
    accessions)](#replicates-analysis-20-accessions)
    - [PCA analysis](#pca-analysis)
  - [Analysis without replicates (97
    accessions)](#analysis-without-replicates-97-accessions)
    - [DEG analysis by VIM2 allele](#deg-analysis-by-vim2-allele)
- [Analysis TE expression](#analysis-te-expression)
  - [R Libraries](#r-libraries-1)
  - [TE annotation](#te-annotation)
  - [Reference genome](#reference-genome)
  - [Bowtie2 index for each TE](#bowtie2-index-for-each-te)
  - [STEP I: Mapping RNA-seq data as first step of
    RepEnrich2](#step-i-mapping-rna-seq-data-as-first-step-of-repenrich2)
  - [STEP II: Split uniquely mapped and multimapping
    reads](#step-ii-split-uniquely-mapped-and-multimapping-reads)
  - [STEP III: Mapping to TEs](#step-iii-mapping-to-tes)
  - [Summarize data in read count
    matrix](#summarize-data-in-read-count-matrix)
  - [Analysis in R](#analysis-in-r-1)
    - [Replicates analysis](#replicates-analysis)
    - [Analysis on single accessions and long
      TEs](#analysis-on-single-accessions-and-long-tes)
    - [Analysis by CMT2 allele](#analysis-by-cmt2-allele)
    - [Analysis by FBX5 allele](#analysis-by-fbx5-allele)
    - [Permutation](#permutation)
- [TE transposition using read coverage as
  proxy](#te-transposition-using-read-coverage-as-proxy)
  - [CMT2](#cmt2)
  - [FBX5](#fbx5)
- [Marginal genealogical tree with
  RELATE](#marginal-genealogical-tree-with-relate)
  - [VIM2](#vim2)
  - [CMT2](#cmt2-1)
  - [FBX5](#fbx5-1)
- [Inference of selection
  coefficient](#inference-of-selection-coefficient)
  - [VIM2](#vim2-1)
  - [CMT2](#cmt2-2)
  - [FBX5](#fbx5-2)
- [Selective sweep analysis (figure
  7)](#selective-sweep-analysis-figure-7)
- [Authors](#authors)
- [License](#license)

# Overview

This documentation explains step by step how was performed the analysis
on whole-genome bisulfite sequencing (WGBS) data on African *Arabidopsis
thaliana* accessions (Morocco and Cape Verde) and the reanalysis of a
subset of the WGBS data from the 1001 Genome Project (1001GP)
([Kawakatsu et al.,
2016](http://www.sciencedirect.com/science/article/pii/S0092867416308522))
and the RNA-seq performed on Cape Verde accessions.

# Softwares required

- Bismark (v0.19.0)
- Python3.5
- GEMMA (v0.94)
- vcftools (v0.1.14)
- bcftools (v1.2)
- bwa (v0.7.15)
- R (\>3.3.0)
- SRA tool kit (v3.2.1)
- FastQC (v0.11.9)
- multiqc (v1.6)
- cutadapt (v2.9
- methylKit (v1.14.2)
- HISAT2 (v11.4.0)

# WGBS library preparation

Libraries were prepared as described previously in [Urich et
al. 2012](http://www.nature.com/nprot/journal/v10/n3/full/nprot.2014.114.html)
with minor modifications.

# Sequencing

Libraries were pooled based on the 24 NEXTFlex Bisulfite-Seq barcodes
(BiooScientific) for multiplex sequencing on the HiSeq3000 sequencer
(Illumina) in 150 bp single-end mode. 1 Gb (7 M reads) of data were
ordered for each library (minimum required by the sequencing facility).
Reads were trimmed from adapters by the sequencing facility using
Cutadapt ([Martin et al.,
2011](http://journal.embnet.org/index.php/embnetjournal/article/view/200)).
Visual inspection on graphics produced by
[fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was
used to visually determine the quality of the reads.

# Reanalysis of the 1001 GP data

## Get information samples

Working directory:
`/srv/netscratch/dep_coupland/grp_hancock/johan/bs-seq_data_1001`

Go to NCBI SRA selector for project PRJNA187927 (SALK La Jolla) and
PRJNA236110 (GMI Vienna) Web:
<https://www.ncbi.nlm.nih.gov/Traces/study/?go=home> enter PRJNA187927,
Download the SraRunTable.txt (478 samples)

`mv SraRunTable.txt SraRunTable_SALK.txt`

Web: <https://www.ncbi.nlm.nih.gov/Traces/study/?go=home> enter
PRJNA236110, Download the SraRunTable.txt (2215 samples)

`mv SraRunTable.txt SraRunTable_GMI.txt`

Go to <http://1001genomes.org/accessions.html> and download the csv file
(link at the bottom of the page)

Convert CSV to tab-separated file

``` bash
sed 's/,/\t/g' query.txt > accessions_1001genome.txt
```

A subset of 526 accessions were selected, discarding accessions from USA
and replicates. A dataframe containing details about each accessions was
created and put in `accessions_1001GP_figure1.txt`.

## Download data

Select SRR names of each the 526 accessions to download

``` bash
cut -f2 accessions_1001GP_figure1.txt > to_download.txt
```

The script `download_sra.sh` retrieves fastq file for each SRR number.

``` bash
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
```

Download the data

``` bash

# Launch script in bsub
while read i; do
  bash download_sra.sh $i
done < to_download.txt
```

## Mapping with Bismark

Reads were mapped on *A. thaliana* TAIR10 reference [fasta
file](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas).
Note that the fasta file contains the 5 chromosomes and the 2 plastids
(chloroplast and mitochondria). In order to be processed by the R
package methylKit ([Akalin et al.,
2011](https://doi.org/10.1186/gb-2012-13-10-r87)), the cytosine report
files from Bismark were generated for each chromosome and each
methylation context. In order to perform the analysis on many
accessions, the bash script [run_bismark.sh](scripts/run_bismark.sh)
performs Bismark analysis step-by-step.

Note that absolute paths can be given, the output files will be
generated in the specified output directory (argument `-o`) or by
default in the directory containing the input fastq file if no -o
argument is specified. While running, the script will echo each step
performed, which can be redirected to a log file.

The script will take care of: \* Building the bismark reference genome
\* Perform the alignment \* Remove duplicate reads from the bam file \*
Extract the methylation status and generate coverage and bedGraph files
(visualization in SeqMonk and IGV) \* Generate cytosine reports (used as
input file for methylKit R package) \* Calculate conversion efficiency
(based on spurious non-converted cytosines from the chloroplast genome)

To get more information on how run_bismark.sh is working:

``` bash
bash run_bismark.sh -h
```

The code itself contains comments for each step so have look at it and
tweak it to your needs.

### Get reference genome

Download the fasta file for *A. thaliana* TAIR10 reference:

``` bash
# Download fasta file
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

# Rename chromosomes to add prefix "Chr" (easier to retrieve information, e.g. Chr1 rather than 1)
sed -i 's/^>\([1-5]\)/>Chr\1/g' TAIR10_chr_all.fas 

# Move the file in "/path/to/dir_fasta/"
mv TAIR10_chr_all.fas /path/to/dir_fasta/
```

### Mapping single-end data

The data from the 1001 GP have a mix of SE and PE data, put fastq files
in separate folders.

``` bash
# Split SE and PE data in different folders
mkdir PE_data
mkdir SE_data

# Move paired-end data in PE_data
mv *_1.fastq.gz PE_data/
mv *_2.fastq.gz PE_data/

# Move the rest (SE data) into SE_data
mv *fastq.gz SE_data
```

Map the reads with the light mode on (remove unnecessary intermediary
files) with the flag `-l`

``` bash
while read i; do
bash run_bismark.sh -l -r </path/to/dir_fasta/> -1 ${i} -o </path/to/output/> 
done < <(ls *fastq.gz)
```

#### Assess mapping efficiency

``` bash

for i in *bismark_bt2_SE_report.txt; do
  library=$(echo $i | cut -d'_' -f1,2)
  map=$(grep "Mapping efficiency" $i | cut -d':' -f2 -)
  echo -e "${library}\t${map}" >> mapping_efficiency.txt
done
```

Put data in excel and calculate average and SD average: 64.89% Stdev:
9.34%

#### Assess conversion efficiency

``` bash
for i in *_report.conversion_efficiency.txt; do 
  name=$(echo $i | cut -d'_' -f1,2) 
  line=$(tail -n1 $i | cut -f3)
  echo $name $line
done >> conversion_efficiency.txt
```

Put data in excel and calculate average and SD

average: 99.52% Stdev: 0.38%

### Mapping paired-end data

``` bash
# Get single name for each pair data
ls *fastq.gz | cut -d'_' -f1,2 | uniq > list_fastq_files.txt

while read i; do
          fastq1=${i}_1.fastq.gz
          fastq2=${i}_2.fastq.gz
            bash run_bismark.sh -l -r </path/to/dir_fasta/> -1 $fastq1 -2 $fastq2 -o </path/to/output/> 
done < list_fastq_files.txt
```

#### Assess mapping efficiency

``` bash
cd /srv/netscratch/dep_coupland/grp_hancock/johan/bs-seq_data_1001/fastq_files/1001/PE_data

for i in *bismark_bt2_PE_report.txt.gz; do
  library=$(echo $i | cut -d'_' -f1,2)
  map=$(zgrep "Mapping efficiency" $i | cut -d':' -f2 -)
  echo -e "${library}\t${map}" >> mapping_efficiency.txt
done
```

average: 51.82% Stdev: 8.21%

Interestingly, the SE end data map at higher efficiency than PE data
(about 10% more uniquely mapped reads) and I observed the same trend for
the data of project GC_4050. This is probably due to the fact that more
reads are unlinked in SE mode and therefore reads mapping in repetitive
regions in PE mode are 2 times more numerous as they belong to the same
DNA fragment.

#### Assess conversion efficiency

``` bash
for i in *_report.conversion_efficiency.txt; do 
  name=$(echo $i | cut -d'_' -f1,2) 
  line=$(tail -n1 $i | cut -f3)
  echo $name $line
done >> conversion_efficiency.txt
```

# Analysis of the accessions from Cape Verde and Morocco

We generated WGBS data for 83 accessions from Cape Verde - Santo Antao,
20 from Morocco and diverse mutants (for FBX5 and CMT2). We also
included as control the accessions Col-0, Col-3, Doer-10, and UKID116.
The data for the fastq files for these samples can be downloaded in the
NCBI depository PRJNA612437.

## Download fastq files

1.  Go on website <https://www.ncbi.nlm.nih.gov/sra> and type
    PRJNA612437.
2.  Click to “SRA Experiments”
3.  Click on “Send Results to Run selector”
4.  Download the SRR list by clicking ‘Accession List’ as
    SRR_Acc_List.txt

*NB: The data are all single-end reads*

In bash, download SRA files and convert them in fastq files:

``` bash
while read name in list; do
    fastq-dump --split-spot $name
done < SRR_Acc_List.txt
```

## Run Bismark

For each fastq file, run the following command:

``` bash
bash run_bismark.sh -1 <filename.fastq> -r </path/to/dir_fasta/> -o </name/output/directory/>
```

## Methylation call with methylKit

Once Bismark has been run, cytosine report files for each methylation
contexts are imported in methylKit R packages for further analysis.

Typical format of the bismark output file imported in methylKit:

The suffix of the file is
`*_bismark_bt2.deduplicated.bismark.cov.gz.CHG_report_only_chr.txt` so
it contains all the call for cytosines in CHG context, excluding the
calls in organelles (chloroplast and mitochondria).

    Chr4    1004    +       1       1       CHG     CAG
    Chr4    1006    -       2       3       CHG     CTG
    Chr4    1009    -       0       5       CHG     CCG
    Chr4    1020    +       1       3       CHG     CCG
    Chr4    1023    -       4       5       CHG     CCG
    Chr4    1052    -       9       3       CHG     CCG
    Chr4    1071    +       3       4       CHG     CAG
    Chr4    1073    -       6       6       CHG     CTG
    Chr4    1090    +       4       5       CHG     CCG
    Chr4    1096    +       4       5       CHG     CTG

The columns represent chromosome, base position, strand position of the
cytosine, number of methylated Cs, number unmethylated Cs, methylation
context, and trinucleotide context.

Different functions were created to run the functions of methylKit in
batch. These functions can be found in
[scripts/functions_methylkit.R](functions_methylkit.R).

Here is described the pipeline used to process the methylation data in
methylKit.

## Annotation Araport11

We need the last annotation of genes and TEs for Col-0. Use Araport11
annotation.

Files downloaded from
<https://www.arabidopsis.org/download/list?dir=Genes%2FAraport11_genome_release>

``` bash
# Get genes
grep -P "\tgene\t" Araport11_GFF3_genes_transposons.201606.gff  > Araport11_GFF3_genes_only.gff

# Convert gff to bed format
gff2bed < Araport11_GFF3_genes_only.gff > Araport11_GFF3_genes_only_full.bed

# Keep only first 4 columns
cut -f1,2,3,4 Araport11_GFF3_genes_only_full.bed > Araport11_GFF3_genes_only.bed

# For TEs
grep "transposable_element" Araport11_GFF3_genes_transposons.201606.gff | wc -l
35090

# In comparison, there were 35082 transposable_element feature in TAIR10 annotation

grep "transposable_element" Araport11_GFF3_genes_transposons.201606.gff  > Araport11_GFF3_transposons.gff

# Convert gff to bed format
gff2bed < Araport11_GFF3_transposons.gff > Araport11_GFF3_transposons_full.bed

# Keep only first 4 columns
cut -f1,2,3,4 Araport11_GFF3_transposons_full.bed > Araport11_GFF3_transposons.bed

wc -l Araport11_GFF3_transposons.bed
31189 Araport11_GFF3_transposons.bed

# Keep TEs bigger than 4 kb
cat Araport11_GFF3_transposons.bed | awk -F'\t' '$3-$2 >= 4000 {print $0}' > Araport11_GFF3_transposons_longer_4kb.bed

wc -l Araport11_GFF3_transposons_longer_4kb.bed
1235 Araport11_GFF3_transposons_longer_4kb.bed
```

We have a total of 31,189 TEs, including 1,235 TEs bigger than 4 kb.

## R libraries and functions

``` r
####################################
#       Libraries and functions    #
####################################

library(scatterpie)
library(plyr)
library(ggmap)

# Load the R script functions_methylkit.R which contains wrap up functions to run in batch several
# methylKit functions
source("scripts/functions_methylkit.R")

####################################
#      Paths to DB directories     #
####################################

# Paths to database for the output files of methylKit
path_DB_CpG <- "F:/NETSCRATCH/methylKit_DB_files/pooled_data/methylDB_CpG"
path_DB_CHG <- "F:/NETSCRATCH/methylKit_DB_files/pooled_data/methylDB_CHG"
path_DB_CHH <- "F:/NETSCRATCH/methylKit_DB_files/pooled_data/methylDB_CHH"

# Create a list of these 3 paths
list_DB_paths <- list(path_DB_CpG, path_DB_CHG, path_DB_CHH)

# Path containing cytosine report and bam files from bismark pipeline
path_bismark_files <- paste("/path/to/bismark/output/files/", sep = "")

####################################################
################# BED FILES ########################
####################################################

path_bed <- "data/bed_files/"

# Path to bed files for region analysis
bed_genes <- paste(path_bed, "Araport11_GFF3_genes_only.bed", sep = "")

# bed_genes_annotate <- paste(workdir, "Arabidopsis_thaliana.TAIR10.39.bed", sep="")
# Version that was made from GTF (works with readTranscriptFeatures)

# All TEs
bed_TEs <- paste(path_bed, "Araport11_GFF3_transposons.bed", sep = "")

# TEs longer than 4 kb
bed_TEs_4kb <- paste(path_bed, "Araport11_GFF3_transposons_longer_4kb.bed", sep = "")

# TEs shorter than 500 bp
bed_TEs_500bp <- paste(path_bed, "Araport11_GFF3_transposons_smaller_500bp.bed", sep = "")


####################################################
################# ACCESSIONS FILES ########################
####################################################

# Path to file with accession information (several were 
# used in the different analyses and they are all available in GitHub)
path_df_accessions <- "data/df_accessions_83.txt"

# Get information of the accessions and generate a table
# Order first the file to export so that the elements are ordered as the fastq files (3542_AA, 3542_AB, ...)
df_accessions <- read.table(path_df_accessions, header = TRUE, stringsAsFactors = TRUE)

# Order the accession as list.files() list the bismark cytosine report files
df_accessions <- order_df_accessions(df_accessions)

# I need to create an hybrid name otherwise the loading of the 
# file won't respect the original order of the input bismark file
df_accessions$sample <- paste(df_accessions$library, df_accessions$name, sep = "_")

# Make a list of samples
list_samples <- as.list(as.vector(df_accessions$sample))

# Get list of treatments and reformat so that the first 
# treatment is 0 (control should be 0 optimally)
# Here I put as example CMT2 allele but the variable used 
# as treatment differ in different analysis
list_treatments <- as.numeric(df_accessions$CMT2)

# Vector of the 3 contexts analyzed
context <- c("CpG", "CHG", "CHH")
```

## Create methylKit objects

## Create methylRawListDB objects

The function `import_bismark_cytosine_report` will retrieve
automatically the different cytosine report files generated by Bismark
and will create flat database files, allowing to reduce RAM usage.

``` r
import_bismark_cytosine_report(path_bismark_files, list_DB_paths, list_samples, list_treatments)
```

## Load methylRawListDB objects

Once created, load methylRawListDB objects. The files won’t actually be
loaded but accessed in real time when needed.

``` r
list_methylRawLists <- load_methylRawListDB(list_DB_paths, type = "raw", list_samples, list_treatments)
```

## Filter methylRawList raw

Keep only cytosine positions that have a define minimum coverage. This
threshold is usually set at around 5 in most WGBS analyses but since our
samples were sequenced at the minimum depth allowed by the sequencing
facility, we defined a lower threshold (minimum 2). This approach is
valid considering that we look at pattern across large genomic regions.
We assumed we would catch any strong signal if any.

``` r
filter_methylRawList(list_methylRawLists_raw)
```

## Load filtered methylRawListDB objects

``` r
list_methylRawLists <- load_methylRawListDB(list_DB_paths, type = "filtered", list_samples, list_treatments)
```

## Create methlBaseDB objects

Create methylBase objects based on given list_methylRawLists object.
Create DB files if not existing.

``` r
list_methylBases <- merged_methylRawList(list_methylRawLists, suffix="filtered_83")
```

## Load methylBaseDB objects

``` r
# Make lists of objects
list_methylBases <- load_methylBaseDB(list_DB_paths, list_samples, list_treatments, suffix="filtered_83")
```

## Subset genomic regions

We want now to analyze methylation patterns is specific genomic regions.
For this, we need to subset our data and generate new DB flat files for
the different regions.

``` r
# Create subset for methylRawList
subset_methylObject(list_methylRawLists, list_DB_paths, bed_genes, "genes", "methylRaw")

subset_methylObject(list_methylRawLists, list_DB_paths, bed_TEs, "TEs", "methylRaw")

subset_methylObject(list_methylRawLists, list_DB_paths, bed_TEs_4kb, "TEs_4kb", "methylRaw")
```

## Load methylRaw subset data

``` r
# Load subset data

# Load methylRawListDB objects (without filtering)
list_methylRawLists_genes <- load_methylRawListDB(list_DB_paths, type = "genes", list_samples, list_treatments)

list_methylRawLists_TEs <- load_methylRawListDB(list_DB_paths, type = "TEs", list_samples, list_treatments)

list_methylRawLists_TEs_4kb <- load_methylRawListDB(list_DB_paths, type = "TEs_4kb", list_samples, list_treatments)
```

## Methylation levels

We want first to visualize the methylation levels in different genomic
regions. For this, we extract the weighted methylation levels

### Whole genome

``` r
df_name <- "df_mean_filtered"
title <- "Weighted Methylation Level for genes"

get_df_wml(list_methylRawLists, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
ggplot_all(get(df_name), title = title)
```

### Genes

``` r
df_name <- "df_mean_genes"
title <- "Weighted Methylation Level for genes"

get_df_wml(list_methylRawLists_genes, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
ggplot_all(df_mean_genes, title = title)
```

### Plot gene body methylation all world (figure 1b)

``` r
require(gghighlight)

df_accessions <- read.table("data/methylKit_files/1001GP/df_accessions_1001_CPV_MOR.txt",
                            header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings="")  

path_DB <- "data/methylKit_files/1001GP"

df_name <- "df_mean_genes"
title <- "Weighted Methylation Level for genes"

get_df_wml(list_methylRawLists_genes, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
df <- merge(df_mean_genes, df_accessions, by="sample")

# Use dplyr prefix as desc function is also in IRange, which creates a conflict
df_subset <- df %>% group_by(country_code) %>% 
  filter(context=="CpG", !(country_code %in% c("USA","CPV-FO")), name!="SRR771702") %>% filter(n() > 7) %>% 
  arrange(-dplyr:::desc(Latitude))

# Change CPV-SA to CPV
levels(df_subset$country_code)[levels(df_subset$country_code)=="CPV-SA"] <- "CPV"

order_country_code <- unique(as.vector(df_subset$country_code))

df_subset$country_code <- factor(df_subset$country_code, levels=order_country_code, ordered=TRUE)

median_meth <- median(df_subset$percent_methylation)

# Check row number for Cvi-0 from 1001
which(grepl(6911, df_subset$seq_ID))
1

# Check our Cvi-0
which(grepl("4073_M", df_subset$seq_ID))
# 84

# The two Moroccan from the 1001
which(grepl(9606, df_subset$seq_ID))
# 90
which(grepl(9939, df_subset$seq_ID))
# 91

# UKID 116 (5822)
which(grepl(5822, df_subset$seq_ID))
# 549

# Dor-10 (5856)
which(grepl(5856, df_subset$seq_ID))
# 625

# I get overlying grey and black dots for CPV and FO
mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggplot(data=df_subset, aes(x=country_code, y=percent_methylation)) + geom_boxplot(outlier.shape = NA) + ggtitle("Gene body methylation") + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + ylab("% of methylated cytosines (CG)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Populations sorted by increasing latitude") +
  scale_y_continuous(limits=c(0,25), breaks=seq(0,25,5)) +
  gghighlight(country_code %in% c("CPV","MAR"), use_direct_label = FALSE) +
  geom_hline(yintercept=median(df_subset$percent_methylation), linetype="dashed", color = "black", size=0.5) + geom_jitter(data=df_subset[107:629,], colour="#d2d2d2ff", size=0.5, width=0.1, height=0) + geom_jitter(data=df_subset[1:106,], colour="black", size=0.5, width=0.1, height=0) +
  geom_jitter(data=df_subset[c(1,84,90,91),], colour=mycol, size=0.5, width=0.1, height=0) +
  geom_jitter(data=df_subset[1,],size=2, fill="#d2d2d2ff", shape=21, width=0.1, height=0) +
  geom_jitter(data=df_subset[84,],size=2, fill="red", shape=23, width=0.1, height=0) +
  geom_jitter(data=df_subset[90,],size=2, fill="#d2d2d2ff", shape=21, width=0.1, height=0) +
  geom_jitter(data=df_subset[91,],size=2, fill="#d2d2d2ff", shape=21, width=0.1, height=0) +
  geom_jitter(data=df_subset[549,],size=0.5, fill="black", shape=21, width=0.1, height=0)+
  geom_jitter(data=df_subset[625,],size=0.5, fill="black", shape=21, width=0.1, height=0)
```

![](images/figure1.png)

### Map of the world (Figure 1a)

``` r
df_accessions <- read.table("data/accessions_analysis_1001_Africa.txt",
                            header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings="")  

df_accessions$Latitude <- as.double(levels(df_accessions$Latitude)[df_accessions$Latitude])
df_accessions$Longitude <- as.double(levels(df_accessions$Longitude)[df_accessions$Longitude])

df_accessions_sub <- df_accessions %>% filter(Latitude > 15 & Latitude < 65) %>% filter(Longitude > -45 & Longitude < 80)

world <- map_data("world")

EU <- world[world$long > -45 & world$long < 80 & world$lat > 15 & world$lat < 65, ]

p <- ggplot(EU, aes(long, lat)) +
  geom_map(map = world, aes(map_id = region), fill = NA, color = "black") +
  coord_quickmap() + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p +  geom_point(data = df_accessions_sub, aes(x = Longitude, y = Latitude), color = "black", size = 1)
```

![](images/map_accessions_eurasia_bw.png)

``` r
df_accessions <- read.table("data/accessions_analysis_1001_Africa.txt",
                            header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings="")  

df_accessions$Latitude <- as.double(levels(df_accessions$Latitude)[df_accessions$Latitude])
df_accessions$Longitude <- as.double(levels(df_accessions$Longitude)[df_accessions$Longitude])

range <- c(left = -45, bottom=10, right=85, top=65)

# Create free account
register_stadiamaps("46ab633a-fe0d-4a09-b69e-7b6cfc4cce1f", write = FALSE)

#p <- get_stadiamap(range, zoom=5, maptype ="stamen_terrain") %>% ggmap()
p1 <- get_stadiamap(range, zoom=5, maptype ="stamen_terrain_background") %>% ggmap() 

p1 + geom_point(data = df_accessions, aes(x = Longitude, y = Latitude), color = "black", size = 1)
```

![](images/map_accessions_eurasia.png)

### Map of Santo Antao (Figure 1a)

``` r
world <- map_data("world")

SA <- world %>% filter(subregion =="Santo Antao") 

# Remove Cvi as GPS coordinates are not within Santo Antao
df_accessions_sub <-  df_accessions %>% filter(country_code=="CPV-SA") %>% filter(seq_ID!="6911")

p <- ggplot(SA, aes(long, lat)) +
  geom_map(map = world, aes(map_id = region), fill = NA, color = "black") +
  coord_quickmap() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p +  geom_point(data = df_accessions_sub, aes(x = Longitude, y = Latitude), color = "black", size = 1)
```

![](images/map_accessions_SA_bw.png)

``` r
df_accessions <- read.table("data/accessions_analysis_1001_Africa.txt",
  header = TRUE, sep = "\t", stringsAsFactors = TRUE, na.strings = ""
)

df_accessions$Latitude <- as.double(levels(df_accessions$Latitude)[df_accessions$Latitude])
df_accessions$Longitude <- as.double(levels(df_accessions$Longitude)[df_accessions$Longitude])

df_accessions_sub <- df_accessions %>%
  filter(country_code == "CPV-SA") %>%
  filter(seq_ID != "6911")

range <- c(left = -25.10, bottom = 17.09, right = -25.01, top = 17.135)

SA <- world %>% filter(subregion == "Santo Antao")

range <- c(left = min(SA$long) - 0.1, bottom = min(SA$lat) - 0.05, right = max(SA$long) + 0.1, top = max(SA$lat) + 0.05)

# Create free account on https://stadiamaps.com/ and get an API 
register_stadiamaps("46ab633a-fe0d-4a09-b69e-7b6cfc4cce1f", write = FALSE)

p1 <- get_stadiamap(range, zoom = 12, maptype = "stamen_terrain_background") %>% ggmap()

p1 + geom_point(data = df_accessions_sub, aes(x = Longitude, y = Latitude), color = "black", size = 1.5)
```

![](images/map_accessions_SA.png)

### All TEs

``` r
df_name <- "df_mean_TEs"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs, path_DB, df_name)

load_df_wml(path_DB, df_name)
```

### Long TEs

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)

load_df_wml(path_DB, df_name)
```

## TE methylation compared to world-wide accessions (figure 3)

``` r
require(gghighlight)

df_accessions <- read.table("data/df_accessions_1001_CPV_MOR.txt", 
                            header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings="")  

path_DB <- "F:/NETSCRATCH/methylKit_DB_files/1001_project"
```

### mCG all TEs

``` r
df_name <- "df_mean_TEs"

get_df_wml(list_methylRawLists_TEs, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
# ggplot_all(df_mean_genes, title = title)

df <- merge(df_mean_TEs, df_accessions, by = "sample")

# Use dplyr prefix as desc function is also in IRange, which creates a conflict
df_subset <- df %>%
  group_by(country_code) %>%
  filter(context == "CpG", country_code != "USA", name != "SRR771702") %>%
  filter(n() > 7) %>%
  arrange(-dplyr:::desc(Latitude))

df_subset <- df %>%
  group_by(country_code) %>%
  filter(context == "CpG", !(country_code %in% c("USA", "CPV-FO")), name != "SRR771702") %>%
  filter(n() > 7) %>%
  arrange(-dplyr:::desc(Latitude))

# Change CPV-SA to CPV
levels(df_subset$country_code)[levels(df_subset$country_code) == "CPV-SA"] <- "CPV"


order_country_code <- unique(as.vector(df_subset$country_code))

df_subset$country_code <- factor(df_subset$country_code, levels = order_country_code, ordered = TRUE)


median_meth <- median(df_subset$percent_methylation)

# Check row number for Cvi-0 from 1001
which(grepl(6911, df_subset$seq_ID))
1

# Check our Cvi-0
which(grepl("4073_M", df_subset$seq_ID))
# 84

# The two Moroccan from the 1001
which(grepl(9606, df_subset$seq_ID))
# 90
which(grepl(9939, df_subset$seq_ID))
# 91

# UKID 116 (5822)
which(grepl(5822, df_subset$seq_ID))
# 549

# Dor-10 (5856)
which(grepl(5856, df_subset$seq_ID))
# 625

# I get overlying grey and black dots for CPV and FO
mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggplot(data = df_subset, aes(x = country_code, y = percent_methylation)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("mCG all TEs") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("% of methylated cytosines (CG)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Populations sorted by increasing latitude") +
  scale_y_continuous(limits = c(40, 70), breaks = seq(40, 70, 5)) +
  gghighlight(country_code %in% c("CPV", "MAR"), use_direct_label = FALSE) +
  geom_hline(yintercept = median(df_subset$percent_methylation), linetype = "dashed", color = "black", size = 0.5) +
  geom_jitter(data = df_subset[107:629, ], colour = "#d2d2d2ff", size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[1:106, ], colour = "black", size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[c(1, 84, 90, 91), ], colour = mycol, size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[1, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[84, ], size = 2, fill = "red", shape = 23, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[90, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[91, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0)
```

![](images/mCG_all_TEs_worldwide.png)

### mCHG all TEs

``` r
df_name <- "df_mean_TEs"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
# ggplot_all(df_mean_genes, title = title)

library(gghighlight)

df <- merge(df_mean_TEs, df_accessions, by = "sample")

# Use dplyr prefix as desc function is also in IRange, which creates a conflict
df_subset <- df %>%
  group_by(country_code) %>%
  filter(context == "CHG", country_code != "USA", name != "SRR771702") %>%
  filter(n() > 7) %>%
  arrange(-dplyr:::desc(Latitude))

df_subset <- df %>%
  group_by(country_code) %>%
  filter(context == "CHG", !(country_code %in% c("USA", "CPV-FO")), name != "SRR771702") %>%
  filter(n() > 7) %>%
  arrange(-dplyr:::desc(Latitude))

# Change CPV-SA to CPV
levels(df_subset$country_code)[levels(df_subset$country_code) == "CPV-SA"] <- "CPV"

order_country_code <- unique(as.vector(df_subset$country_code))

df_subset$country_code <- factor(df_subset$country_code, levels = order_country_code, ordered = TRUE)

median_meth <- median(df_subset$percent_methylation)

# Check row number for Cvi-0 from 1001
which(grepl(6911, df_subset$seq_ID))
1

# Check our Cvi-0
which(grepl("4073_M", df_subset$seq_ID))
# 84

# The two Moroccan from the 1001
which(grepl(9606, df_subset$seq_ID))
# 90
which(grepl(9939, df_subset$seq_ID))
# 91

# UKID 116 (5822)
which(grepl(5822, df_subset$seq_ID))
# 549

# Dor-10 (5856)
which(grepl(5856, df_subset$seq_ID))
# 625

# I get overlying grey and black dots for CPV and FO

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggplot(data = df_subset, aes(x = country_code, y = percent_methylation)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("mCHG all TEs") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("% of methylated cytosines (CHG)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Populations sorted by increasing latitude") +
  scale_y_continuous(limits = c(10, 40), breaks = seq(10, 40, 5)) +
  gghighlight(country_code %in% c("CPV", "MAR"), use_direct_label = FALSE) +
  geom_hline(yintercept = median(df_subset$percent_methylation), linetype = "dashed", color = "black", size = 0.5) +
  geom_jitter(data = df_subset[107:629, ], colour = "#d2d2d2ff", size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[1:106, ], colour = "black", size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[c(1, 84, 90, 91), ], colour = mycol, size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[1, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[84, ], size = 2, fill = "red", shape = 23, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[90, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[91, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0)
```

![](images/mCHG_all_TEs_worldwide.png)

### mCHH all TEs

``` r
require(gghighlight)

df_name <- "df_mean_TEs"
title <- "Weighted Methylation Level for genes"

get_df_wml(list_methylRawLists_TEs, path_DB, df_name)

load_df_wml(path_DB, df_name)

df <- merge(df_mean_TEs, df_accessions, by = "sample")

# Use dplyr prefix as desc function is also in IRange, which creates a conflict
df_subset <- df %>%
  group_by(country_code) %>%
  filter(context == "CHH", country_code != "USA", name != "SRR771702") %>%
  filter(n() > 7) %>%
  arrange(-dplyr:::desc(Latitude))

df_subset <- df %>%
  group_by(country_code) %>%
  filter(context == "CHH", !(country_code %in% c("USA", "CPV-FO")), name != "SRR771702") %>%
  filter(n() > 7) %>%
  arrange(-dplyr:::desc(Latitude))

# Change CPV-SA to CPV
levels(df_subset$country_code)[levels(df_subset$country_code) == "CPV-SA"] <- "CPV"

order_country_code <- unique(as.vector(df_subset$country_code))

df_subset$country_code <- factor(df_subset$country_code, 
                                 levels = order_country_code, ordered = TRUE)


median_meth <- median(df_subset$percent_methylation)

# Check row number for Cvi-0 from 1001
which(grepl(6911, df_subset$seq_ID))

# Check our Cvi-0
which(grepl("4073_M", df_subset$seq_ID))
# 84

# The two Moroccan from the 1001
which(grepl(9606, df_subset$seq_ID))
# 90
which(grepl(9939, df_subset$seq_ID))
# 91

# UKID 116 (5822)
which(grepl(5822, df_subset$seq_ID))
# 549

# Dor-10 (5856)
which(grepl(5856, df_subset$seq_ID))
# 625

# I get overlying grey and black dots for CPV and FO

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggplot(data = df_subset, aes(x = country_code, y = percent_methylation)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("mCHH all TEs") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("% of methylated cytosines (CHH)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Populations sorted by increasing latitude") +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  gghighlight(country_code %in% c("CPV", "MAR"), use_direct_label = FALSE) +
  geom_hline(yintercept = median(df_subset$percent_methylation), linetype = "dashed", color = "black", size = 0.5) +
  geom_jitter(data = df_subset[107:629, ], colour = "#d2d2d2ff", size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[1:106, ], colour = "black", size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[c(1, 84, 90, 91), ], colour = mycol, size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[1, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[84, ], size = 2, fill = "red", shape = 23, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[90, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[91, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0)
```

![](images/mCHH_all_TEs_worldwide.png)

### mCG long TEs

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
# ggplot_all(df_mean_genes, title = title)

df <- merge(df_mean_TEs_4kb, df_accessions, by = "sample")

# Use dplyr prefix as desc function is also in IRange, which creates a conflict
df_subset <- df %>%
  group_by(country_code) %>%
  filter(context == "CpG", country_code != "USA", name != "SRR771702") %>%
  filter(n() > 7) %>%
  arrange(-dplyr:::desc(Latitude))

df_subset <- df %>%
  group_by(country_code) %>%
  filter(context == "CpG", !(country_code %in% c("USA", "CPV-FO")), name != "SRR771702") %>%
  filter(n() > 7) %>%
  arrange(-dplyr:::desc(Latitude))

# Change CPV-SA to CPV
levels(df_subset$country_code)[levels(df_subset$country_code) == "CPV-SA"] <- "CPV"

order_country_code <- unique(as.vector(df_subset$country_code))

df_subset$country_code <- factor(df_subset$country_code, levels = order_country_code, ordered = TRUE)

median_meth <- median(df_subset$percent_methylation)

# Check row number for Cvi-0 from 1001
which(grepl(6911, df_subset$seq_ID))
1

# Check our Cvi-0
which(grepl("4073_M", df_subset$seq_ID))
# 84

# The two Moroccan from the 1001
which(grepl(9606, df_subset$seq_ID))
# 90
which(grepl(9939, df_subset$seq_ID))
# 91

# UKID 116 (5822)
which(grepl(5822, df_subset$seq_ID))
# 549

# Dor-10 (5856)
which(grepl(5856, df_subset$seq_ID))
# 625

# I get overlying grey and black dots for CPV and FO

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggplot(data = df_subset, aes(x = country_code, y = percent_methylation)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("mCG long TEs") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("% of methylated cytosines (CG)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Populations sorted by increasing latitude") +
  scale_y_continuous(limits = c(70, 90), breaks = seq(70, 90, 5)) +
  gghighlight(country_code %in% c("CPV", "MAR"), use_direct_label = FALSE) +
  geom_hline(yintercept = median(df_subset$percent_methylation), linetype = "dashed", color = "black", size = 0.5) +
  geom_jitter(data = df_subset[107:629, ], colour = "#d2d2d2ff", size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[1:106, ], colour = "black", size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[c(1, 84, 90, 91), ], colour = mycol, size = 0.5, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[1, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[84, ], size = 2, fill = "red", shape = 23, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[90, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0) +
  geom_jitter(data = df_subset[91, ], size = 2, fill = "#d2d2d2ff", shape = 21, width = 0.1, height = 0)
```

![](images/mCG_long_TEs_worldwide.png)

### mCHG long TEs

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
#ggplot_all(df_mean_genes, title = title)

df <- merge(df_mean_TEs_4kb, df_accessions, by="sample")

# Use dplyr prefix as desc function is also in IRange, which creates a conflict
df_subset <- df %>% group_by(country_code) %>% filter(context=="CHG", country_code!="USA", name!="SRR771702") %>% filter(n() > 7) %>%  arrange(-dplyr:::desc(Latitude))

df_subset <- df %>% group_by(country_code) %>% filter(context=="CHG", !(country_code %in% c("USA","CPV-FO")), name!="SRR771702") %>% filter(n() > 7) %>%  arrange(-dplyr:::desc(Latitude))

# Change CPV-SA to CPV
levels(df_subset$country_code)[levels(df_subset$country_code)=="CPV-SA"] <- "CPV"

order_country_code <- unique(as.vector(df_subset$country_code))

df_subset$country_code <- factor(df_subset$country_code, levels=order_country_code, ordered=TRUE)

median_meth <- median(df_subset$percent_methylation)

# Check row number for Cvi-0 from 1001
which(grepl(6911, df_subset$seq_ID))
1

# Check our Cvi-0
which(grepl("4073_M", df_subset$seq_ID))
# 84

# The two Moroccan from the 1001
which(grepl(9606, df_subset$seq_ID))
# 90
which(grepl(9939, df_subset$seq_ID))
# 91

# UKID 116 (5822)
which(grepl(5822, df_subset$seq_ID))
# 549

# Dor-10 (5856)
which(grepl(5856, df_subset$seq_ID))
# 625

# I get overlying grey and black dots for CPV and FO

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggplot(data=df_subset, aes(x=country_code, y=percent_methylation)) + geom_boxplot(outlier.shape = NA) + ggtitle("mCHG long TEs") + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + ylab("% of methylated cytosines (CHG)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Populations sorted by increasing latitude") +
  scale_y_continuous(limits=c(30,70), breaks=seq(30,70,5)) +
  gghighlight(country_code %in% c("CPV","MAR"), use_direct_label = FALSE) +
  geom_hline(yintercept=median(df_subset$percent_methylation), linetype="dashed", color = "black", size=0.5) + geom_jitter(data=df_subset[107:629,], colour="#d2d2d2ff", size=0.5, width=0.1, height=0) + geom_jitter(data=df_subset[1:106,], colour="black", size=0.5, width=0.1, height=0) +
  geom_jitter(data=df_subset[c(1,84,90,91),], colour=mycol, size=0.5, width=0.1, height=0) +
  geom_jitter(data=df_subset[1,],size=2, fill="#d2d2d2ff", shape=21, width=0.1, height=0) +
  geom_jitter(data=df_subset[84,],size=2, fill="red", shape=23, width=0.1, height=0) +
  geom_jitter(data=df_subset[90,],size=2, fill="#d2d2d2ff", shape=21, width=0.1, height=0) +
  geom_jitter(data=df_subset[91,],size=2, fill="#d2d2d2ff", shape=21, width=0.1, height=0)
```

![](images/mCHG_long_TEs_worldwide.png)

### mCHH long TEs

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
#ggplot_all(df_mean_genes, title = title)
df <- merge(df_mean_TEs_4kb, df_accessions, by="sample")

# Use dplyr prefix as desc function is also in IRange, which creates a conflict
df_subset <- df %>% group_by(country_code) %>% filter(context=="CHH", country_code!="USA", name!="SRR771702") %>% filter(n() > 7) %>%  arrange(-dplyr:::desc(Latitude))

df_subset <- df %>% group_by(country_code) %>% filter(context=="CHH", !(country_code %in% c("USA","CPV-FO")), name!="SRR771702") %>% filter(n() > 7) %>%  arrange(-dplyr:::desc(Latitude))

# Change CPV-SA to CPV
levels(df_subset$country_code)[levels(df_subset$country_code)=="CPV-SA"] <- "CPV"

order_country_code <- unique(as.vector(df_subset$country_code))

df_subset$country_code <- factor(df_subset$country_code, levels=order_country_code, ordered=TRUE)

median_meth <- median(df_subset$percent_methylation)

# Check row number for Cvi-0 from 1001
which(grepl(6911, df_subset$seq_ID))
1

# Check our Cvi-0
which(grepl("4073_M", df_subset$seq_ID))
# 84

# The two Moroccan from the 1001
which(grepl(9606, df_subset$seq_ID))
# 90
which(grepl(9939, df_subset$seq_ID))
# 91

# UKID 116 (5822)
which(grepl(5822, df_subset$seq_ID))
# 549

# Dor-10 (5856)
which(grepl(5856, df_subset$seq_ID))
# 625

# I get overlying grey and black dots for CPV and FO

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggplot(data=df_subset, aes(x=country_code, y=percent_methylation)) + geom_boxplot(outlier.shape = NA) + ggtitle("mCHH long TEs") + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + ylab("% of methylated cytosines (CHH)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Populations sorted by increasing latitude") +
  scale_y_continuous(limits=c(0,25), breaks=seq(0,25,5)) +
  gghighlight(country_code %in% c("CPV","MAR"), use_direct_label = FALSE) +
  geom_hline(yintercept=median(df_subset$percent_methylation), linetype="dashed", color = "black", size=0.5) + geom_jitter(data=df_subset[107:629,], colour="#d2d2d2ff", size=0.5, width=0.1, height=0) + geom_jitter(data=df_subset[1:106,], colour="black", size=0.5, width=0.1, height=0) +
  geom_jitter(data=df_subset[c(1,84,90,91),], colour=mycol, size=0.5, width=0.1, height=0) +
  geom_jitter(data=df_subset[1,],size=2, fill="#d2d2d2ff", shape=21, width=0.1, height=0) +
  geom_jitter(data=df_subset[84,],size=2, fill="red", shape=23, width=0.1, height=0) +
  geom_jitter(data=df_subset[90,],size=2, fill="#d2d2d2ff", shape=21, width=0.1, height=0) +
  geom_jitter(data=df_subset[91,],size=2, fill="#d2d2d2ff", shape=21, width=0.1, height=0)
```

![](images/mCHH_long_TEs_worldwide.png)

# GWAS analysis

## Download raw fastq reads whole-genome sequencing

The paired-end sequencing data of the 190 CVI accessions from Santo
Antao were published in Fulgione et al., 2022 (DOI:
10.1038/s41467-022-28800-z) and are available in ENA (accession code
PRJEB39079 <https://www.ebi.ac.uk/ena/browser/view/PRJEB39079>).

Data can be downloaded using fasterq-dump from SRA toolkit:

``` bash

# Install SRA toolkit if needed (v3.2.1 used)
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

tar -vxzf sratoolkit.tar.gz
# Add path to binaries to .bashrc
#export PATH="${PATH}:/usr/users/zicola/bin/sratoolkit.3.2.1-ubuntu64/bin"

fasterq-dump --version
fasterq-dump : 3.2.1
```

``` bash
#!/bin/bash
#
#SBATCH --job-name=download_fastq
#SBATCH --partition=medium
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=6
#SBATCH --output=slurm/%x.%A.%a.out
#SBATCH --error=slurm/%x.%A.%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=johan.zicola@uni-goettingen.de
#SBATCH --array=1-190

accession_ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p WGS_ENA_data_190_SA.txt | cut -f1)

fasterq-dump $accession_ID --threads 6 --split-files -O raw_fastq/
```

``` bash
sbatch run_download.sh
```

-S: skip technical reads -e 6: use 6 threads –split-files: splits into
\_1.fastq and \_2.fastq -O output_dir/: specify output directory

It took about 10 to 40 minutes per sample to download.

## SNP calling

The SNP calling is described in the GitHub repository
<https://github.com/johanzi/SNP_calling_Arabidopsis>

## SNP annotation in CVI

Use SNPeff to annotate variants present within Santo Antao island. We
have WGS data for 190 accessions.

## Prepare VCF file for the 83 CVI accession

The VCF file `superVcf_19-07-04_cvis.vcf.b.gz_snps.vcf.b.gz` (matching
the `merged.vcf.gz` described before) was generated with the Shore
pipeline ([Ossowski et al.,
2008](http://genome.cshlp.org/content/18/12/2024)). The file is 11 Gb
and therefore available only upon request.

``` bash
# Reorder seqID
cut -f2 data/df_accessions_83.txt | grep -v "seq_ID" | sort > data/df_accessions_83_seqID_sorted.txt

VCF="/srv/biodata/irg/grp_hancock/VCF/superVcf_19-07-04_cvis.vcf.b.gz_snps.vcf.b.gz"

# Keep the 83 accessions from Santo Antao and keep only chromosomes
bcftools view -S df_accessions_83_seqID_sorted.txt -r Chr1,Chr2,Chr3,Chr4,Chr5 $VCF > subset_83_only_chr.vcf

# Remove positions without alternative alleles and non-biallelic SNPs
bcftools view --min-ac=1 --max-alleles 2  subset_83_only_chr.vcf > subset_83_only_chr_biallelic_only_alt.vcf

# Filter SNP quality with DP>=3, GQ>=25
vcftools --vcf subset_83_only_chr_biallelic_only_alt.vcf  \
            --minDP 3 --minGQ 25 --recode --recode-INFO-all \
            --out subset_83_only_chr_biallelic_only_alt_DP3_GQ25

# Mark singletons
vcftools --singletons --vcf subset_83_only_chr_biallelic_only_alt_DP3_GQ25.recode.vcf

# Remove singletons
vcftools --vcf subset_83_only_chr_biallelic_only_alt_DP3_GQ25.recode.vcf \
            --exclude-positions out.singletons --recode --recode-INFO-all \
            --out subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons
            
# Compress and tabix file
bgzip subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf \
  && tabix subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf.gz
```

The output file is be
`subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf.gz`.

## Prepare phenotype

Check
<https://github.com/johanzi/gwas_gemma?tab=readme-ov-file#section-id-139>

### Generate phenotype files

In bash, the rda files generated by the function `load_df_wml` is used
to generate an equivalent txt files. Here is the pipeline in bash

The script `generate_phenotype_gemma.sh` takes as arguments an rda
object (generated by `load_df_wml`), the name of the methylation
context, the region (gene, TE, …) and a file containing the accession
order as present in the vcf file. Check GitHub
[documentation](https://github.com/HancockLab/WGBS_African_Arabidopsis)
for more information.

``` bash

# Generate phenotypes

# Working directory
cd /srv/netscratch/irg/grp_hancock/johan/GWAS/dna_methylation/GWAS_83_SA

# Need to add 'sample' (BS-seq library underscore accession name) variable to df_accessions_83.txt
#awk -v FS='\t' -v OFS='\t' '{print $0,$3"_"$1}' df_accessions_83.txt > df_accessions_83_sample.txt
#sed -i 's/library_name/sample/' df_accessions_83_sample.txt
#mv df_accessions_83_sample.txt df_accessions_83.txt

# Path variables
path_script="/home/zicola/SCRIPTS/bismark_pipeline/WGBS_African_Arabidopsis"

path_rda="/srv/netscratch/irg/grp_hancock/johan/methylKit_DB_files/GC_3427_3542_3599_4050_4220_4373_TAIR10"

# Generate phenotype for whole genome
bash scripts/generate_phenotype_gemma.sh data/df_mean_filtered_83.rda CpG whole_genome data/df_accessions_83_seqID_sorted.txt data/df_accessions_83.txt

bash scripts/generate_phenotype_gemma.sh data/df_mean_filtered_83.rda CHG whole_genome data/df_accessions_83_seqID_sorted.txt data/df_accessions_83.txt

bash scripts/generate_phenotype_gemma.sh data/df_mean_filtered_83.rda CHH whole_genome data/df_accessions_83_seqID_sorted.txt data/df_accessions_83.txt


# Generate phenotype for genes
bash scripts/generate_phenotype_gemma.sh data/df_mean_genes_83.rda CpG genes data/df_accessions_83_seqID_sorted.txt 

bash scripts/generate_phenotype_gemma.sh data/df_mean_genes_83.rda CHG genes data/df_accessions_83_seqID_sorted.txt 

bash scripts/generate_phenotype_gemma.sh data/df_mean_genes_83.rda CHH genes data/df_accessions_83_seqID_sorted.txt

# Generate phenotype for TEs
bash scripts/generate_phenotype_gemma.sh data/df_mean_TEs_83.rda CpG TEs data/df_accessions_83_seqID_sorted.txt

bash scripts/generate_phenotype_gemma.sh data/df_mean_TEs_83.rda CHG TEs data/df_accessions_83_seqID_sorted.txt

bash scripts/generate_phenotype_gemma.sh data/df_mean_TEs_83.rda CHH TEs data/df_accessions_83_seqID_sorted.txt

# Generate phenotype for long TEs
bash scripts/generate_phenotype_gemma.sh data/df_mean_TEs_4kb_83.rda CpG data/TEs_4kb df_accessions_83_seqID_sorted.txt

bash scripts/generate_phenotype_gemma.sh data/df_mean_TEs_4kb_83.rda CHG data/TEs_4kb df_accessions_83_seqID_sorted.txt

bash scripts/generate_phenotype_gemma.sh data/df_mean_TEs_4kb_83.rda CHH data/TEs_4kb df_accessions_83_seqID_sorted.txt
```

## Run Gemma

GWAS was performed as described in this repository:
<https://github.com/johanzi/gwas_gemma>

The script should be executed for each context and each genomic region:

``` bash
cd data

VCF="subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf.gz"

# Whole genome
bash ../scripts/run_gwas_gemma.sh CpG_whole_genome.tsv $VCF

bash ../scripts/run_gwas_gemma.sh CHG_whole_genome.tsv $VCF

bash ../scripts/run_gwas_gemma.sh CHH_whole_genome.tsv $VCF

# Genes
bash ../scripts/run_gwas_gemma.sh CpG_genes.tsv $VCF

bash ../scripts/run_gwas_gemma.sh CHG_genes.tsv $VCF

bash ../scripts/run_gwas_gemma.sh CHH_genes.tsv $VCF

# TEs
bash ../scripts/run_gwas_gemma.sh CpG_TEs.tsv $VCF

bash ../scripts/run_gwas_gemma.sh CHG_TEs.tsv $VCF

bash ../scripts/run_gwas_gemma.sh CHH_TEs.tsv $VCF

# long TEs
bash ../scripts/run_gwas_gemma.sh CpG_TEs_4kb.tsv $VCF

bash ../scripts/run_gwas_gemma.sh CHG_TEs_4kb.tsv $VCF

bash ../scripts/run_gwas_gemma.sh CHH_TEs_4kb.tsv $VCF
```

### GWAS whole genome

Load R functions for GWAS (originally from
<https://raw.githubusercontent.com/johanzi/gwas_gemma/refs/heads/master/GWAS_run.R>)

``` r
source("scripts/functions_gwas.R")
```

#### mCG

``` r
dir_file <- "data/output/"

file.name <- "CpG_whole_genome.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")

SNP_to_BED(SNP_significant, "data/mCG_whole_genome_GWAS_SNPs.bed")
```

![](images/GWAS_mCG_whole_genome.png)

#### mCHG

``` r
dir_file="data/output/"

file.name <- "CHG_whole_genome.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCHG_whole_genome.png)

#### mCHH

``` r
dir_file="data/output/"

file.name <- "CHH_whole_genome.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCHH_whole_genome.png)

### GWAS genes

#### mCG

``` r
dir_file="data/output/"

file.name <- "CpG_genes.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")

SNP_to_BED(SNP_significant, "data/gbM_GWAS_SNPs.bed")
```

![](images/GWAS_mCG_genes.png)

##### mCG with VIM annotation

``` r
source("scripts/functions_gwas.R") 

file.name <- "CpG_genes.assoc.clean.txt"

dir_file="data/output/"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

# Check code https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R

# Won't work if no SNPs within the gene region, make a larger range to get some SNPs to be displayed

VIM1 <- c(21413981,21418115)
VIM2 <- c(24589343,24592780)
VIM3 <- c(15837178,15840678)


# Add n bp before and after the gene coordinates to widen the regions to find. I played by adding 0 until I could see the 4 genes on the plot
n = 50000
VIM1 <- c(VIM1[[1]]-n, VIM1[[2]]+n)
VIM2 <- c(VIM2[[1]]-n, VIM2[[2]]+n)
VIM3 <- c(VIM3[[1]]-n, VIM3[[2]]+n)

ann<-rep(1, length(gwas.results$P))
ann[with(gwas.results, CHR==1 & BP>=VIM1[[1]] & BP<VIM1[[2]])]<-2
ann[with(gwas.results, CHR==1 & BP>=VIM2[[1]] & BP<VIM2[[2]])]<-3
ann[with(gwas.results, CHR==5 & BP>=VIM3[[1]] & BP<VIM3[[2]])]<-4
ann<-factor(ann, levels=1:4, labels=c("","VIM1","VIM2", "VIM3"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title="GWAS gbM")
```

![](images/GWAS_gbM_VIMs.png)

Top SNP on Chr5 is near VIM3 (790 kb upstream) and falls within the gene
AT5G37810.

![](images/chr5_15047549.png)

#### mCHG

``` r
dir_file="data/output/"

file.name <- "CHG_genes.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCHG_genes.png)

#### mCHH

``` r
dir_file="data/output/"

file.name <- "CHH_genes.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCHH_genes.png)

### GWAS all TEs

#### mCG

``` r
dir_file <- "data/output/"

file.name <- "CpG_TEs.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCG_TEs.png)

#### mCHG

``` r
dir_file="data/output/"

file.name <- "CHG_TEs.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCHG_TEs.png)

``` r
source("scripts/functions_gwas.R") 

dir_file <- "data/output/"

file.name <- "CHG_TEs_4kb.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

SUVH4 <- c(4501362,4506250)
AGO9 <- c(7192240,7198381)
MET1 <- c(19932114,19938470)
FBX5 <- c(18513626)
DRM1 <- c(4991347,4994924)

# Add n bp before and after the gene coordinates to widen the regions to find. I played by adding 0 until I could see the 4 genes on the plot
n = 50000
SUVH4 <- c(SUVH4[[1]]-n, SUVH4[[2]]+n)
AGO9 <- c(AGO9[[1]]-n, AGO9[[2]]+n)
MET1 <- c(MET1[[1]]-n, MET1[[2]]+n)
DRM1 <- c(DRM1[[1]]-n, DRM1[[2]]+n)

ann<-rep(1, length(gwas.results$P))
ann[with(gwas.results, CHR==5 & BP>=SUVH4[[1]] & BP<SUVH4[[2]])]<-2
ann[with(gwas.results, CHR==5 & BP>=AGO9[[1]] & BP<AGO9[[2]])]<-3
ann[with(gwas.results, CHR==5 & BP>=MET1[[1]] & BP<MET1[[2]])]<-4
ann[with(gwas.results, CHR==2 & BP==FBX5)]<-5
ann[with(gwas.results, CHR==5 & BP>=DRM1[[1]] & BP<DRM1[[2]])]<-6
ann<-factor(ann, levels=1:6, labels=c("","SUVH4","AGO9","MET1","FBX5","DRM1"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title="GWAS mCHG at long TEs")
```

#### mCHH

``` r
dir_file="data/output/"

file.name <- "CHH_TEs.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCHH_TEs.png)

### GWAS long TEs

#### mCG

``` r
dir_file <- "data/output/"

file.name <- "CpG_TEs_4kb.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCG_long_TEs.png)

Manhattan plot with FBX5 SNP highlighted

``` r
source("scripts/functions_gwas.R") 

dir_file <- "data/output/"

file.name <- "CpG_TEs_4kb.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

FBX5 <- c(18513626)

ann<-rep(1, length(gwas.results$P))
ann[with(gwas.results, CHR==2 & BP==FBX5)]<-2
ann<-factor(ann, levels=1:2, labels=c("","FBX5"))

#png("GWAS_test.png", width = 3000, height = 500, units = "px", res=300)

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title="GWAS for mCG at long TEs")

#dev.off()
```

![](images/GWAS_mCG_long_TEs_FBX5_SNP.png)

``` r
manhattan.plot(chr = gwas.results$CHR, 
               pos=gwas.results$BP, pvalue=gwas.results$P, 
               annotate=ann, sig.level=bonferroni_threshold, title="GWAS for mCG at long TEs")

ggsave(filename = "GWAS_test.png", width = 13, height = 4, units = "cm", dpi=300)
```

#### mCHG

``` r
dir_file="data/output/"

file.name <- "CHG_TEs_4kb.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCHG_long_TEs.png)

##### Highlight FBX5 SNP

``` r
source("scripts/functions_gwas.R")

dir_file="data/output/"

file.name <- "CHG_TEs_4kb.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

FBX5 <- c(18513626)

ann<-rep(1, length(gwas.results$P))
ann[with(gwas.results, CHR==2 & BP==FBX5)]<-2
ann<-factor(ann, levels=1:2, labels=c("","FBX5"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title="GWAS mCHG at long TEs")
```

![](images/GWAS_mCHG_long_TEs_FBX5.png)

Recover SNPs above -log10(p) = 3

``` r
SNP_mCHG_long_TEs <- GWAS_run(path.file, threshold_pvalue = 10e-4)
# 50 SNPs recovered

SNP_to_BED(SNP_mCHG_long_TEs, "data/mCHG_long_TEs_GWAS_SNPs.bed")
```

The peak at Chr5 seems to be close to SUVH4. Most significant SNP is
5:5843297

##### Highlight FBX5, SUVH4, AGO9, DRM1, and MET1

``` r
source("scripts/functions_gwas.R") 

dir_file="data/output/"

file.name <- "CHG_TEs_4kb.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

SUVH4 <- c(4501362,4506250)
AGO9 <- c(7192240,7198381)
MET1 <- c(19932114,19938470)
FBX5 <- c(18513626)
DRM1 <- c(4991347,4994924)

# Add n bp before and after the gene coordinates to widen the regions to find. I played by adding 0 until I could see the 4 genes on the plot
n = 50000
SUVH4 <- c(SUVH4[[1]]-n, SUVH4[[2]]+n)
AGO9 <- c(AGO9[[1]]-n, AGO9[[2]]+n)
MET1 <- c(MET1[[1]]-n, MET1[[2]]+n)
DRM1 <- c(DRM1[[1]]-n, DRM1[[2]]+n)

ann<-rep(1, length(gwas.results$P))
ann[with(gwas.results, CHR==5 & BP>=SUVH4[[1]] & BP<SUVH4[[2]])]<-2
ann[with(gwas.results, CHR==5 & BP>=AGO9[[1]] & BP<AGO9[[2]])]<-3
ann[with(gwas.results, CHR==5 & BP>=MET1[[1]] & BP<MET1[[2]])]<-4
ann[with(gwas.results, CHR==2 & BP==FBX5)]<-5
ann[with(gwas.results, CHR==5 & BP>=DRM1[[1]] & BP<DRM1[[2]])]<-6
ann<-factor(ann, levels=1:6, labels=c("","SUVH4","AGO9","MET1","FBX5","DRM1"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title="GWAS for mCHG in long TEs")
```

![](images/GWAS_mCHG_long_TEs_genes_highlighted.png)

``` r
SNP_mCHG_long_TEs <- GWAS_run(path.file, threshold_pvalue = 10e-5)
# 50 SNPs recovered

SNP_to_BED(SNP_mCHG_long_TEs, "data/mCHG_long_TEs_GWAS_SNPs.bed")
```

#### mCHH

``` r
dir_file="data/output/"

file.name <- "CHH_TEs_4kb.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

SNP_significant <- GWAS_run(path.file, threshold_pvalue = "bonferroni")
```

![](images/GWAS_mCHH_long_TEs.png)

## Variants at SUVH4, AGO9, DRM1, and MET1

We used SnpEff (<https://pcingola.github.io/SnpEff/>) to define variants
and their effect of the 4 candidate genes in Chr5 (AGO9, SUVH4, DRM1,
and MET1).

### Subsample VCF file for the 190 Santo Antao accessions

Note: We include Cvi-0 as Santo Antao accession as it is genetically
closer to Santo Antao accessions than Fogo accessions. Cvi-0 was
collected in 1983 by Lobin in Cape Verde but the island was not
specified. See
<https://www.arabidopsis.org/ais/1983/lobin-1983-aabme.html>

``` bash

module load bcftools vcftools

# Keep the 190 accessions from Santo Antao and keep only chromosomes
bcftools view -S list_190_SA_accessions.txt -r Chr1,Chr2,Chr3,Chr4,Chr5 \
    superVcf_19-07-04_cvis.vcf.b.gz_snps.vcf.b.gz > subset_190_SA.vcf
    
# Remove positions without alternative alleles and non-biallelic SNPs
bcftools view --min-ac=1 --max-alleles 2  subset_190_SA.vcf > subset_190_SA_alt.vcf

# Filter SNP quality with DP>=3, GQ>=25, remove all SNP with more than 10% missing data
vcftools --vcf subset_190_SA_alt.vcf --minDP 3 --minGQ 25 --max-missing 0.9 \
  --recode --recode-INFO-all --out subset_190_SA_alt_DP3_GQ25

#After filtering, kept 484836 out of a possible 1370456 Sites
#Run Time = 222.00 seconds

# Mark singletons
vcftools --singletons --vcf subset_190_SA_alt_DP3_GQ25.recode.vcf

wc -l out.singletons
# 5743 out.singletons

# Remove singletons
vcftools --vcf subset_190_SA_alt_DP3_GQ25.recode.vcf \
  --exclude-positions out.singletons --recode \
  --recode-INFO-all --out subset_190_SA_alt_DP3_GQ25_wo_singletons

#After filtering, kept 479094 out of a possible 484836 Sites
#Run Time = 110.00 seconds

# Get allele frequency of each SNP on the SA population only
# https://speciationgenomics.github.io/filtering_vcfs/
vcftools --vcf subset_190_SA_alt_DP3_GQ25_wo_singletons.recode.vcf \
  --freq2 --out frequencies_SA

# Return file frequencies_SA.frq
# Remove all variants that are absent in SA (col 6 = 0% AF) and remove variants 
# where no data are present in SA (returns -nan error value when doing the division)
awk '$4!=0 && $6!=0' frequencies_SA.frq | grep -v "\-nan" > list_SNP_ALT_more_than_0.txt

wc -l list_SNP_ALT_more_than_0.txt
442725 list_SNP_ALT_more_than_0.t

# Retrieve SNPs with ALT present in SA
vcftools --vcf subset_190_SA_alt_DP3_GQ25_wo_singletons.recode.vcf \
  --positions list_SNP_ALT_more_than_0.txt --recode \
  --recode-INFO-all --out subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only
#After filtering, kept 442724 out of a possible 479094 Sites
#Run Time = 104.00 seconds

# Compress this file for further use in other analysis
bcftools view subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only.recode.vcf -Oz -o subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only.recode.vcf.gz
bcftools index subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only.recode.vcf.gz
```

Outpuf VCF file is
`subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only.recode.vcf`

### Run SnpEff

``` bash
# Download snpeff (2025-05-15)
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

unzip snpEff_latest_core.zip

cd snpEff

# Check version
java -jar snpEff.jar -version
SnpEff  5.2f    2025-02-07

# Install Arabidopsis database (Java 24.0.1 is used)
java -jar snpEff.jar download Arabidopsis_thaliana

\ls -1 data/Arabidopsis_thaliana/
sequence.1.bin
sequence.2.bin
sequence.3.bin
sequence.4.bin
sequence.5.bin
sequence.bin
snpEffectPredictor.bin

# Run SnpEff
java -Xmx18g -jar /usr/users/zicola/bin/snpEff/snpEff.jar \
    -c /usr/users/zicola/bin/snpEff/snpEff.config \
    Arabidopsis_thaliana subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only.recode.vcf \
    > subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only.recode.ann.vcf
    
```

Took 6 min. Output file is
`subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only.recode.ann.vcf`.

The file `snpEff_genes.txt` summarize the number of variants per
transcript and variant category. I can import the text file in Excel for
easy visualization and sorting.

### Summary variants

``` bash
# Recover only a subset of the output for easier parsing and reading
grep -v "##" subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only.recode.ann.vcf | \
  cut -f1,2,4,5,8 | cut -d'|' -f1,2,3,4,5,10,11 | \
  sed 's/|/\t/g' > summary_SnpEff.txt
```

List of candidate genes at chr5:

    AT5G49160
    AT5G21150
    AT5G13960
    AT5G15380

Put gene list in `genes_chr5_mCHG_long_TEs.txt`

``` bash

while read i; do
  grep $i summary_SnpEff.txt
done > summary_SnpEff_chr5_candidates.txt < genes_chr5_mCHG_long_TEs.txt

# Get frequencies
while read i; do
  pos=$(echo "$i" | cut -f2)
  grep -P -w "Chr5\t$pos" frequencies_SA.frq
done < summary_SnpEff_chr5_candidates.txt | cut -f5,6 > summary_SnpEff_chr5_candidates_freq.txt

# Merge summary and frequencies
paste summary_SnpEff_chr5_candidates.txt summary_SnpEff_chr5_candidates_freq.txt > summary_SnpEff_chr5_candidates_final.txt

rm summary_SnpEff_chr5_candidates.txt summary_SnpEff_chr5_candidates_freq.txt
```

Import text file `summary_SnpEff_chr5_candidates_final.txt` in Excel
(Data \> From Text \> Tab-separated) =\> Supplementary Table 11.

# Allele status in CPV

``` bash

# This script defines the step to retrieve the allele status for
# each CPV-SA accessions of the CMT2, FBX5, and VIM2 variants

################################################################################
# Packages/Sofwares requires
################################################################################

# vcf_melt (install with command "pip install PyVCF --user")

# bcftools

# samtools

################################################################################
########## Get VIM2 deletion, FBX5, and CMT2 status
################################################################################

# Use VCF file used for GWAS 

VCF="subset_190_SA_alt_DP3_GQ25_wo_singletons_SA_only.recode.vcf.gz"

#################################################################
# FBX5
#################################################################

# Get VCF data for SNP in FBX5
bcftools view -r Chr2:18513626 $VCF > Chr2_18513626.vcf

# Convert into vertical
vcf_melt Chr2_18513626.vcf > Chr2_18513626.melted.vcf 

# Keep only line with GQ >= 25 and DP >= 3
awk '$3>=25 && $4>=3 {print $0}' Chr2_18513626.melted.vcf > Chr2_18513626_GQ25_DP3.melted.vcf
 
# Accessions with alternative allele
awk '$2 == "1" {print $0}' Chr2_18513626_GQ25_DP3.melted.vcf | wc -l
84

awk '$2 == "1" {print $0}' Chr2_18513626_GQ25_DP3.melted.vcf | cut -f1 | sort - > FBX5_alt.txt
awk '$2 == "0" {print $0}' Chr2_18513626_GQ25_DP3.melted.vcf | cut -f1 | sort - > FBX5_ref.txt

awk -v OFS='\t' '{print $0,"FBX5_alt"}' FBX5_alt.txt > FBX5_alt_final.txt
awk -v OFS='\t' '{print $0,"FBX5_ref"}' FBX5_ref.txt > FBX5_ref_final.txt

cat FBX5_alt_final.txt FBX5_ref_final.txt > FBX5_allele_status.txt

rm FBX5_alt.txt FBX5_ref.txt FBX5_alt_final.txt FBX5_ref_final.txt 

#################################################################
# CMT2
#################################################################

bcftools view -r Chr4:10420088 $VCF > Chr4_10420088.vcf

# Convert into vertical
vcf_melt Chr4_10420088.vcf > Chr4_10420088.melted.vcf 

# Keep only line with GQ >= 25 and DP >= 3
awk '$3>=25 && $4>=3 {print $0}' Chr4_10420088.melted.vcf  > Chr4_10420088_GQ25_DP3.melted.vcf 

# Accessions with alternative allele
awk '$2 == "1" {print $0}' Chr4_10420088_GQ25_DP3.melted.vcf  | wc -l
65

awk '$2 == "1" {print $0}' Chr4_10420088_GQ25_DP3.melted.vcf  | cut -f1 | sort - > CMT2_alt.txt
awk '$2 == "0" {print $0}' Chr4_10420088_GQ25_DP3.melted.vcf  | cut -f1 | sort - > CMT2_ref.txt

awk -v OFS='\t' '{print $0,"CMT2_alt"}' CMT2_alt.txt > CMT2_alt_final.txt
awk -v OFS='\t' '{print $0,"CMT2_ref"}' CMT2_ref.txt > CMT2_ref_final.txt

cat CMT2_alt_final.txt CMT2_ref_final.txt > CMT2_allele_status.txt

rm CMT2_alt.txt CMT2_ref.txt CMT2_alt_final.txt CMT2_ref_final.txt 

#################################################################
# VIM2 deletion
#################################################################

# Considering that the VIM2 deletion is not present in the VCF file as it is a 
# structural variant and not a SNP, we need to assess the presence of the deletion
# based on read density at the deletion region

# coordinates of the deletion location (based on Cvi-0). The deletion is 2740 bp
# This region was defined by looking at read mapping in Cvi-0 and determine visually
# the beginning and end of the deletion
coordinates="chr1:24586731-24589471"

for i in mappedBAM/*bam; do
    name=$(basename $i | cut -d'.' -f1)
    nb_reads=$(samtools view $i $coordinates | wc -l)
    echo -e "${name}\t${nb_reads}" >> nb_reads_vim2_deletion.txt
done

cut -f2 nb_reads_vim2_deletion.txt | sort -n -
# Looking at the distribution of reads, it seems there is threshold at 14 reads. i
# Let's use 50 reads as the threshold to define that there is indeed a deletion

# Classify each sample based on nb of reads with threshold = 50
# If a given accession has 50 or less reads mapping at the promoter region of VIM2-VIM4, the accession is considered
# as having the VIM2 deletion.
awk -v OFS="\t" '$2 <= 50 {print $1,$2,"deletion"} $2 > 50 {print $1,$2,"no_deletion"}' nb_reads_vim2_deletion.txt > nb_reads_vim2_deletion_status.txt

#################################################################
# "VIM3" SNP

# Get VCF data for SNP
bcftools view -r Chr5:15047549 $VCF > Chr5_15047549.vcf

# Convert into vertical
vcf_melt Chr5_15047549.vcf > Chr5_15047549.melted.vcf

# Keep only line with GQ >= 25 and DP >= 3
awk '$3>=25 && $4>=3 {print $0}' Chr5_15047549.melted.vcf > Chr5_15047549_GQ25_DP3.melted.vcf

# Accessions with alternative allele
awk '$2 == "1" {print $0}' Chr5_15047549_GQ25_DP3.melted.vcf | wc -l
84

awk '$2 == "1" {print $0}' Chr5_15047549_GQ25_DP3.melted.vcf | cut -f1 | sort - > VIM3_alt.txt
awk '$2 == "0" {print $0}' Chr5_15047549_GQ25_DP3.melted.vcf | cut -f1 | sort - > VIM3_ref.txt

awk -v OFS='\t' '{print $0,"VIM3_alt"}' VIM3_alt.txt > VIM3_alt_final.txt
awk -v OFS='\t' '{print $0,"VIM3_ref"}' VIM3_ref.txt > VIM3_ref_final.txt

cat VIM3_alt_final.txt VIM3_ref_final.txt > VIM3_allele_status.txt
```

## Allele distribution by population

Create the python script `scripts/find_accession/find_accession.py` to
easily find back an accession name based on its sequencing ID (seqID).

``` bash

# Make a Python dictionary of the CPV accessions names and their seqID
python ./scripts/find_accession/find_accession.py make_dict \
  ./scripts/find_accession/name_seqID_CPV_190_accessions.txt > \ 
  ./scripts/find_accession/name_seqID_CPV_190_accessions.dict

# Make a second Python dictionary of the seqID and accession names

awk '{print $2,$1}' OFS='\t' name_seqID_CPV_190_accessions.txt > seqID_name_CPV_190_accessions.txt

python ./scripts/find_accession/find_accession.py make_dict \
  ./scripts/find_accession/seqID_name_CPV_190_accessions.txt > \
  ./scripts/find_accession/seqID_name_CPV_190_accessions.dict
```

``` bash
########################################################
# Summary by population in SA for VIM2, FBX5, and CMT2
#########################################################

# 4073_M (Cvi-0 is included) => 190 accessions
clean_file="/srv/biodata/irg/grp_hancock/VCF/santos_clean_2019-07-11.txt"

#########################################################
# For FBX5

while read i; do 
    grep -w $i FBX5_allele_status.txt >> FBX5_allele_status_SA.txt
done < $clean_file

while read i; do
    seqID=$(echo "$i" | cut -f1)
    name=$(python ./scripts/find_accession/find_accession.py ./scripts/find_accession/seqID_name_CPV_190_accessions.dict $seqID | cut -f2)
    population=$(echo $name | cut -d'-' -f1)
    echo -e "${i}\t${name}\t${population}"
done < FBX5_allele_status_SA.txt > FBX5_allele_status_SA_with_names.txt


# Replace FBX5_ref by 0 and FBX5_alt by 1
sed -i 's/FBX5_ref/0/' FBX5_allele_status_SA_with_names.txt
sed -i 's/FBX5_alt/1/' FBX5_allele_status_SA_with_names.txt

# 189 retrieved, accession 12849 has no GT assigned  

cd /srv/netscratch/dep_coupland/grp_hancock/mappedBAM/CVI
samtools tview  12849.sorted.bam  -p chr2:18513626 --reference /home/zicola/TAIR10_chr_Pt_Mt/TAIR10.fasta

# Weird as many reads support the ALT allele A
cd /srv/biodata/dep_coupland/grp_hancock/johan/allele_status

# Also present after filtering but no genotype given (second column
grep "12849" Chr2_18513626_GQ25_DP3.melted.vcf
12849   .                       ['q25'] Chr2    18513626        T       [A]             1       94957   84      2003


#########################################################
# For CMT2

while read i; do 
    grep -w $i CMT2_allele_status.txt >> CMT2_allele_status_SA.txt
done < /srv/biodata/dep_coupland/grp_hancock/VCF/santos_clean_2019-07-11.txt

# 190 accessions retrieved => OK

while read i; do
    seqID=$(echo "$i" | cut -f1)
    name=$(python ./scripts/find_accession/find_accession.py ./scripts/find_accession/seqID_name_CPV_190_accessions.dict $seqID | cut -f2)
    population=$(echo $name | cut -d'-' -f1)
    echo -e "${i}\t${name}\t${population}"
done < CMT2_allele_status_SA.txt > CMT2_allele_status_SA_with_names.txt


# Replace CMT2_ref by 0 and CMT2_alt by 1
sed -i 's/CMT2_ref/0/' CMT2_allele_status_SA_with_names.txt
sed -i 's/CMT2_alt/1/' CMT2_allele_status_SA_with_names.txt


#########################################################
# For VIM2

while read i; do 
    foo=$(grep -w $i nb_reads_vim2_deletion_all_with_name_clean.txt | cut -f1,3,4)
    population=$(echo "$foo" | cut -f3 | cut -d'-' -f1)
    echo -e "${foo}\t${population}"
done < /srv/biodata/dep_coupland/grp_hancock/VCF/santos_clean_2019-07-11.txt > VIM2_allele_status_SA_with_names.txt

# 190 accessions retrieved => OK

# Replace no_deletion by 0 and deletion by 1
sed -i 's/no_deletion/0/' VIM2_allele_status_SA_with_names.txt
sed -i 's/deletion/1/' VIM2_allele_status_SA_with_names.txt

# Integrate that with coordinates data

# Reference files
/netscratch/dep_coupland/grp_hancock/Celia/Experiments/newBronson/coord_genotype.txt

# Get first 3 first rows
cut -f1,2,3 /netscratch/dep_coupland/grp_hancock/Celia/Experiments/newBronson/coord_genotype.txt > coord_populations_SA.txt
```

Note that for S3-9 (12849), reads support 50% of A and 50% of T,
suggesting heterozygosity for the FBX5 SNP Chr2:18,513,626. We excluded
this accession from the analysis based on this uncertainty.

## Map allele distribution by population

``` r
#install.packages("scatterpie")
#install.packages("ggmap")
#install.packages("maps")

library(scatterpie)
library(ggmap)
library(maps)
library(plyr)

# As input file, I need longitude, latitude, population name, frequency of ancestral allele, frequency of derived allele, gene name (optional if only 1 gene), 
# the radius to define the size of the pie charts on the map (use total_individuals  * 0.005 / 10)
```

## Coordinate file

``` r
# Replace Cratera by SCratera in bash
# Modify header
# Data from SantoCoordinates.csv files in Google Drive Hancock lab
# There are coordinates for 31 populations

# I chose arbitrarily S11-rav1  -25.076404  17.114951 as S11 population and S24-1   -25.076925  17.105766 as S24
# population, remove the other S11 and S24 from the file

df_coordinates <- read.table("data/coord_populations_SA.txt", header=TRUE)
names(df_coordinates) <- c("population","long","lat")
```

``` r
# Take coordinates from Supp Table 1 from Fulgione 2022
# Remove duplicates by pop and keep one of the two S24 coordinates

df_coordinates <- read.table("data/coord_populations_SA.txt", header=TRUE)
names(df_coordinates) <- c("population","long","lat")
```

## VIM2 distribution

``` r
df_vim2 <- read.table("data/VIM2_allele_status_SA_with_names.txt", header = FALSE)
names(df_vim2) <- c("seqID", "GT", "accession", "population")

df_vim2$population <- as.factor(df_vim2$population)

# Summarize by population
df_vim2_summary <- plyr::ddply(df_vim2, .(population), summarise, total = length(population), freq_derived = sum(GT) / total, freq_derived = (sum(GT) / total), freq_ancestral = (1 - freq_derived), radius = (total * 0.005 / 10))

# Merge with coordinate, somehow I get duplicated rows. Use unique()
df_vim2_coord <- unique(merge(df_vim2_summary, df_coordinates, by = "population"))

world <- map_data("world")
SA <- world[world$long > -30 & world$long < -20 & world$lat > 10 & world$lat < 20, ]

ggplot(SA, aes(long, lat)) +
  geom_map(map = world, aes(map_id = region), fill = NA, color = "black") +
  coord_quickmap() +
  xlim(-25.1, -25.01) +
  ylim(17.09, 17.135) +
  theme(plot.background = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line.x = element_line(color = "black", size = 0.5), axis.line.y = element_line(color = "black", size = 0.5)) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_scatterpie(aes(x = long, y = lat, r = radius, group = population), data = df_vim2_coord, cols = c("freq_derived", "freq_ancestral"), sorted_by_radius = TRUE, alpha = .5) +
  labs(fill = "VIM2 alleles") +
  scale_fill_manual(values = c("coral2", "turquoise4")) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  geom_scatterpie_legend(df_vim2_coord$radius, x = -25.03, y = 17.124, n = 4, labeller = function(x) 10 * x / 0.005) +  theme(axis.text.x = element_text(color="black"), 
      axis.text.y = element_text(color="black"),
      axis.ticks = element_line(color = "black")) + 
      theme(plot.title = element_text(hjust = 0.5))
```

![](images/VIM2_allele_map.png)

## CMT2 distribution

``` r
df_cmt2 <- read.table("data/CMT2_allele_status_SA_with_names.txt", header=FALSE)
names(df_cmt2) <- c("seqID","GT","accession","population")

# Summarize by population
df_cmt2_summary <- ddply(df_cmt2, .(population), summarise, total=length(population), freq_derived=(sum(GT)/total), freq_ancestral=(1-freq_derived), radius=(total*0.005/10))

# Merge with coordinate, somehow I get duplicated rows. Use unique()
df_cmt2_coord <- merge(df_cmt2_summary, df_coordinates, by="population")

world <- map_data("world")
SA <- world[world$long > -30 & world$long < -20 & world$lat > 10 & world$lat < 20,]

ggplot(SA, aes(long, lat)) +
  geom_map(map = world, aes(map_id = region), fill = NA, color = "black") +
  coord_quickmap() +
  xlim(-25.1, -25.01) +
  ylim(17.09, 17.135) +
  theme(plot.background = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line.x = element_line(color = "black", size = 0.5), axis.line.y = element_line(color = "black", size = 0.5)) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_scatterpie(aes(x = long, y = lat, r = radius, group = population), data = df_cmt2_coord, cols = c("freq_derived","freq_ancestral"), sorted_by_radius = TRUE, alpha = .5) +
  labs(fill = "CMT2 alleles") +
  scale_fill_manual(values = c("coral2","turquoise4")) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  geom_scatterpie_legend(df_cmt2_coord$radius, x = -25.03, y = 17.124, n = 4, labeller = function(x) 10 * x / 0.005) +  theme(axis.text.x = element_text(color="black"), 
      axis.text.y = element_text(color="black"),
      axis.ticks = element_line(color = "black")) + 
      theme(plot.title = element_text(hjust = 0.5)) 
```

![](images/CMT2_allele_map.png)

## FBX5 distribution

``` r
df_fbx5 <- read.table("data/FBX5_allele_status_SA_with_names.txt", header=FALSE)
names(df_fbx5) <- c("seqID","GT","accession","population")

# Summarize by population
df_fbx5_summary <- ddply(df_fbx5, .(population), summarise, total=length(population), freq_derived=(sum(GT)/total), freq_ancestral=(1-freq_derived), radius=(total*0.005/10))

# Merge with coordinate, somehow I get duplicated rows. Use unique()
df_fbx5_coord <- unique(merge(df_fbx5_summary, df_coordinates, by="population"))

world <- map_data("world")
SA <- world[world$long > -30 & world$long < -20 & world$lat > 10 & world$lat < 20,]

ggplot(SA, aes(long, lat)) +
  geom_map(map = world, aes(map_id = region), fill = NA, color = "black") +
  coord_quickmap() +
  xlim(-25.1, -25.01) +
  ylim(17.09, 17.135) +
  theme(plot.background = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line.x = element_line(color = "black", size = 0.5), axis.line.y = element_line(color = "black", size = 0.5)) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_scatterpie(aes(x = long, y = lat, r = radius, group = population), data = df_fbx5_coord, cols = c("freq_derived","freq_ancestral"), sorted_by_radius = TRUE, alpha = .5) +
  labs(fill = "FBX5 alleles") +
  scale_fill_manual(values = c("coral2","turquoise4")) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  geom_scatterpie_legend(df_fbx5_coord$radius, x = -25.03, y = 17.124, n = 4, labeller = function(x) 10 * x / 0.005) +  theme(axis.text.x = element_text(color="black"), 
      axis.text.y = element_text(color="black"),
      axis.ticks = element_line(color = "black")) + 
      theme(plot.title = element_text(hjust = 0.5)) 
```

![](images/FBX5_allele_map.png)

## Plot diagram VIM2

``` r
df_vim2_coord$population_nb <- paste(df_vim2_coord$population, " n=", df_vim2_coord$total, sep="")

# Remove Cvi
df <- df_vim2_coord[!(df_vim2_coord$population=="Cvi"),]

# Order accesstion by longitude
df$population_nb <- factor(df$population_nb, levels = df$population_nb[order(df$long)])

ggplot(data = df, aes(population_nb, freq_derived)) +
  ggtitle("VIM2 deletion frequency") +
  geom_bar(aes(x = population_nb, y = freq_derived), stat = "identity", colour = "black", fill = "grey") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Deletion frequency") +
  xlab("Population") +
  scale_y_continuous(labels = scales::percent)
```

![](images/VIM2_deletion_frequency.png)

## Plot diagram CMT2

``` r
df_cmt2_coord$population_nb <- paste(df_cmt2_coord$population, " n=", df_cmt2_coord$total, sep="")

# Remove Cvi
df <- df_cmt2_coord[!(df_cmt2_coord$population=="Cvi"),]

# Order accesstion by longitude
df$population_nb <- factor(df$population_nb, levels = df$population_nb[order(df$long)])

ggplot(data = df, aes(population_nb, freq_derived)) +
  ggtitle("Derived CMT2 allele frequency") +
  geom_bar(aes(x = population_nb, y = freq_derived), stat = "identity", colour = "black", fill = "grey") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Derived allele frequency") +
  xlab("Population") +
  scale_y_continuous(labels = scales::percent)
```

![](images/CMT2_deletion_frequency.png)

## Plot diagram FBX5

``` r
# Create new variable which concatenat population name and number of accessions
df_fbx5_coord$population_nb <- paste(df_fbx5_coord$population, " n=", df_fbx5_coord$total, sep="")

# Remove Cvi
df <- df_fbx5_coord[!(df_fbx5_coord$population=="Cvi"),]

# Order accesstion by longitude
df$population_nb <- factor(df$population_nb, levels = df$population_nb[order(df$long)])

ggplot(data = df, aes(population_nb, freq_derived)) +
  ggtitle("Derived ARABIDILLO-1 allele frequency") +
  geom_bar(aes(x = population_nb, y = freq_derived), stat = "identity", colour = "black", fill = "grey") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Derived allele frequency") +
  xlab("Population") +
  scale_y_continuous(labels = scales::percent)
```

![](images/FBX5_deletion_frequency.png)

## Plot diagram all alleles

``` r
############# VIM2

# Create new variable which concatenat population name and number of accessions
df_vim2_coord$population_nb <- paste(df_vim2_coord$population, " n=", df_vim2_coord$total, sep="")

# Remove Cvi
df_vim2 <- df_vim2_coord[!(df_vim2_coord$population=="Cvi"),]

# Order accession by longitude
df_vim2$population <- factor(df_vim2$population, levels = df_vim2$population[order(df_vim2$long)])
df_vim2$population_nb <- factor(df_vim2$population_nb, levels = df_vim2$population_nb[order(df_vim2$long)])

df_vim2$allele <- as.factor("VIM2")

############# CMT2

# Remove Cvi
df_cmt2 <- df_cmt2_coord[!(df_cmt2_coord$population=="Cvi"),]

# Order accession by longitude
df_cmt2$population <- factor(df_cmt2$population, levels = df_cmt2$population[order(df_cmt2$long)])
df_cmt2$population_nb <- factor(df_cmt2$population_nb, levels = df_cmt2$population_nb[order(df_cmt2$long)])

df_cmt2$allele <- as.factor("CMT2")

############ FBX5

# Remove Cvi
df_fbx5 <- df_fbx5_coord[!(df_fbx5_coord$population=="Cvi"),]

# Order accession by longitude
df_fbx5$population <- factor(df_fbx5$population, levels = df_fbx5$population[order(df_fbx5$long)])
df_fbx5$population_nb <- factor(df_fbx5$population_nb, levels = df_fbx5$population_nb[order(df_fbx5$long)])

df_fbx5$allele <- as.factor("FBX5")

df_all_alleles <- rbind(df_vim2, df_cmt2, df_fbx5)

df_cmt2_coord$population_nb <- paste(df_cmt2_coord$population, " n=", df_cmt2_coord$total, sep="")

# Remove Cvi
df <- df_cmt2_coord[!(df_cmt2_coord$population=="Cvi"),]

# Order accesstion by longitude
df$population_nb <- factor(df$population_nb, levels = df$population_nb[order(df$long)], ordered=TRUE)

########## Plot

# Colorblind-friendly palette
col <- c("#7570b3","#1b9e77","#d95f02")

ggplot(data = df_all_alleles, aes(population_nb, freq_derived)) +
  ggtitle("Derived allele frequencies per population") +
  geom_bar(aes(x = population_nb, y = freq_derived, fill=allele), position="dodge", stat = "identity") + 
  theme_bw() + scale_fill_manual(values=col) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Derived allele frequency") +
  xlab("Population") +
  scale_y_continuous(labels = scales::percent) +  theme(axis.text.x = element_text(color="black"), 
      axis.text.y = element_text(color="black"),
      axis.ticks = element_line(color = "black"))
```

![](images/deletion_frequency_all_alleles.png)

# Plot gbM by VIM2/4 allele

``` r
df_accessions <- read.table("data/df_accessions_83.txt", header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings="")  

path_DB <- "F:/NETSCRATCH/methylKit_DB_files/1001_project"

df_name <- "df_mean_genes"
title <- "Weighted Methylation Level for genes"

#get_df_wml(list_methylRawLists_genes, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Merge dataframes of methylation levels and df_accessions to get detailed information about dataset
df_mean <- merge(get(df_name), df_accessions, by="sample")

df_mean$FBX5 <- as.factor(df_mean$FBX5)
df_mean$CMT2 <- as.factor(df_mean$CMT2)
df_mean$VIM2 <- as.factor(df_mean$VIM2)

ggplot_per_context <- function(df, context, group){
  
  give.n <- function(x){
    return(c(y = mean(x), label = length(x)))
  }
  
  df_context <- df[(df$context == context),]
  
  ggplot(data=df_context, aes_string(x=group, y="percent_methylation", group=group))+
    ggtitle(paste(context, " methylation", sep="")) + theme_bw()+
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(data=df_context[-82, ], height=.05, width=.05, size=0.8) +
    geom_point(data=df_context[82,],size=2, fill="red", shape=23, width=0.1, height=0) +
    stat_summary(fun.data = give.n, geom = "text", position = position_stack(vjust = 0.7)) +
    scale_y_continuous(breaks=c(seq(0,15,2)))
}

ggplot_per_context(df_mean, "CpG", "VIM2") + xlab("") + 
  ylab("% of methylated cytosines (CG)") + ggtitle("Gene body methylation")
```

![](images/gbM_VIM2ref_VIM2del_SA.png)

#### Statistical test VIM2ref vs VIM2del

``` r
library(onewaytests)

# Test if variances equals
bartlett.test(percent_methylation~VIM2, data=df_mean[(df_mean$context=="CpG"),])
# Variances not equal => use Welch's test

welch.test(percent_methylation~VIM2, data = df_mean[(df_mean$context=="CpG"),], alpha=0.05)

with(df_mean[(df_mean$context=="CpG"),], t.test(percent_methylation[VIM2==0], percent_methylation[VIM2==1], var.equal=FALSE))

# Either 1 or 2-sided, p-value is below 2.2e-16
```

    Bartlett test of homogeneity of variances

    data:  percent_methylation by VIM2
    Bartlett's K-squared = 21.233, df = 1, p-value = 4.066e-06


      Welch's Heteroscedastic F Test (alpha = 0.05) 
    ------------------------------------------------------------- 
      data : percent_methylation and VIM2 

      statistic  : 257.2669 
      num df     : 1 
      denom df   : 60.19284 
      p.value    : 2.09147e-23 

      Result     : Difference is statistically significant. 
    ------------------------------------------------------------- 

        Welch Two Sample t-test

    data:  percent_methylation[VIM2 == 0] and percent_methylation[VIM2 == 1]
    t = 16.04, df = 60.193, p-value < 2.2e-16
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     2.580910 3.316306
    sample estimates:
    mean of x mean of y 
    10.528693  7.580085 

# Plot TE methylation by FBX5 allele

## mCG in long TEs

``` r
df_accessions <- read.table("data/df_accessions_83.txt", header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings="")  

path_DB <- "F:/NETSCRATCH/methylKit_DB_files/1001_project"

df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for genes"

#get_df_wml(list_methylRawLists_genes, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Merge dataframes of methylation levels and df_accessions to get detailed information about dataset
df_mean <- merge(get(df_name), df_accessions, by="sample")

df_mean$FBX5 <- as.factor(df_mean$FBX5)
df_mean$CMT2 <- as.factor(df_mean$CMT2)
df_mean$VIM2 <- as.factor(df_mean$VIM2)


ggplot_per_context <- function(df, context, group){
  
  give.n <- function(x){
    return(c(y = mean(x), label = length(x)))
  }
  
  df_context <- df[(df$context == context),]

  ggplot(data=df_context, aes_string(x=group, y="percent_methylation", group=group))+
    ggtitle(paste(context, " methylation", sep="")) + theme_bw()+
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(data=df_context[-82, ], height=.05, width=.05, size=0.8) +
    geom_point(data=df_context[82,],size=2, fill="red", shape=23, width=0.1, height=0) +
   stat_summary(fun.data = give.n, geom = "text", position = position_stack(vjust = 0.7)) +
    coord_cartesian(ylim=c(72,84)) +
    scale_y_continuous(breaks=c(seq(72,94,2)))
}

ggplot_per_context(df_mean, "CpG", "FBX5") + xlab("") + 
  ylab("% of methylated cytosines (CG)") + ggtitle("mCG at long TEs")
```

![](images/mCG_long_TEs_FBX5.png) \### Statistics

``` r
library(onewaytests)

# Test if variances equals
bartlett.test(percent_methylation~FBX5, data=df_mean[(df_mean$context=="CpG"),])
# Variances not equal => use Welch's test

welch.test(percent_methylation~FBX5, data = df_mean[(df_mean$context=="CpG"),], alpha=0.05)

#with(df_mean[(df_mean$context=="CpG"),], t.test(percent_methylation[FBX5==0], percent_methylation[FBX5==1], var.equal=FALSE))
```


        Bartlett test of homogeneity of variances

    data:  percent_methylation by FBX5
    Bartlett's K-squared = 4.812, df = 1, p-value = 0.02826

      Welch's Heteroscedastic F Test (alpha = 0.05) 
    ------------------------------------------------------------- 
      data : percent_methylation and FBX5 

      statistic  : 172.1605 
      num df     : 1 
      denom df   : 80.77658 
      p.value    : 1.015067e-21 

      Result     : Difference is statistically significant. 
    ------------------------------------------------------------- 

## mCHG in long TEs

``` r
ggplot_per_context <- function(df, context, group){
  
  give.n <- function(x){
    return(c(y = mean(x), label = length(x)))
  }
  
  df_context <- df[(df$context == context),]

  ggplot(data=df_context, aes_string(x=group, y="percent_methylation", group=group))+
    ggtitle(paste(context, " methylation", sep="")) + theme_bw()+
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(data=df_context[-82, ], height=.05, width=.05, size=0.8) +
    geom_point(data=df_context[82,],size=2, fill="red", shape=23, width=0.1, height=0) +
   stat_summary(fun.data = give.n, geom = "text", position = position_stack(vjust = 0.7)) +
    coord_cartesian(ylim=c(38,56)) +
    scale_y_continuous(breaks=c(seq(38,56,2)))
}

ggplot_per_context(df_mean, "CHG", "FBX5") + xlab("") + 
  ylab("% of methylated cytosines (CHG)") + ggtitle("mCHG at long TEs")
```

![](images/mCHG_long_TEs_FBX5.png)

### Statistics

``` r
library(onewaytests)

# Test if variances equals
bartlett.test(percent_methylation~FBX5, data=df_mean[(df_mean$context=="CHG"),])
# Variances equal => use T test

#welch.test(percent_methylation~FBX5, data = df_mean[(df_mean$context=="CpG"),], alpha=0.05)

with(df_mean[(df_mean$context=="CHG"),], t.test(percent_methylation[FBX5==0], percent_methylation[FBX5==1], var.equal=TRUE))
```

        Bartlett test of homogeneity of variances

    data:  percent_methylation by FBX5
    Bartlett's K-squared = 1.4562, df = 1, p-value = 0.2275

        Two Sample t-test

    data:  percent_methylation[FBX5 == 0] and percent_methylation[FBX5 == 1]
    t = -7.5565, df = 81, p-value = 5.582e-11
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -4.918492 -2.868187
    sample estimates:
    mean of x mean of y 
     45.36060  49.25394 

## mCHH in long TEs

``` r
ggplot_per_context <- function(df, context, group){
  
  give.n <- function(x){
    return(c(y = mean(x), label = length(x)))
  }
  
  df_context <- df[(df$context == context),]

  ggplot(data=df_context, aes_string(x=group, y="percent_methylation", color="CMT2"))+
    theme_bw()+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data=df[-82,], position=position_jitterdodge(), size=1) +
    geom_jitter(data=df[82, ],size=3, fill="white", shape=23, width=0.1, height=0) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_y_continuous(limits=c(0,14), breaks=seq(0,14,2))
  
}

ggplot_per_context(df_mean, "CHH", "FBX5") + xlab("FBX5 allele") + ylab("% of methylated cytosines (CHH)")
```

![](images/mCHH_long_TEs_FBX5_CMT2.png) We have 2% of mCHH in FBX5stop
vs 1.2% in FBX5ref.

``` r
df_mean %>% filter(context=="CHH") %>% group_by(FBX5, CMT2) %>% summarize(meth=mean(percent_methylation), sd=sd(percent_methylation)) %>% knitr::kable()
```

| FBX5 | CMT2 |      meth |        sd |
|:-----|:-----|----------:|----------:|
| 0    | 0    | 10.618035 | 1.3241299 |
| 0    | 1    |  2.562190 | 0.4001650 |
| 1    | 0    | 11.741630 | 0.9626544 |
| 1    | 1    |  2.865667 | 0.8650467 |

### Statistics

### Test difference mCHH for FBX5 in CMT2ref/stop background

Note that all these stats are probably better explained by a linear
model integrating the different genes as factor (see part linear model)

Test difference in CMT2ref background

``` r
require(onewaytests)

df_mean_formated_CMT2stop <- df_mean_formated %>% filter(df_mean_formated$CMT2==0) %>% filter(context=="CHH")

bartlett.test(percent_methylation~FBX5, data=df_mean_formated_CMT2stop)
# Variances equal

with(df_mean_formated_CMT2stop[(df_mean_formated_CMT2stop$context=="CHH"),], t.test(percent_methylation[FBX5==0], percent_methylation[FBX5==1], var.equal=TRUE))
```

    data:  percent_methylation by FBX5
    Bartlett's K-squared = 3.0462, df = 1, p-value = 0.08093


        Two Sample t-test

    data:  percent_methylation[FBX5 == 0] and percent_methylation[FBX5 == 1]
    t = -3.1466, df = 54, p-value = 0.002688
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.4622313 -0.1024430
    sample estimates:
    mean of x mean of y 
     1.919552  2.201889 

Test difference in CMT2stop background

``` r
df_mean_formated_CMT2ref <- df_mean_formated %>% filter(df_mean_formated$CMT2==1) %>% filter(context=="CHH")

bartlett.test(percent_methylation~FBX5, data=df_mean_formated_CMT2ref)
# Variances equal

with(df_mean_formated_CMT2ref[(df_mean_formated_CMT2ref$context=="CHH"),], t.test(percent_methylation[FBX5==0], percent_methylation[FBX5==1], var.equal=TRUE))
```

        Bartlett test of homogeneity of variances

    data:  percent_methylation by FBX5
    Bartlett's K-squared = 1.3303, df = 1, p-value = 0.2487


        Two Sample t-test

    data:  percent_methylation[FBX5 == 0] and percent_methylation[FBX5 == 1]
    t = -0.95461, df = 25, p-value = 0.3489
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.5672242  0.2079338
    sample estimates:
    mean of x mean of y 
     1.166271  1.345917 

# GWAS gbM with VIM2/4 insertion as covariate

Perform GWAS using a new VCF file which contains the VIM2 deletion.

I need to add manually the status of the deletions for the 83 samples in
the VCF file
`subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf.gz`

``` bash
cd /srv/netscratch/irg/grp_hancock/johan/GWAS/dna_methylation/GWAS_83_SA

# Copy the original VCF file
cp subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf.gz subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons_plus_VIM_del.recode.vcf.gz

# Uncompress the file to be able to add the deletion
gunzip subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons_plus_VIM_del.recode.vcf.gz

# Keep the 7th columns, which contains VIM2 status (the df_accessions_83.txt is not ordered as in the VCF file)
# Note I use sort as that is how the accessions seq_ID is sorted in the VCF file
cut -f2,7 df_accessions_83.txt | sort | grep -v "seq_ID" - | cut -f2  > accessions_83_order_genotype_VIM2_del.txt

# Add fake GQ and DP value (let's set to 30 and 5 => for instance 1:30:5)
sed 's/$/:30:5/' accessions_83_order_genotype_VIM2_del.txt > accessions_83_order_genotype_VIM2_del_GQ30_DP5.txt

# Row to column conversion
tr -s '\n'  '\t' < accessions_83_order_genotype_VIM2_del_GQ30_DP5.txt > accessions_83_order_genotype_VIM2_del_GQ30_DP5_column.txt

# The deletion is located in Chr1 24,586,731..24,589,471
# Looking manually, the line should be added after line 213419 (chr1:24586704) and before the line 213420  (chr1 24586759)

# I can invent the SNP (for instance T>C). Keep, NS, DP, AN, AC as for Chr1:24589515 => it become lines 213420
# Add at the beginning of accessions_83_order_genotype_VIM2_del_GQ30_DP5_column.txt this (check tab and no space)
Chr1    24586731    .   T   C   40  q25 NS=1;DP=90360;AN=84;AC=1    GT:GQ:DP

# Open file and add the first part of the VCF line (check no spaces but tabs)

# Insert at line 213419 => the new line will become 73728
sed -i -e "213419r accessions_83_order_genotype_VIM2_del_GQ30_DP5_column.txt" subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons_plus_VIM_del.recode.vcf


# Check the new VCF file manually to see if no spaces instead of tabs. Correct insertion at line 213420
# Remove manually a tab at the end of the line

# Compress and tabix
bgzip subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons_plus_VIM_del.recode.vcf && tabix subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons_plus_VIM_del.recode.vcf.gz

# Working directory
cd /srv/netscratch/irg/grp_hancock/johan/GWAS/dna_methylation/GWAS_83_SA

# No need to create a new output directory as the output file with have a "covariate_VIM2" suffix.

# Create file withs first column of 1s, second column with the VIM2 deletion genotype
awk -v OFS="\t" '{print "1",$1}' accessions_83_order_genotype_VIM2_del.txt > covariate_VIM2.txt

VCF="subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons_plus_VIM_del.recode.vcf.gz"

# Run GWAS for gbM
bash ../scripts/run_gwas_gemma.sh CpG_genes.tsv vcf_files/$VCF covariate_VIM2.txt
```

``` r
source("scripts/functions_gwas.R")

dir_file <- "data/output/"

file.name <- "CpG_genes_covariate_VIM2.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

# Check code https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
#source("manhattan_plot_fun.R")

# Remove dots below 0.5 to reduce size of PDF and lag in Inkscape
#gwas.results <- gwas.results %>% filter(-log(P) > 0.5)

# Won't work if no SNPs within the gene region, make a larger range to get some SNPs to be displayed
VIM2 <- c(24589343,24592780)

# Add n bp before and after the gene coordinates to widen the regions to find. I played by adding 0 until I could see the 4 genes on the plot
n = 50000
VIM2 <- c(VIM2[[1]]-n, VIM2[[2]]+n)

ann<-rep(1, length(gwas.results$P))
ann[with(gwas.results, CHR==1 & BP>=VIM2[[1]] & BP<VIM2[[2]])]<-2
ann<-factor(ann, levels=1:2, labels=c("","VIM2"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title=file.name)
```

![](images/GWAS_mCG_genes_VIM2del_covariate.png)

# GWAS mCG genome-wide with VIM2del as covariate

``` bash
VCF="subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons_plus_VIM_del.recode.vcf.gz"

# Run GWAS for gbM
bash ../scripts/run_gwas_gemma.sh CpG_whole_genome.tsv vcf_files/$VCF covariate_VIM2.txt
```

``` r
source("scripts/functions_gwas.R")

dir_file <- "data/output/"

file.name <- "CpG_whole_genome_covariate_VIM2.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

# Check code https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
#source("manhattan_plot_fun.R")

# Remove dots below 0.5 to reduce size of PDF and lag in Inkscape
#gwas.results <- gwas.results %>% filter(-log(P) > 0.5)

# Won't work if no SNPs within the gene region, make a larger range to get some SNPs to be displayed
VIM2 <- c(24589343,24592780)

# Add n bp before and after the gene coordinates to widen the regions to find. I played by adding 0 until I could see the 4 genes on the plot
n = 50000
VIM2 <- c(VIM2[[1]]-n, VIM2[[2]]+n)

ann<-rep(1, length(gwas.results$P))
ann[with(gwas.results, CHR==1 & BP>=VIM2[[1]] & BP<VIM2[[2]])]<-2
ann<-factor(ann, levels=1:2, labels=c("","VIM2"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title=file.name)
```

![](images/GWAS_mCG_whole_genome_VIM2_covariate.png)

# GWAS mCHH long TEs with CMT2 as covariate

``` bash

# Note I use sort as that is how the accessions seq_ID is sorted in the VCF file
cut -f2,6 df_accessions_83.txt | sort | grep -v "seq_ID" - | cut -f2  > accessions_83_order_genotype_CMT2.txt

# Create file withs first column of 1s, second column with the CMT2 genotype
awk -v OFS="\t" '{print "1",$1}' accessions_83_order_genotype_CMT2.txt > covariate_CMT2.txt

VCF="subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf.gz"

bash ../scripts/run_gwas_gemma.sh CHH_TEs_4kb.tsv vcf_files/$VCF covariate_CMT2.txt
```

``` r
source("scripts/functions_gwas.R")

dir_file <- "data/output/"

file.name <- "CHH_TEs_4kb_covariate_CMT2.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

ann<-rep(1, length(gwas.results$P))
#ann<-factor(ann, levels=1:2, labels=c("","VIM2"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title="GWAS for mCHH at long TEs (>4kb) with CMT2stop SNP as covariate")
```

![](images/GWAS_mCHH_long_TEs_CMT2_covariate.png)

# GWAS mCHH genome-wide CMT2 as covariate

``` bash

# Note I use sort as that is how the accessions seq_ID is sorted in the VCF file
cut -f2,6 df_accessions_83.txt | sort | grep -v "seq_ID" - | cut -f2  > accessions_83_order_genotype_CMT2.txt

# Create file withs first column of 1s, second column with the CMT2 genotype
awk -v OFS="\t" '{print "1",$1}' accessions_83_order_genotype_CMT2.txt > covariate_CMT2.txt


VCF="subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf.gz"

# Whole genome
bash ../scripts/run_gwas_gemma.sh CHH_whole_genome.tsv vcf_files/$VCF covariate_CMT2.txt
```

``` r
source("scripts/functions_gwas.R")

dir_file <- "data/output/"

file.name <- "CHH_whole_genome_covariate_CMT2.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

ann<-rep(1, length(gwas.results$P))
#ann<-factor(ann, levels=1:2, labels=c("","VIM2"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title="GWAS for mCHH genome-wide with CMT2stop SNP as covariate")
```

![](images/GWAS_mCHH_whole_genome_CMT2_covariate.png)

# GWAS mCHG genome-wide CMT2 as covariate

``` bash
cd data

# Note I use sort as that is how the accessions seq_ID is sorted in the VCF file
cut -f2,6 df_accessions_83.txt | sort | grep -v "seq_ID" - | cut -f2  > accessions_83_order_genotype_CMT2.txt

# Create file withs first column of 1s, second column with the CMT2 genotype
awk -v OFS="\t" '{print "1",$1}' accessions_83_order_genotype_CMT2.txt > covariate_CMT2.txt


VCF="subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf.gz"

bash ../scripts/run_gwas_gemma.sh CHG_whole_genome.tsv vcf_files/$VCF covariate_CMT2.txt
```

``` r
source("scripts/functions_gwas.R")

dir_file <- "data/output/"

file.name <- "CHG_whole_genome_covariate_CMT2.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold

ann<-rep(1, length(gwas.results$P))
#ann<-factor(ann, levels=1:2, labels=c("","VIM2"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title="GWAS for mCHG genome-wide with CMT2stop SNP as covariate")
```

![](images/GWAS_mCHG_whole_genome_CMT2_covariate.png)

# GWAS mCG long TEs with FBX5 as covariate

``` bash
cd data

# Sort by seqID to be in the order used in the VCF file
sort -k2 df_accessions_83.txt > df_accessions_83_seqID_sorted.txt 

# Remove header
sed -i '/name/d' df_accessions_83_seqID_sorted.txt

# No need to create a new output directory as the output file with have a "covariate_VIM2" suffix.
# Create file withs first column of 1s, second column with the ARA1 allele. Get 5th column (ARA1)
awk -v OFS="\t" '{print "1",$5}' df_accessions_83_seqID_sorted.txt > covariate_FBX5.txt

VCF="subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons.recode.vcf.gz"

bash ${path_run_gwas}/run_gwas_gemma.sh CHH_TEs_4kb.tsv $VCF covariate_FBX5.txt
```

``` r
source("scripts/functions_gwas.R")

dir_file <- "data/output/"

file.name <- "CpG_TEs_4kb_covariate_FBX5.assoc.clean.txt"

path.file <- paste(dir_file, file.name, sep="")

gwas.results <- read.delim(path.file, sep="\t")

nb_snps <- dim(gwas.results)[[1]]

## Calculate Bonferroni corrected P-value threshold
bonferroni_threshold <- 0.05/nb_snps

threshold_pvalue <- bonferroni_threshold


ann<-rep(1, length(gwas.results$P))
#ann<-factor(ann, levels=1:2, labels=c("","VIM2"))

manhattan.plot(chr = gwas.results$CHR, pos=gwas.results$BP, pvalue=gwas.results$P, annotate=ann, sig.level=bonferroni_threshold, title="GWAS for mCG at long TEs (>4kb) with FBX5 SNP as covariate")
```

![](images/GWAS_mCG_long_TEs_FBX5_covariate.png)

# kmersGWAS

kmersGWAS software GitHub repository:
<https://github.com/voichek/kmersGWAS> (version 2.0 used
<https://github.com/voichek/kmersGWAS/releases/tag/v0.2-beta>). Download
zip file `v0.2-beta`.

``` bash
mkdir ~/kmersGWAS

mv v0.2-beta.zip ~/kmersGWAS

cd ~/kmersGWAS && unzip v0.2-beta

make
```

The fastq files for all 83 accessions used for GWAS are needed (see
Supplementary Table 14 for ENA accession numbers).

## Running KMC

`scripts/kmersGWAS/kmers_methylome.SA.sh`

``` bash

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
```

Run in LSF array `kmers_methylome.SA.sh` across all 83 accessions:

``` bash
bsub -q normal -R "rusage[mem=15000]" \
  -M 17000 -n 4 -J "kmc[1-83]%10" < kmerGWAS_methylome_SA.sh
```

## Create k-mers list to be used in GWAS from all individuals

``` bash
# Make a list of all output from KMC
DIR='~/output_kmersGWAS/samples'
cat ./samples.txt | awk '{printf "'$DIR'/%s/kmers_with_strand\t%s\n",$1,$1}' > kmers_list_paths.txt
```

## Filter k-mers from separate lists to one list with all k-mers to use

`scripts/kmersGWAS/combine_kmers_methylome_SA.sh`

``` bash

BIN='~/kmersGWAS/bin/'

${BIN}list_kmers_found_in_multiple_samples \
  -l kmers_list_paths.txt \
    -k 31 --mac 5 -p 0.2 -o kmers_to_use

${BIN}build_kmers_table -l kmers_list_paths.txt -k 31 \
  -a kmers_to_use -o kmers_table

${BIN}emma_kinship_kmers -t kmers_table -k 31 --maf 0.05 > kmers_table.kinship
```

``` bash
bsub -q normal -R "rusage[mem=5000]" -M 7000 combine_kmers_methylome_SA.sh
```

## Run the GWAS

``` bash
mkdir output_dir
bsub -q multicore20 -R "rusage[mem=2000]" -M 5000 -n 15 -J "kmers[1-20]%3" < /home/tergemina/CVI_Project/computer/kmersGWAS/runGWASkmers.sh
```

## Prepare output from GWAS for Bowtie2

Convert the *pass_threshold_10per* file of each phenotype into a
multifasta file for bowtie2 alignment.

``` bash
python convertKmers4mapping.py
```

## Map kmers to TAIR10

``` bash
python map_kmers_methylome_SA.py
```

## Manhattan plots

The plotting was performed on a Jupyter notebook.

``` python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
```

``` python
genome=[30427671,19698289,23459830,18585056,26975502]
space=2000000
chr1=30427671 / 2 
chr2=30427671 + 19698289 / 2 + space
chr3=30427671 + 19698289 + 23459830 / 2 + space * 2
chr4=30427671 + 19698289 + 23459830 + 18585056 / 2 + space * 3
chr5=30427671 + 19698289 + 23459830 + 18585056 + 26975502 / 2 + space * 4
```

### Figure 2b

``` python
start = 24580000
end = 24625000

def parse_classical_GWAS(file):
    df = pd.read_table(file)
    df['P'] = -(np.log10(df['p_lrt']))
    df['P_bi'] = (np.log10(df['p_lrt']))
    space=2000000
    df['position'] = np.select([df['chr'] == 1, 
                                                df['chr'] == 2, 
                                                df['chr'] == 3, 
                                                df['chr'] == 4, 
                                                df['chr'] == 5],
                        [df['ps'], 
                        df['ps'] + 30427671 + space, 
                        df['ps'] + 30427671 + 19698289 + space * 2,
                        df['ps'] + 30427671 + 19698289 + 23459830 + space * 3,
                        df['ps'] + 30427671 + 19698289 + 23459830 + 18585056 + space * 4])
    return(df)

def zoom_classical_GWAS(df):    
    return(df.loc[(df['chr'] == 1) &(df['ps'] > start) & (df['ps'] < end)])

def zoom_kmer_GWAS(file):
    df = pd.read_table(file)
    df['Chr_S1-1'] = df['Chr_S1-1'].str.lower()
    df = df[df['MP_TAIR10'] >= 24]
    df['P'] = -(np.log10(df['p_lrt']))
    space=2000000
    df['position'] = np.select([df['Chr_TAIR10'] == 'chr1', df['Chr_TAIR10'] == 'chr2', df['Chr_TAIR10'] == 'chr3', df['Chr_TAIR10'] == 'chr4', df['Chr_TAIR10'] == 'chr5'],
                        [df['Pos_TAIR10'], 
                        df['Pos_TAIR10'] + 30427671 + space, 
                        df['Pos_TAIR10'] + 30427671 + 19698289 + space * 2,
                        df['Pos_TAIR10'] + 30427671 + 19698289 + 23459830 + space * 3,
                        df['Pos_TAIR10'] + 30427671 + 19698289 + 23459830 + 18585056 + space * 4])
    return(df.loc[(df['Chr_TAIR10'] == 'chr1') & (df['Pos_TAIR10'] > start) & (df['Pos_TAIR10'] < end)])

gw_snp=parse_classical_GWAS('CpG_genes.assoc.txt')
snp=zoom_classical_GWAS(gw_snp)
kmer=zoom_kmer_GWAS('CpG_genes.kmers.assoc.txt')

kmer['CVI-0_mapped']=np.select([kmer['Chr_Cvi-0'] == '*'],['No'],'Yes')


def major_formatter(x, pos):
    y = x /1000000
    return "%.2f" % y
cm = 1/2.54
fig = plt.figure(figsize=(7*cm,4*cm))
ax=sns.scatterplot(x='Pos_TAIR10',
                y='P',
                data=kmer,
                color='0',
                marker='+',
                hue='CVI-0_mapped',
                hue_order=['Yes','No'],
                palette=["gold", "red"],
                s=20,
                legend=False)
ax=sns.scatterplot(x='ps',
                y='P',
                data=snp,
                palette='cividis',
                marker='o',
                color='0',
                s=20,
                linewidth=0)
plt.hlines(y=-np.log10(0.05/len(gw_snp)),
             xmin=min(gw_snp['position']),
             xmax=max(gw_snp['position']),
             color='0',
             linestyle='--',
             linewidth=1,
             zorder=1)
ax.annotate(r'$\itVIM4$',xy=(24586760, 21),xytext=(24586760-8000,21.5))
ax.annotate(r'$\itVIM2$',xy=(24589343, 21),xytext=(24589343,21.5))
ax.xaxis.set_major_formatter(major_formatter)
ax.xaxis.set_major_locator(ticker.MultipleLocator(10000))
plt.xlabel('Position in Mb',fontsize=10)
plt.ylabel('-log$_{10}$($\itP$)',fontsize=10)
arrow1 = mpatches.FancyArrowPatch((24586760, 20), (24583740, 20),mutation_scale=10,color='0')
ax.add_patch(arrow1)
arrow2 = mpatches.FancyArrowPatch((24589343, 20), (24592780, 20),mutation_scale=10,color='0')
ax.add_patch(arrow2)
plt.xlim(start,end+5000)
plt.ylim(0,21)
sns.despine(trim=True)
plt.tick_params(labelsize=8)
plt.savefig("Fig2b_20260809.pdf", format="pdf",bbox_inches="tight")
```

![](images/output_3_0.png)

### Figure S1a

``` python
def parse_classical_GWAS(file):
    df = pd.read_table(file)
    df['P'] = -(np.log10(df['p_lrt']))
    df['P_bi'] = (np.log10(df['p_lrt']))
    space=2000000
    df['position'] = np.select([df['chr'] == 1, 
                                                df['chr'] == 2, 
                                                df['chr'] == 3, 
                                                df['chr'] == 4, 
                                                df['chr'] == 5],
                        [df['ps'], 
                        df['ps'] + 30427671 + space, 
                        df['ps'] + 30427671 + 19698289 + space * 2,
                        df['ps'] + 30427671 + 19698289 + 23459830 + space * 3,
                        df['ps'] + 30427671 + 19698289 + 23459830 + 18585056 + space * 4])
    
    return(df)

def parse_kmer_GWAS(file):
    df = pd.read_table(file)
    df['Chr_S1-1'] = df['Chr_S1-1'].str.lower()
    df = df[df['MP_TAIR10'] >= 24]
    df['P'] = -(np.log10(df['p_lrt']))
    space=2000000
    df['position'] = np.select([df['Chr_TAIR10'] == 'chr1', df['Chr_TAIR10'] == 'chr2', df['Chr_TAIR10'] == 'chr3', df['Chr_TAIR10'] == 'chr4', df['Chr_TAIR10'] == 'chr5'],
                        [df['Pos_TAIR10'], 
                        df['Pos_TAIR10'] + 30427671 + space, 
                        df['Pos_TAIR10'] + 30427671 + 19698289 + space * 2,
                        df['Pos_TAIR10'] + 30427671 + 19698289 + 23459830 + space * 3,
                        df['Pos_TAIR10'] + 30427671 + 19698289 + 23459830 + 18585056 + space * 4])
#     return(df)
    return(df)

snp=parse_classical_GWAS('CpG_genes.assoc.txt')
kmer=parse_kmer_GWAS('CpG_genes.kmers.assoc.txt')
top_SNP=snp.loc[(snp['chr'] == 1) & (snp['ps'] == 24622585)]

        
def major_formatter(x, pos):
    y = x /1000000
    return "%.1f" % y
cm = 1/2.54
fig = plt.figure(figsize=(14*cm,4*cm))
ax = sns.scatterplot(x='position',
                y='P',
                data=snp,
                hue='chr',
                palette=['dimgray','silver','dimgray','silver','dimgray'],
                marker='o',
                s=12,
                linewidth=0,
                legend=False)
ax=sns.scatterplot(x='position',
                y='P',
                data=kmer,
                color='0',
                marker='+',
                s=20,
                legend=False)
ax = sns.scatterplot(x='ps',
                y='P',
                data=top_SNP,
                color='purple',
                marker='D',
                s=20,
                linewidth=0)
ax.hlines(y=-np.log10(0.05/len(snp)),
             xmin=min(snp['position']),
             xmax=max(snp['position']),
             color='0',
             linestyle='--',
             linewidth=1,
#              alpha=0.4,
             zorder=1)
ax.vlines(24583740,
          0,
          20,
          colors='0',
          zorder=1,
          linestyle='--',
          linewidth=1)
ax.annotate(r'$\itVIM2/4$',xy=(24586760, 20),xytext=(24586760-9500000,20.5))
ax.annotate(r'$\itVIM3$',xy=(15837178 + 30427671 + 19698289 + 23459830 + 18585056 + space * 4, 20),xytext=(15837178 + + 30427671 + 19698289 + 23459830 + 18585056 + space * 4 -5000000,20.5))
ax.vlines(15837178 + 30427671 + 19698289 + 23459830 + 18585056 + space * 4,
          0,
          20,
          colors='0',
          zorder=1,
          linestyle='--',
          linewidth=1)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.xlabel('Chromosome',fontsize=10)
plt.ylim(0,20)
plt.ylabel('-log$_{10}$($\itP$)',fontsize=10)
plt.xticks([chr1,chr2,chr3,chr4,chr5],["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"])
sns.despine(trim=True)
plt.tick_params(labelsize=8)
plt.savefig("FigSup1a_20260809.pdf", format="pdf",bbox_inches="tight")
```

![](images/output_5_0.png)

### Figure S1b

``` python
start = 14073198
end = 16043404

def parse_classical_GWAS(file):
    df = pd.read_table(file)
    df['P'] = -(np.log10(df['p_lrt']))
    df['P_bi'] = (np.log10(df['p_lrt']))
    space=2000000
    df['position'] = np.select([df['chr'] == 1, 
                                                df['chr'] == 2, 
                                                df['chr'] == 3, 
                                                df['chr'] == 4, 
                                                df['chr'] == 5],
                        [df['ps'], 
                        df['ps'] + 30427671 + space, 
                        df['ps'] + 30427671 + 19698289 + space * 2,
                        df['ps'] + 30427671 + 19698289 + 23459830 + space * 3,
                        df['ps'] + 30427671 + 19698289 + 23459830 + 18585056 + space * 4])
    return(df)

def zoom_classical_GWAS(df):    
    return(df.loc[(df['chr'] == 5) &(df['ps'] > start) & (df['ps'] < end)])

def zoom_kmer_GWAS(file):
    df = pd.read_table(file)
    df['Chr_S1-1'] = df['Chr_S1-1'].str.lower()
    df = df[df['MP_TAIR10'] >= 24]
    df['P'] = -(np.log10(df['p_lrt']))
    space=2000000
    df['position'] = np.select([df['Chr_TAIR10'] == 'chr1', df['Chr_TAIR10'] == 'chr2', df['Chr_TAIR10'] == 'chr3', df['Chr_TAIR10'] == 'chr4', df['Chr_TAIR10'] == 'chr5'],
                        [df['Pos_TAIR10'], 
                        df['Pos_TAIR10'] + 30427671 + space, 
                        df['Pos_TAIR10'] + 30427671 + 19698289 + space * 2,
                        df['Pos_TAIR10'] + 30427671 + 19698289 + 23459830 + space * 3,
                        df['Pos_TAIR10'] + 30427671 + 19698289 + 23459830 + 18585056 + space * 4])
    return(df.loc[(df['Chr_TAIR10'] == 'chr5') & (df['Pos_TAIR10'] > start) & (df['Pos_TAIR10'] < end)])

gw_snp=parse_classical_GWAS('CpG_genes.assoc.txt')
snp=zoom_classical_GWAS(gw_snp)
kmer=zoom_kmer_GWAS('CpG_genes.kmers.assoc.txt')
Chr5_15047549_ld=pd.read_table('Chr5_15047549.ld', sep = '\s+')

zoom=pd.merge(snp,Chr5_15047549_ld,left_on='ps',right_on='BP_B')
zoom.head()

        
def major_formatter(x, pos):
    y = x /1000000
    return "%.1f" % y
cm = 1/2.54
fig = plt.figure(figsize=(15*cm,10*cm))
ax=sns.scatterplot(x='Pos_TAIR10',
                y='P',
                data=kmer,
                color='0',
                marker='+',
                s=20,
                legend=False)
ax=sns.scatterplot(x='ps',
                y='P',
                data=zoom,
                palette='cividis',
                hue='R2',
                marker='o',
                s=20,
                linewidth=0)
plt.hlines(y=-np.log10(0.05/len(gw_snp)),
             xmin=min(gw_snp['position']),
             xmax=max(gw_snp['position']),
             color='0',
             linestyle='--',
             linewidth=1,
             zorder=1)
ax.set_xlim(start,end)
ax.xaxis.set_major_formatter(major_formatter)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.xlabel('Position in Mb',fontsize=10)
plt.ylabel('-log$_{10}$($\itP$)',fontsize=10)
ax.annotate(r'VIM3',xy=(15837178, 19),xytext=(15837178-300000,19),arrowprops=dict(facecolor='black',arrowstyle="-|>"),zorder=7)
ax.vlines(15837178,0,19,colors='r',zorder=0)
sns.despine(trim=True)
plt.tick_params(labelsize=8)
norm = plt.Normalize(zoom['R2'].min(), zoom['R2'].max())
sm = plt.cm.ScalarMappable(norm=norm,cmap='cividis')
ax.get_legend().remove()
cbar=ax.figure.colorbar(sm,pad=0.05,shrink=0.5)
cbar.ax.set_title("LD (D')", fontsize=10, pad=10, x=1.2)
plt.savefig("FigSup1b_20260809.pdf", format="pdf",bbox_inches="tight")
```

![](images/output_7_0.png)

## Kmer mapping in the four CVI assemblies

### Identify coordinates VIM2/4 and VIM3 in the four assemblies

We first identified the positions of VIM2/4 and VIM3 in the four CVI
assemblies we generated: Cvi-0, S1-1, S5-10, and S15-3

Download the fasta files of these accession on NCBI (BioProject
PRJNA1112558).

For this, we used blastn.

Get TAIR10 assembly from www.arabidopsis.org
(<https://www.arabidopsis.org/download/file?path=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz>)

``` bash
# Rename chromosome
gunzip TAIR10_chr_all.fas.gz

# Replace chromosome names to add prefix Chr
sed 's/^>\([1-5]\)/>Chr\1/g' TAIR10_chr_all.fas > TAIR10.fa
```

Build blastn indexes for the five assemblies:

``` bash
# Make blast DB
for i in *fasta; do 
  makeblastdb -in $i -dbtype nucl -out ${i%.*}
done
```

#### VIM4 coordinates

AT1G66040

Get full gDNA sequence from Col-0 on
<https://www.arabidopsis.org/sequence?key=1002494865> (fasta file in
`data/fasta_files/VIM4_gDNA_TAIR10.fa`)

``` bash

blastn -version
blastn: 2.2.27+

# Blast gDNA sequence
for i in *fasta; do
 blastn -query VIM4_gDNA_TAIR10.fa -db ${i%.*} -outfmt 6 > vim4_TAIR10_${i%.*}.blast
done

# Start coordinates
for i in vim4_*_unmasked.blast; do grep "Chr1" $i  | sort -k7 | head -n1 | cut -f2,9; done
Chr1    25129235
Chr1    25173075
Chr1    25149867
Chr1    25131477



# Stop coordinates => tricky because VIM2 5' is very close to VIM4 5'
# Choose value manually based on VIM2 start positions
25129235
25173075
25149867
25131477
```

#### VIM2 coordinates

AT1G66050

Get gDNA sequence from Col-0 on
<https://www.arabidopsis.org/sequence?key=2002965408> (fasta file in
`data/fasta_files/VIM2_gDNA_TAIR10.fa`)

``` bash

# Blast
for i in *fasta; do
 blastn -query VIM2_gDNA_TAIR10.fa -db ${i%.*} -outfmt 6 > vim2_TAIR10_${i%.*}.blast
done

# Start coordinates
for i in vim2_*_unmasked.blast; do grep "Chr1" $i  | sort -k7 | head -n1 | cut -f2,9; done
Chr1    25132496
Chr1    25178990
Chr1    25155780
Chr1    25134739


# Stop coordinates
for i in vim2_*_unmasked.blast; do grep "Chr1" $i  | sort -k10 | tail -n1 | cut -f2,10; done
Chr1    25135833
Chr1    25182474
Chr1    25159265
Chr1    25138076
```

#### VIM3 coordinates

AT5G39550

Get gDNA sequence from Col-0 in
<https://www.arabidopsis.org/sequence?key=2503953419> (fasta file in
`data/fasta_files/VIM3_gDNA_TAIR10.fa`).

``` bash

for i in *fasta; do
 blastn -query VIM3_gDNA_TAIR10.fa -db ${i%.*} -outfmt 6 > vim3_TAIR10_${i%.*}.blast
done

for i in vim3_TAIR10_*.blast; do
echo $i
head -n 1 $i
done

for i in *fasta; do echo $i;  blastn -query $VIM3_TAIR10 -db ${i%.*} -outfmt 6 | head -n 1 | cut -f2,9,10; done

CVI-0_unmasked.fasta
Chr5    15556114        15559612
S1-1_unmasked.fasta
Chr5    15834646        15838144
S15-3_unmasked.fasta
Chr5    15742749        15746247
S5-10_unmasked.fasta
Chr5    15837351        15840849
TAIR10.fasta
chr5    15837177        15840678
```

#### Summary positions VIMs in the five assemblies

VIM4

| assembly | chr  | start      | end        | length |
|----------|------|------------|------------|--------|
| Cvi-0    | chr1 | 25,129,235 | 25,132,504 | 3,269  |
| S1-1     | chr1 | 25,173,075 | 25,176,360 | 3,285  |
| S15-3    | chr1 | 25,149,867 | 25,153,150 | 3,283  |
| S5-10    | chr1 | 25,131,477 | 25,134,747 | 3,270  |
| TAIR10   | chr1 | 24,583,740 | 24,586,760 | 3,020  |

VIM2

| assembly | chr  | start      | end        | length |
|----------|------|------------|------------|--------|
| Cvi-0    | chr1 | 25,132,496 | 25,135,833 | 3,337  |
| S1-1     | chr1 | 25,178,990 | 25,182,474 | 3,484  |
| S15-3    | chr1 | 25,155,780 | 25,159,265 | 3,485  |
| S5-10    | chr1 | 25,134,739 | 25,138,076 | 3,337  |
| TAIR10   | chr1 | 24,589,342 | 24,592,780 | 3,438  |

Note that VIM4 gDNA regions are longer in CVI due to some
TAIR10-specific deletions present in the first intron.

VIM3

| assembly | chr  | start      | end        | length |
|----------|------|------------|------------|--------|
| Cvi-0    | chr5 | 15,556,114 | 15,559,612 | 3,498  |
| S1-1     | chr5 | 15,834,646 | 15,838,144 | 3,498  |
| S15-3    | chr5 | 15,742,749 | 15,746,247 | 3,498  |
| S5-10    | chr5 | 15,837,351 | 15,840,849 | 3,498  |
| TAIR10   | chr5 | 15,837,177 | 15,840,678 | 3,501  |

### Blast kmer on the five assemblies

``` bash

# Build bowtie index for the four CVI assemblies and TAIR10 Col-0
for i in *fasta; do
  bowtie-build -f $i ${i%%.*} &
done

# Map kmer 2571 (most significant mapping on chr5 with one mismatch in TAIR10)
for i in *fasta; do
  bowtie ${i%%.*} -c AAACAGTTATTTTCCATTAAACCACCACTAA --all > kmer_2571_mapping_all_${i%%.*}.sam
done

# reads processed: 1
# reads with at least one reported alignment: 1 (100.00%)
# reads that failed to align: 0 (0.00%)
Reported 1 alignments
# reads processed: 1
# reads with at least one reported alignment: 1 (100.00%)
# reads that failed to align: 0 (0.00%)
Reported 2 alignments
# reads processed: 1
# reads with at least one reported alignment: 1 (100.00%)
# reads that failed to align: 0 (0.00%)
Reported 2 alignments
# reads processed: 1
# reads with at least one reported alignment: 1 (100.00%)
# reads that failed to align: 0 (0.00%)
Reported 1 alignments
# reads processed: 1
# reads with at least one reported alignment: 1 (100.00%)
# reads that failed to align: 0 (0.00%)
Reported 2 alignments

for i in kmer_2571_mapping_all*; do
  echo $i; cat $i | cut -f3,4,8
done

kmer_2571_mapping_all_CVI-0_unmasked.sam
Chr5    15559705        16:G>A
kmer_2571_mapping_all_S1-1_unmasked.sam
Chr1    25176565
Chr5    15838237        16:G>A
kmer_2571_mapping_all_S15-3_unmasked.sam
Chr1    25153355
Chr5    15746340        16:G>A
kmer_2571_mapping_all_S5-10_unmasked.sam
Chr5    15840942        16:G>A
kmer_2571_mapping_all_TAIR10.sam
chr5    15840771        16:G>A
chr1    24586963        5:A>C,16:T>A
```

kmer at VIM2/4 promoter region:

| assembly | VIM4 5’ | VIM2 5’ | distance | kmer location | in between? | distance from VIM4 |
|----|----|----|----|----|----|----|
| Cvi-0 | 25,132,504 | 25,132,496 | -8 |  |  |  |
| S1-1 | 25,176,360 | 25,178,990 | 2,630 | 25,176,565 | Yes | 205 |
| S15-3 | 25,153,150 | 25,155,780 | 2,630 | 25,153,355 | Yes | 205 |
| S5-10 | 25,134,747 | 25,134,739 | -8 |  |  |  |
| TAIR10 | 24,586,760 | 24,589,342 | 2,582 | 24,586,964 |  | 204 |

kmer at VIM3 promoter region:

| assembly | 5’ end     | kmer location | distance |
|----------|------------|---------------|----------|
| Cvi-0    | 15,559,612 | 15,559,705    | 93       |
| S1-1     | 15,838,144 | 15,838,237    | 93       |
| S15-3    | 15,746,247 | 15,746,340    | 93       |
| S5-10    | 15,840,849 | 15,840,942    | 93       |
| TAIR10   | 15,840,678 | 15,840,771    | 93       |

# DMR analysis for VIM2

To call for differential methylated regions, we pooled the WGBS data of
the accessions sharing the same alleles. Our sequencing depth was not
enough to use accessions as replicates as we end up with few cytosines
that have enough coverage in each sample to perform proper statistics.

## Pooling of the data

``` bash
# Get WGBS library ID of the samples that segregate for VIM2/4 deletion
awk '$7==0 {print $3}' df_accessions_83.txt > libraries_VIM2wt.txt
awk '$7==1 {print $3}' df_accessions_83.txt > libraries_VIM2del.txt

cd /srv/netscratch/irg/grp_hancock/johan/pooled_by_cmt2_ara1_vim2

nohup xargs zcat < libraries_VIM2del.txt > merged_libraries_VIM2del.fq &
nohup xargs zcat < libraries_VIM2wt.txt > merged_libraries_VIM2wt.fq &

# Check read numbers
echo $(cat merged_libraries_VIM2del.fq | wc -l)/4 | bc
275834214
# 275 M reads

echo $(cat merged_libraries_VIM2wt.fq | wc -l)/4 | bc
308855669

wget https://raw.githubusercontent.com/mojones/random_scripts/14218de511d24b6450df4dc98ca15752626b6797/sample_fastq.py

# Subsample the VIM2wt fastq to have a similar coverage in each group
python sample_fastq.py -n 275834214 merged_libraries_VIM2wt.fq merged_libraries_VIM2wt_276M.fq
# Bug, just run bismark on the full dataset
```

## Run Bismark

``` bash
cd /srv/netscratch/irg/grp_hancock/johan/pooled_by_cmt2_ara1_vim2

ref="/srv/netscratch/irg/grp_hancock/johan/GC_3542_3599/TAIR10_mapping/TAIR10"
output="/srv/netscratch/irg/grp_hancock/johan/pooled_by_cmt2_ara1_vim2/output_bismark_VIM2"

nohup bash /home/zicola/SCRIPTS/WGBS/WGBS_African_Arabidopsis/curated_scripts/run_bismark.sh -r $ref -1 merged_libraries_VIM2del.fq -o $output &

nohup bash /home/zicola/SCRIPTS/WGBS/WGBS_African_Arabidopsis/curated_scripts/run_bismark.sh -r $ref -1 merged_libraries_VIM2wt.fq -o $output &
```

## Create methylKit objects

``` r
# Paths to database
path_DB_CpG="E:/analysis_VIM2/methylKit_DB_files/methylDB_CpG"
path_DB_CHG="E:/analysis_VIM2/methylKit_DB_files/methylDB_CHG"
path_DB_CHH="E:/analysis_VIM2/methylKit_DB_files/methylDB_CHH"

list_DB_paths <- list(path_DB_CpG, path_DB_CHG, path_DB_CHH)

# Add general path to save dataframe (DNA methylation levels for instance
path_DB <- "E:/analysis_VIM2/methylKit_DB_files"

# Path containing cytosine report and bam files from bismark pipeline
path_bismark_files <- "E:/analysis_VIM2/output_bismark"

df_accessions <- data.frame(library=c("merged_libraries_VIM2del", "merged_libraries_VIM2wt"), 
                            name = c("VIM2del", "VIM2wt") , allele = c("VIM2del", "VIM2wt") )   
df_accessions[] <- lapply(df_accessions, as.factor)

# Order the accession as list.files() list the bismark cytosine report files
# Important otherwise the list_samples vector won't contain samples in order
#df_accessions <- order_df_accessions(df_accessions)

list_samples <- as.list(as.vector(df_accessions$library))

# Get list of treatments and reformat so that the first treatment is 0 (control should be 0 optimally)
# Actually, control should be first, treatment second
list_treatments <- c(1,0)

context <- c("CpG","CHG", "CHH")

path_bed <- "data/bed_files/"

# Path to bed files for region analysis
bed_genes <- paste(path_bed, "Araport11_GFF3_genes_only.bed", sep = "")
```

### Create methylRawListDB objects

``` r
import_bismark_cytosine_report(path_bismark_files, list_DB_paths, list_samples, list_treatments)
```

### Load methylRawListDB objects

``` r
# Load methylRawListDB objects (without filtering)
list_methylRawLists_raw <- load_methylRawListDB(list_DB_paths, type="raw", list_samples, list_treatments)
```

### Filter methylRawList raw

``` r
filter_methylRawList(list_methylRawLists_raw)
```

### Load filtered methylRawListDB objects

``` r
list_methylRawLists <- load_methylRawListDB(list_DB_paths, type="filtered", list_samples, list_treatments)
```

## Gene methylation

``` r
# Create subset
subset_methylObject(list_methylRawLists, list_DB_paths, bed_genes, "genes", "methylRaw")

# Load
list_methylRawLists_genes <- load_methylRawListDB(list_DB_paths, type="genes", list_samples, list_treatments)

df_name <- "df_mean_genes"
title <- "Weighted Methylation Level in genes"

get_df_wml(list_methylRawLists_genes, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
ggplot_all(get(df_name), title = title)
```

![](images/gbM_pooled_samples_VIM2.png)

The decrease of gbM in the deletion background is clear.

## Create methlyBaseDB objects

``` r
list_methylBases <- merged_methylRawList(list_methylRawLists)
```

## Load methylBaseDB objects

``` r
list_methylBases <- load_methylBaseDB(list_DB_paths, list_samples, list_treatments)
```

## Create tiles objects

``` r
# Define window and step sizes
tiles_list_methylBases <- tileMethylCounts_batch(list_methylBases, win_size=300, step_size=100, suffix="", list_dbdir = list_DB_paths)
```

## Load tiles objects

``` r
tiles_list_methylBases <- load_methylBaseDB(list_DB_paths, list_samples, list_treatments, suffix="tiled__300_100")
```

## Create methylDiffDB objects

``` r
list_methylDiffs <- lapply(tiles_list_methylBases, calculateDiffMeth, 
                     overdispersion="MN", 
                     test="fast.fisher",mc.cores=1, 
                     save.db = TRUE)
```

## Load methylDiffDB objects

``` r
# Make lists of objects
list_methylDiffs <- load_methylDiffDB(list_DB_paths, list_samples, list_treatments, suffix="tiled__300_100")


# Put data in RAM to avoid trouble when using getMethylDiff
list_methylDiffs_RAM <- lapply(list_methylDiffs, DB_to_RAM_conversion, type="methylDiff")
```

## Create DMRs

Select only DMRs with more than 25% of difference between the 2
treatments. Use 5% as threshold for the qvalue

``` r
# Get all DMRs
list_DMRs_25p <- lapply(list_methylDiffs_RAM, getMethylDiff, difference=25,
                        qvalue=0.05, type="all", save.db=FALSE)

# Keep only mCG DMRs
list_DMRs_25p_mCG <- list_DMRs_25p[[1]]
```

## Distribution of mCG DMRs

``` r
diffMethPerChr(list_DMRs_25p_mCG, plot=TRUE, qvalue.cutoff=0.05, meth.cutoff=25)

# Get list of DMRs distribution per chromosome.
# I tried to put a name to the each plot but it does not work
diffMethPerChr_data <- diffMethPerChr(list_DMRs_25p_mCG, plot=FALSE, qvalue.cutoff=0.05, meth.cutoff=25)

# Number of hyperDMRs
sum(diffMethPerChr_data$diffMeth.per.chr$number.of.hypermethylated)
# 2155

# Number of hypoDMRs
sum(diffMethPerChr_data$diffMeth.per.chr$number.of.hypomethylated)
# 40319

# Total DMRs
sum(diffMethPerChr_data$diffMeth.per.chr$number.of.hypermethylated)+ sum(diffMethPerChr_data$diffMeth.per.chr$number.of.hypomethylated)
```

## Merge mCG DMRs

Not to loose information about methylation direction, I created a
function that average across overlapping DMRs and return a tibble in bed
format (see merge_DMRs in `scripts/functions_methylkit.R`)

``` r
list_DMRs_merged <- lapply(list_DMRs_25p, merge_DMRs)

# Check number of DEGs
# 12996 DMRs_mCG_25p_VIM2.bed

DMRs_merged_mCG <- list_DMRs_merged[[1]]

DMRs_merged_mCG %>% summarize(hyperDMRs=sum(score>0), hypoDMRs=sum(score<0))

DMRs_merged_mCG %>% group_by(chrom) %>% summarize(hyperDMRs=sum(score>0), hypoDMRs=sum(score<0))

saveRDS(DMRs_merged_mCG, "data/DMRs_merged_mCG_VIM2.Rds")

DMRs_merged_mCG <- readRDS("data/DMRs_merged_mCG_VIM2.Rds")

# Also export as txt file
write.table(DMRs_merged_mCG, "data/DMRs_merged_mCG_VIM2.txt", quote=F, row.names = F, sep = "\t")

sum(DMRs_merged_mCG$score < 0)
sum(DMRs_merged_mCG$score > 0)
```

\[1\] 12243 \[1\] 753

## Overlap with genes

``` r
# Get gene list and coordinate in bed format
coordinate_genes <- valr::read_bed(bed_genes, n_fields = 4)

overlap_mCG_DMR_gene_analysis_VIM2 <- valr::bed_intersect(DMRs_merged_mCG, coordinate_genes, suffix=c("_DMRs","_genes"))

# Number of unique DMRs
length(unique(overlap_mCG_DMR_gene_analysis_VIM2$name_DMRs))
# 11948

length(unique(overlap_mCG_DMR_gene_analysis_VIM2$name_genes))
# 7541

saveRDS(overlap_mCG_DMR_gene_analysis_VIM2, "data/overlap_mCG_DMR_gene_analysis_VIM2.rds")

overlap_mCG_DMR_gene_analysis_VIM2 <- readRDS("data/overlap_mCG_DMR_gene_analysis_VIM2.rds")

# Write list 11948 genes
list_DMR_genes <- unique(overlap_mCG_DMR_gene_analysis_VIM2$name_genes)
write(list_DMR_genes, "data/list_DMR_genes_VIM2.txt")
```

11,948 DMRs overlap 7541 unique genes

``` r
hist(overlap_mCG_DMR_gene_analysis_VIM2$score_DMRs)
```

Most DMRs at genes are hypoDMRs (follow distribution of global mCG DMRs)

#### Protein-coding genes

Overlap with protein-coding genes only

``` r
bed_protein_coding_genes <- valr::read_bed(paste(path_bed, "Araport11_GFF3_protein_coding_genes_only.bed", sep=""), n_fields = 6)

overlap_mCG_DMR_gene_protein_coding_analysis_VIM2 <- valr::bed_intersect(DMRs_merged_mCG, bed_protein_coding_genes, suffix=c("_DMRs","_genes"))

length(unique(overlap_mCG_DMR_gene_protein_coding_analysis_VIM2$name_genes))
# 7200

length(unique(overlap_mCG_DMR_gene_protein_coding_analysis_VIM2$name_DMRs))
# 11841
```

7200 protein-coding genes overlap DMRs out of 27655 protein-coding genes

``` r
list_DMR_PC_genes <- unique(overlap_mCG_DMR_gene_protein_coding_analysis_VIM2$name_genes)
write(list_DMR_PC_genes, "data/list_DMR_PC_genes_VIM2.txt")
```

## GO analysis

The list of the 7200 PC genes in `list_DMR_PC_genes_VIM2.txt` were
copied in the DAVID Functional annotation webtool
(<https://davidbioinformatics.nih.gov/>) (Step 1), using `TAIR_ID` as
identifier (Step 2: Select Identifier), then “Gene List” (Step 3: List
Type). Gene Ontology results were exported and put in Excel format
(Supplementary Table 5).

DAVID output were put in text files to be plotted in R:

``` r
CC <-  read.delim("data/CC_DAVID_DMR_VIM2.txt")

# Arrange by fold enrichment

CC %>% mutate(Name = fct_reorder(Name, Fold.Enrichment)) %>%  ggplot(aes(x=Fold.Enrichment, y=Name, color=-log10(FDR))) +
  geom_point(aes(size = gene))  +  scale_size_continuous(trans="sqrt",range = c(1, 10),breaks=c(100,500, 1000,2000, 3000)) + xlab("Fold enrichment") + ylab("GO term") + ggtitle("GO - Cellular components") + theme_bw()+ theme(axis.text.x = element_text(color="black"),axis.text.y = element_text(color="black"),axis.ticks = element_line(color = "black")) + theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(limits=c(1,2.5), breaks=seq(1,2.5,0.5))
```

![](images/GO_CC_DMRs.png)

``` r
BP <- read.delim("data/BP_DAVID_DMR_VIM2.txt")

BP %>% mutate(Name = fct_reorder(Name, Fold.Enrichment)) %>%  ggplot(aes(x=Fold.Enrichment, y=Name, color=-log10(FDR))) +
  geom_point(aes(size = gene)) + xlab("Fold enrichment") + ylab("GO term") + ggtitle("GO - Biological processes") + theme_bw()+ theme(axis.text.x = element_text(color="black"),axis.text.y = element_text(color="black"),axis.ticks = element_line(color = "black")) + theme(plot.title = element_text(hjust = 0.5))  + scale_size_area(max_size = 10) + scale_x_continuous(limits=c(1.5,3), breaks=seq(1.5,3,0.5))
```

![](images/GO_BP_DMRs.png)

## Overlap of mCG DMRs with promoter regions

Check if promoter regions are also affected in VIM2del and which are the
associated genes. Then define the overlap with DEGs called between
VIM2del and VIM2ref.

``` r
library(rtracklayer)
library(GenomicRanges)

DMRs_merged_mCG_VIM2 <- read.delim("data/DMRs_merged_mCG_VIM2.txt")

DMR_gr=GRanges(seqnames=DMRs_merged_mCG_VIM2$chrom,
           ranges=IRanges(DMRs_merged_mCG_VIM2$start,
                          end=DMRs_merged_mCG_VIM2$end),
                          score=DMRs_merged_mCG_VIM2$score)

path_bed <- "data/bed_files/"

# Path to bed files for region analysis
bed_genes <- paste(path_bed, "Arabidopsis_thaliana.TAIR10.39.bed", sep = "")

# Get gene list and coordinate in bed format
coordinate_genes <- valr::read_bed(bed_genes)

# Convert to GRange
coordinate_genes_gr <- rtracklayer::import(bed_genes)

# Get 500 bp promoter
coordinate_500bp_promoter_gr <- GenomicRanges::flank(coordinate_genes_gr, width=500, start=TRUE)

overlap_mCG_DMR_gene_analysis_VIM2 <- valr::bed_intersect(DMR_gr, coordinate_500bp_promoter_gr, suffix=c("_DMRs","_genes"))

# Find DMR that overlap at least 400
overlap_DMR_promoters <- IRanges::subsetByOverlaps(coordinate_500bp_promoter_gr, DMR_gr, minoverlap=500)

# Check if overlap with VIM2 DEGs

significant_94_DEGs_VIM2 <- read.delim("data/significant_94_DEGs_VIM2.txt")

# Get gene ID without transcript suffix and put all in uppercase
gene_DMR_promoter <- overlap_DMR_promoters$name
gene_DMR_promoter <- gsub("\\..*","", gene_DMR_promoter)
gene_DMR_promoter <- toupper(gene_DMR_promoter)

intersect(gene_DMR_promoter, significant_94_DEGs_VIM2$geneID)
# 3 out of 94
```

- AT2G10608
- AT5G03090
- AT5G01055

AT2G10608 and AT5G03090 are protein-coding genes while AT5G01055 encodes
a long non-coding RNA.

Check all genes with DMRs and the odds to overlap with 3 of the 94 DEGs.

``` r
genes_with_DMRs <- read.delim("data/list_DMR_genes_VIM2.txt", header=F)

# Number of genes with DMRs = 751
# Number of genes with DMRs in their promoter = 101
# Number of DEGs = 94
# Number of genes in the overlap = 3
# Number of genes used for DEG analysis = 37686

# Using number of genes as universe
library(GeneOverlap)
overlap <- newGeneOverlap(gene_DMR_promoter,significant_94_DEGs_VIM2$geneID, genome.size=37686)
testGeneOverlap(overlap )
```

    GeneOverlap object:
    listA size=101
    listB size=94
    Intersection size=3
    Overlapping p-value=2.1e-03
    Jaccard Index=0.0

p-value = 0.2%

``` r
# Using number of DMRs as universe
library(GeneOverlap)
overlap <- newGeneOverlap(gene_DMR_promoter,significant_94_DEGs_VIM2$geneID, genome.size=12996)
testGeneOverlap(overlap )
```

    GeneOverlap object:
    listA size=101
    listB size=94
    Intersection size=3
    Overlapping p-value=0.037
    Jaccard Index=0.0

p-value = 3.7%

The intersect is significant.

# Analysis DNA methylation in vim mutants

We reused the WGBS data generated for WT Col-0 and the *vim* mutants
from Stroud et al 2013 (10.1016/j.cell.2012.10.054). BioProject
PRJNA172021.

SRA files:

- vim1 SRR534263
- vim2 SRR534264
- vim3 SRR534265
- vim1/2/3 SRR534266
- WT rep 2 SRR534177
- WT rep 3 SRR534193

Download fastq files:

``` bash
fastq-dump --split-spot SRR534263
fastq-dump --split-spot SRR534264
fastq-dump --split-spot SRR534265
fastq-dump --split-spot SRR534266
fastq-dump --split-spot SRR534177
fastq-dump --split-spot SRR534193
```

WT rep 1 matches the BioSample SAMN03765231 and is splitted into 96
different fastq files (SE 51 bp). Download the SraRunTable on NCBI (Send
results to Run selector) and retrieve all SRA that match WT rep1.

``` bash
grep "Wild" SraRunTable_PRJNA172021.txt | grep -v "replicate 2" - | grep -v "replicate 3" - | cut -d',' -f1 - > SRR_WT_rep1.txt

while read i; do
    bsub "fastq-dump --split-spot $i"
done < SRR_WT_rep1.txt

# concatenate all files for WT Rep1
cat *fastq.gz > WT_rep1.fastq.gz
```

## Run Bismark

For each fastq file, run the following command:

``` bash
bash run_bismark.sh -1 <filename.fastq> -r </path/to/dir_fasta/> -o </name/output/directory/>
```

## Create methylKit objects

``` r
SRA <- c("SRR534263","SRR534264","SRR534265","SRR534266","WT_rep1","SRR534177","SRR534193")
list_samples <- as.list(c("vim1","vim2","vim3","vim123","WT_rep1","WT_rep2","WT_rep3"))
list_treatments <- as.numeric(rep(0,7))

#list_samples <- paste(SRA,"_",sample, sep="")
```

## Create methylRawListDB objects

``` r
import_bismark_cytosine_report(path_bismark_files, list_DB_paths, list_samples, list_treatments)
```

## Load methylRawListDB objects

``` r
# Load methylRawListDB objects (without filtering)
list_methylRawLists_raw <- load_methylRawListDB(list_DB_paths, type="raw", list_samples, list_treatments)
```

## Filter methylRawList raw

``` r
filter_methylRawList(list_methylRawLists_raw)
```

## Load filtered methylRawListDB objects

``` r
list_methylRawLists <- load_methylRawListDB(list_DB_paths, type="filtered", list_samples, list_treatments)
```

## Subset genomic regions

## Subset data

Get data for each regions. In this case we select TEs, genic and
intergenic regions.

``` r
# Create subset for methylRawList
subset_methylObject(list_methylRawLists, list_DB_paths, bed_genes, "genes", "methylRaw")

subset_methylObject(list_methylRawLists, list_DB_paths, bed_TEs, "TEs", "methylRaw")

subset_methylObject(list_methylRawLists, list_DB_paths, bed_TEs_4kb, "TEs_4kb", "methylRaw")
```

## Load methylRaw subset data

``` r
# Load methylRawListDB objects (without filtering)
list_methylRawLists_genes <- load_methylRawListDB(list_DB_paths, type="genes", list_samples, list_treatments)

list_methylRawLists_TEs <- load_methylRawListDB(list_DB_paths, type="TEs", list_samples, list_treatments)

list_methylRawLists_TEs_4kb <- load_methylRawListDB(list_DB_paths, type="TEs_4kb", list_samples, list_treatments)
```

## Genome-wide methylation

``` r
path_DB <- "data/methylKit_DB_files/stroud_2013/"

df_name <- "df_mean_filtered"
title <- "Weighted Methylation Genome-wide"

get_df_wml(list_methylRawLists, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
ggplot_all(get(df_name), title = title)
```

![](images/genome_wide_methylation_vim_mutants.png)

## Genes

``` r
df_name <- "df_mean_genes"
title <- "Weighted Methylation Level for genes"

get_df_wml(list_methylRawLists_genes, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
ggplot_all(get(df_name), title = title)
```

![](images/gene_methylation_vim_mutants.png)

``` r
df_mean_genes %>% filter(context=="CpG") %>% aggregate(percent_methylation ~ sample,  
          FUN = function(x) c(mean = mean(x), sd = sd(x)))
```

``` r
aggregate(percent_methylation ~ sample, data = df_mean_genes, FUN = mean)
```

## All TEs

``` r
df_name <- "df_mean_TEs"
title <- "Weighted Methylation Level for all TEs"

get_df_wml(list_methylRawLists_TEs, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
ggplot_all(get(df_name), title = title)
```

![](images/all_TEs_methylation_vim_mutants.png)

## Long TEs

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Plot (use get() to pass the string name of the dataframe as a R object)
ggplot_all(get(df_name), title = title)
```

![](images/long_TEs_methylation_vim_mutants.png)

## Whole genome methylation compared to SA

``` r
df_accessions <- read.table("data/df_accessions_all.txt", header = TRUE, stringsAsFactors = TRUE)

path_DB <- "data/methylKit_DB_files/stroud_2013/"
#path_DB <- "F:/NETSCRATCH/methylKit_DB_files/GC_3427_3542_3599_4050_4220_4373_TAIR10"

df_name <- "df_mean_filtered_134"
title <- "Weighted Methylation Level genome-wide"

load_df_wml(path_DB, df_name)

# Merge dataframes of methylation levels and df_accessions to get detailed information about dataset
df_mean <- merge(get(df_name), df_accessions, by="sample")

# Keep only CpG context and remove cmt2-5
df_mean_CpG <- df_mean %>% filter(context=="CpG") %>% filter(location %in% c("Col-0","Col-0_Stroud","vim1","vim2","vim3","vim123","SA"))

# Arrange by name on the plot
df_mean_CpG$location <- factor(df_mean_CpG$location, levels = c("Col-0","Col-0_Stroud","vim1","vim2","vim3","vim123","SA"), ordered = TRUE)

df_mean_CpG$VIM2 <- as.factor(df_mean_CpG$VIM2)


# Highlight Cvi-0 within SA (84th samples in dataframe)
ggplot(data = df_mean_CpG, aes(location, percent_methylation, fill = VIM2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = df_mean_CpG[-84, ], position = position_jitterdodge(), size = 0.5) +
  geom_jitter(data = df_mean_CpG[84, ], size = 2, fill = "white", shape = 23, width = 0.1, height = 0) +
  ylab("% of methyated cytosines (CG)") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 2)) +
  scale_x_discrete(labels = c("Col-0" = "Col-0", "Col-0_Stroud" = "Col-0", "vim1" = "vim1-2", "vim2" = "vim2-1", "vim3" = "vim3-1", "vim123" = "vim1-2/vim2-1/vim3-1", "SA" = "CPV")) +
  scale_fill_manual(values = c("#ffffffff", "#949494ff")) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::comma_format())
```

![](images/genome_wide_methylation_vim_mutants_with_SA.png)

## Gene body methylation compared to SA

``` r
df_accessions <- read.table("data/df_accessions_all.txt", header = TRUE, stringsAsFactors = TRUE)

path_DB <- "data/methylKit_DB_files/stroud_2013/"
#path_DB <- "F:/NETSCRATCH/methylKit_DB_files/GC_3427_3542_3599_4050_4220_4373_TAIR10"

df_name <- "df_mean_genes_134"
title <- "Weighted Methylation Level for genes"

load_df_wml(path_DB, df_name)

# Merge dataframes of methylation levels and df_accessions to get detailed information about dataset
df_mean <- merge(get(df_name), df_accessions, by = "sample")

# Keep only CpG context and remove cmt2-5
df_mean_CpG <- df_mean %>%
  filter(context == "CpG") %>%
  filter(location %in% c("Col-0", "Col-0_Stroud", "vim1", "vim2", "vim3", "vim123", "SA"))

# Arrange by name on the plot
df_mean_CpG$location <- factor(df_mean_CpG$location, levels = c("Col-0", "Col-0_Stroud", "vim1", "vim2", "vim3", "vim123", "SA"), ordered = TRUE)

df_mean_CpG$VIM2 <- as.factor(df_mean_CpG$VIM2)


# Highlight Cvi-0 within SA (84th samples in dataframe)
ggplot(data = df_mean_CpG, aes(location, percent_methylation, fill = VIM2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = df_mean_CpG[-84, ], position = position_jitterdodge(), size = 0.5) +
  geom_jitter(data = df_mean_CpG[84, ], size = 2, fill = "white", shape = 23, width = 0.1, height = 0) +
  ylab("% of methyated cytosines (CG)") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 2)) +
  scale_x_discrete(labels = c("Col-0" = "Col-0", "Col-0_Stroud" = "Col-0", "vim1" = "vim1-2", "vim2" = "vim2-1", "vim3" = "vim3-1", "vim123" = "vim1-2/vim2-1/vim3-1", "SA" = "CPV")) +
  scale_fill_manual(values = c("#ffffffff", "#949494ff")) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::comma_format())
```

![](images/gbM_methylation_vim_mutants_with_SA.png)

## TEs

``` r
df_accessions <- read.table("data/df_accessions_all.txt", header = TRUE, stringsAsFactors = TRUE)

path_DB <- "data/methylKit_DB_files/stroud_2013/"
#path_DB <- "F:/NETSCRATCH/methylKit_DB_files/GC_3427_3542_3599_4050_4220_4373_TAIR10"

df_name <- "df_mean_TEs_134"
title <- "Weighted Methylation Level all TEs"

load_df_wml(path_DB, df_name)

# Merge dataframes of methylation levels and df_accessions to get detailed information about dataset
df_mean <- merge(get(df_name), df_accessions, by = "sample")

# Keep only CpG context and remove cmt2-5
df_mean_CpG <- df_mean %>%
  filter(context == "CpG") %>%
  filter(location %in% c("Col-0", "Col-0_Stroud", "vim1", "vim2", "vim3", "vim123", "SA"))

# Arrange by name on the plot
df_mean_CpG$location <- factor(df_mean_CpG$location, levels = c("Col-0", "Col-0_Stroud", "vim1", "vim2", "vim3", "vim123", "SA"), ordered = TRUE)

df_mean_CpG$VIM2 <- as.factor(df_mean_CpG$VIM2)

# Highlight Cvi-0 within SA (84th samples in dataframe)
ggplot(data = df_mean_CpG, aes(location, percent_methylation, fill = VIM2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = df_mean_CpG[-84, ], position = position_jitterdodge(), size = 0.5) +
  geom_jitter(data = df_mean_CpG[84, ], size = 2, fill = "white", shape = 23, width = 0.1, height = 0) +
  ylab("% of methyated cytosines (CG)") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 2)) +
  scale_x_discrete(labels = c("Col-0" = "Col-0", "Col-0_Stroud" = "Col-0", "vim1" = "vim1-2", "vim2" = "vim2-1", "vim3" = "vim3-1", "vim123" = "vim1-2/vim2-1/vim3-1", "SA" = "CPV")) +
  scale_fill_manual(values = c("#ffffffff", "#949494ff")) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::comma_format())
```

![](images/all_TEs_methylation_vim_mutants_with_SA.png)

# Analysis DNA methylation in SAIL and SALK fbx5 and cmt2 mutants

## Samples

Note that we use here ARA1 (short for ARABIDILLO-1) in place of FBX5. We
changed this in the manuscript since ARA1 is already allocated to
another gene (AT4G16130) and FBX5 was the first reference for the gene
AT2G44900.

ara12 stands for the double mutant ARABIDILLO-1/FBX5 (AT2G44900) and
ARABIDILLO-2 ( AT3G60350). Details on the SRA for each library is
available on Supplementary table 1 of the paper.

| name          | library | sample               | pool      |
|---------------|---------|----------------------|-----------|
| ara12_I       | 4373_A  | 4373_A_ara12_I       | ara12     |
| ara12_II      | 4373_B  | 4373_B_ara12_II      | ara12     |
| ara12_III     | 4373_C  | 4373_C_ara12_III     | ara12     |
| ara1-SALK_I   | 4373_D  | 4373_D_ara1-SALK_I   | ara1-SALK |
| ara1-SALK_II  | 4373_E  | 4373_E_ara1-SALK_II  | ara1-SALK |
| ara1-SALK_III | 4373_F  | 4373_F_ara1-SALK_III | ara1-SALK |
| ARA1-OE_I     | 4373_G  | 4373_G_ARA1-OE_I     | ARA1-OE   |
| ARA1-OE_II    | 4373_H  | 4373_H_ARA1-OE_II    | ARA1-OE   |
| ARA1-OE_III   | 4373_I  | 4373_I_ARA1-OE_III   | ARA1-OE   |
| ara1-SAIL_I   | 4373_J  | 4373_J_ara1-SAIL_I   | ara1-SAIL |
| ara1-SAIL_II  | 4373_K  | 4373_K_ara1-SAIL_II  | ara1-SAIL |
| ara1-SAIL_III | 4373_L  | 4373_L_ara1-SAIL_III | ara1-SAIL |
| cmt2-5_I      | 4373_M  | 4373_M_cmt2-5_I      | cmt2-5    |
| cmt2-5_II     | 4373_N  | 4373_N_cmt2-5_II     | cmt2-5    |
| cmt2-5_III    | 4373_O  | 4373_O_cmt2-5_III    | cmt2-5    |
| Col-3_I       | 4373_P  | 4373_P_Col-3_I       | Col-3     |
| Col-3_II      | 4373_Q  | 4373_Q_Col-3_II      | Col-3     |
| Col-3_III     | 4373_R  | 4373_R_Col-3_III     | Col-3     |
| Col-0_I       | 4373_S  | 4373_S_Col-0_I       | Col-0     |
| Col-0_II      | 4373_T  | 4373_T_Col-0_II      | Col-0     |
| Col-0_III     | 4373_U  | 4373_U_Col-0_III     | Col-0     |
| Cvi-0_I       | 4373_V  | 4373_V_Cvi-0_I       | Cvi-0     |
| Cvi-0_II      | 4373_W  | 4373_W_Cvi-0_II      | Cvi-0     |
| Cvi-0_III     | 4373_X  | 4373_X_Cvi-0_III     | Cvi-0     |

## Run Bismark

For each fastq file, run the following command:

``` bash
bash run_bismark.sh -1 <filename.fastq> -r </path/to/dir_fasta/> -o </name/output/directory/>
```

## Create methylKit objects

## Create methylRawListDB objects

``` r
import_bismark_cytosine_report(path_bismark_files, list_DB_paths, list_samples, list_treatments)
```

## Load methylRawListDB objects

``` r
# Load methylRawListDB objects (without filtering)
list_methylRawLists_raw <- load_methylRawListDB(list_DB_paths, type="raw", list_samples, list_treatments)
```

## Filter methylRawList raw

``` r
filter_methylRawList(list_methylRawLists_raw)
```

## Load filtered methylRawListDB objects

``` r
list_methylRawLists <- load_methylRawListDB(list_DB_paths, type="filtered", list_samples, list_treatments)
```

## Subset genomic regions

## Subset data

Get data for each regions. In this case we select TEs, genic and
intergenic regions.

``` r
# Create subset for methylRawList
subset_methylObject(list_methylRawLists, list_DB_paths, bed_genes, "genes", "methylRaw")

subset_methylObject(list_methylRawLists, list_DB_paths, bed_TEs, "TEs", "methylRaw")

subset_methylObject(list_methylRawLists, list_DB_paths, bed_TEs_4kb, "TEs_4kb", "methylRaw")
```

## Load methylRaw subset data

``` r
# Load methylRawListDB objects (without filtering)
list_methylRawLists_genes <- load_methylRawListDB(list_DB_paths, type="genes", list_samples, list_treatments)

list_methylRawLists_TEs <- load_methylRawListDB(list_DB_paths, type="TEs", list_samples, list_treatments)

list_methylRawLists_TEs_4kb <- load_methylRawListDB(list_DB_paths, type="TEs_4kb", list_samples, list_treatments)
```

## mCG in long TEs for FBX5 mutants

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)

load_df_wml(path_DB, df_name)


# Merge dataframes of methylation levels and df_accessions to get detailed information about dataset
df_mean <- merge(get(df_name), df_accessions, by = "sample")


df_mean_subset <- df_mean %>%
  filter(context == "CpG") %>%
  filter(!pool %in% c("cmt2-5", "cmt2", "S27-204", "Cvi-0"))

order_pool <- c("Col-0", "ara1-SALK", "ARA1-OE", "Col-3", "ara1-SAIL", "ara12")

df_mean_subset$pool <- factor(df_mean_subset$pool, levels = order_pool, ordered = TRUE)

# Rename samples

ggplot(data = df_mean_subset, aes(x = pool, y = percent_methylation, group = pool)) +
  geom_boxplot(outlier.size = 1) +
  geom_jitter(height = .05, width = .05, size = 1) +
  theme_bw() +
  ylab("% methylated cytosines (CG)") +
  ggtitle("CpG methylation in long TEs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  scale_y_continuous(limits = c(78, 90), breaks = seq(78, 90, 2)) +
  scale_x_discrete(labels = c("Col-0", "SALK_082977", "FBX5-OE", "Col-3", "SAIL_190_D02", "SAILD_190_D02/arabidillo-2"))
```

![](images/mCG_long_TEs_fbx5_mutants.png)

## mCHG in long TEs for FBX5 mutants

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Merge dataframes of methylation levels and df_accessions to get detailed information about dataset
df_mean <- merge(get(df_name), df_accessions, by = "sample")

df_mean_subset <- df_mean %>%
  filter(context == "CHG") %>%
  filter(!pool %in% c("cmt2-5", "cmt2", "S27-204", "Cvi-0"))

order_pool <- c("Col-0", "ara1-SALK", "ARA1-OE", "Col-3", "ara1-SAIL", "ara12")

df_mean_subset$pool <- factor(df_mean_subset$pool, levels = order_pool, ordered = TRUE)

# Rename samples

ggplot(data = df_mean_subset, aes(x = pool, y = percent_methylation, group = pool)) +
  geom_boxplot(outlier.size = 1) +
  geom_jitter(height = .05, width = .05, size = 1) +
  theme_bw() +
  ylab("% methylated cytosines (CHG)") +
  ggtitle("CHG methylation in long TEs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  scale_y_continuous(limits = c(40, 48), breaks = seq(40, 48, 2)) +
  scale_x_discrete(labels = c("Col-0", "SALK_082977", "FBX5-OE", "Col-3", "SAIL_190_D02", "SAILD_190_D02/arabidillo-2"))
```

![](images/fbx5_mCHG_longTEs.png)

## mCHH in long TEs for FBX5 mutants

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)
load_df_wml(path_DB, df_name)

# Merge methylation levels and df_accessions for detailed information
df_mean <- merge(get(df_name), df_accessions, by = "sample")

df_mean_subset <- df_mean %>%
  filter(context == "CHH") %>%
  filter(!pool %in% c("cmt2-5", "cmt2", "S27-204", "Cvi-0"))

order_pool <- c("Col-0", "ara1-SALK", "ARA1-OE", "Col-3", "ara1-SAIL", "ara12")

df_mean_subset$pool <- factor(df_mean_subset$pool, levels = order_pool, ordered = TRUE)

ggplot(data = df_mean_subset, aes(x = pool, y = percent_methylation, group = pool)) +
  geom_boxplot(outlier.size = 1) +
  geom_jitter(height = .05, width = .05, size = 1) +
  theme_bw() +
  ylab("% methylated cytosines (CHH)") +
  ggtitle("CHH methylation in long TEs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  scale_y_continuous(limits = c(8, 14), breaks = seq(8, 14, 2)) +
  scale_x_discrete(labels = c(
    "Col-0", "SALK_082977", "FBX5-OE", "Col-3",
    "SAIL_190_D02", "SAILD_190_D02/arabidillo-2"
  ))
```

![](images/fbx5_mCHH_long_TEs.png)

## mCHH in whole genome for FBX5 mutant

``` r
df_name <- "df_mean_filtered"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)
load_df_wml(path_DB, df_name)

# Merge methylation levels and df_accessions for detailed information
df_mean <- merge(get(df_name), df_accessions, by = "sample")

df_mean_subset <- df_mean %>%
  filter(context == "CHH") %>%
  filter(!pool %in% c("cmt2-5", "cmt2", "S27-204", "Cvi-0"))

order_pool <- c("Col-0", "ara1-SALK", "ARA1-OE", "Col-3", "ara1-SAIL", "ara12")

df_mean_subset$pool <- factor(df_mean_subset$pool, levels = order_pool, ordered = TRUE)

ggplot(data = df_mean_subset, aes(x = pool, y = percent_methylation, group = pool)) +
  geom_boxplot(outlier.size = 1) +
  geom_jitter(height = .05, width = .05, size = 1) +
  theme_bw() +
  ylab("% methylated cytosines (CHH)") +
  ggtitle("CHH methylation in long TEs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  scale_y_continuous(limits = c(8, 14), breaks = seq(8, 14, 2)) +
  scale_x_discrete(labels = c(
    "Col-0", "SALK_082977", "FBX5-OE", "Col-3",
    "SAIL_190_D02", "SAILD_190_D02/arabidillo-2"
  ))
```

## Statistical tests

Check significance of difference between Col-0 and ara1-SALK and
ARA1-OE, and Col-3 and ara12 and ara1-SAIL

### Col-0 background

#### mCG

``` r
require(onewaytests)
require(multcomp)

df <- df_mean %>% filter(context=="CpG") %>% filter(pool %in% c("Col-0","ara1-SALK","ARA1-OE"))

bartlett.test(percent_methylation~pool, data=df)

lm1 <- aov(percent_methylation ~ pool, data=df)
summary(glht(lm1, linfct=mcp(pool="Tukey"), alternative="two.sided"))
```

Bartlett test of homogeneity of variances

data: percent_methylation by pool Bartlett’s K-squared = 2.4966, df = 2,
p-value = 0.287

     Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Tukey Contrasts

Fit: aov(formula = percent_methylation ~ pool, data = df)

Linear Hypotheses: Estimate Std. Error t value Pr(\>\|t\|)  
ara1-SALK - ARA1-OE == 0 5.0455 0.5667 8.903 \< 0.001 \*** Col-0 -
ARA1-OE == 0 1.2969 0.5667 2.288 0.13379  
Col-0 - ara1-SALK == 0 -3.7486 0.5667 -6.614 0.00127 **

#### mCHG

``` r
require(onewaytests)
require(multcomp)

df <- df_mean %>% filter(context=="CHG") %>% filter(pool %in% c("Col-0","ara1-SALK","ARA1-OE"))

bartlett.test(percent_methylation~pool, data=df)

lm1 <- aov(percent_methylation ~ pool, data=df)
summary(glht(lm1, linfct=mcp(pool="Tukey"), alternative="two.sided"))
```

Bartlett test of homogeneity of variances

data: percent_methylation by pool Bartlett’s K-squared = 0.29114, df =
2, p-value = 0.8645

     Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Tukey Contrasts

Fit: aov(formula = percent_methylation ~ pool, data = df)

Linear Hypotheses: Estimate Std. Error t value Pr(\>\|t\|)  
ara1-SALK - ARA1-OE == 0 1.8133 0.7679 2.361 0.1218  
Col-0 - ARA1-OE == 0 -0.2800 0.7679 -0.365 0.9302  
Col-0 - ara1-SALK == 0 -2.0933 0.7679 -2.726 0.0768 .

#### mCHH

``` r
require(onewaytests)
require(multcomp)

df <- df_mean %>% filter(context=="CHH") %>% filter(pool %in% c("Col-0","ara1-SALK","ARA1-OE"))

bartlett.test(percent_methylation~pool, data=df)

lm1 <- aov(percent_methylation ~ pool, data=df)
summary(glht(lm1, linfct=mcp(pool="Tukey"), alternative="two.sided"))
```

    Bartlett test of homogeneity of variances

data: percent_methylation by pool Bartlett’s K-squared = 0.28809, df =
2, p-value = 0.8658

     Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Tukey Contrasts

Fit: aov(formula = percent_methylation ~ pool, data = df)

Linear Hypotheses: Estimate Std. Error t value Pr(\>\|t\|)  
ara1-SALK - ARA1-OE == 0 -0.5501 0.3617 -1.521 0.3469  
Col-0 - ARA1-OE == 0 -1.1019 0.3617 -3.046 0.0511 . Col-0 - ara1-SALK ==
0 -0.5518 0.3617 -1.525 0.3450

### Col-3 background

#### mCG

``` r
require(onewaytests)
require(multcomp)

df <- df_mean %>% filter(context=="CpG") %>% filter(pool %in% c("Col-3","ara1-SAIL","ara12"))

bartlett.test(percent_methylation~pool, data=df)

lm1 <- aov(percent_methylation ~ pool, data=df)
summary(glht(lm1, linfct=mcp(pool="Tukey"), alternative="two.sided"))
```

    Bartlett test of homogeneity of variances

data: percent_methylation by pool Bartlett’s K-squared = 0.77576, df =
2, p-value = 0.6785

     Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Tukey Contrasts

Fit: aov(formula = percent_methylation ~ pool, data = df)

Linear Hypotheses: Estimate Std. Error t value Pr(\>\|t\|)  
ara12 - ara1-SAIL == 0 1.3206 0.2585 5.108 0.0052 \*\* Col-3 - ara1-SAIL
== 0 -4.3065 0.2585 -16.657 \<0.001 *** Col-3 - ara12 == 0 -5.6272
0.2585 -21.765 \<0.001 ***

#### mCHG

``` r
require(onewaytests)
require(multcomp)

df <- df_mean %>% filter(context=="CHG") %>% filter(pool %in% c("Col-3","ara1-SAIL","ara12"))

bartlett.test(percent_methylation~pool, data=df)

lm1 <- aov(percent_methylation ~ pool, data=df)
summary(glht(lm1, linfct=mcp(pool="Tukey"), alternative="two.sided"))
```

    Bartlett test of homogeneity of variances

data: percent_methylation by pool Bartlett’s K-squared = 0.056846, df =
2, p-value = 0.972

     Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Tukey Contrasts

Fit: aov(formula = percent_methylation ~ pool, data = df)

Linear Hypotheses: Estimate Std. Error t value Pr(\>\|t\|)  
ara12 - ara1-SAIL == 0 0.3367 0.7415 0.454 0.8945  
Col-3 - ara1-SAIL == 0 -2.6067 0.7415 -3.516 0.0295 * Col-3 - ara12 == 0
-2.9433 0.7415 -3.970 0.0173 *

#### mCHH

``` r
require(onewaytests)
require(multcomp)

df <- df_mean %>% filter(context=="CHH") %>% filter(pool %in% c("Col-3","ara1-SAIL","ara12"))

bartlett.test(percent_methylation~pool, data=df)

lm1 <- aov(percent_methylation ~ pool, data=df)
summary(glht(lm1, linfct=mcp(pool="Tukey"), alternative="two.sided"))
```

Bartlett test of homogeneity of variances

data: percent_methylation by pool Bartlett’s K-squared = 0.46264, df =
2, p-value = 0.7935

     Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Tukey Contrasts

Fit: aov(formula = percent_methylation ~ pool, data = df)

Linear Hypotheses: Estimate Std. Error t value Pr(\>\|t\|) ara12 -
ara1-SAIL == 0 -0.9421 0.7153 -1.317 0.437 Col-3 - ara1-SAIL == 0
-1.0532 0.7153 -1.472 0.367 Col-3 - ara12 == 0 -0.1111 0.7153 -0.155
0.987

## mCHH in long TEs for cmt2-5 mutants

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs (>4 kb)"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Merge dataframes of methylation levels and df_accessions to get detailed information about dataset
df_mean <- merge(get(df_name), df_accessions, by="sample")

df_mean_subset <- df_mean %>% filter(pool %in% c("cmt2-5","Col-3")) %>% filter(context=="CHH")

order_pool <- c("Col-3","cmt2-5")

df_mean_subset$pool <- factor(df_mean_subset$pool , levels=order_pool, ordered=TRUE)


ggplot(data=df_mean_subset[(df_mean_subset$context=="CHH"),], aes(x=pool, y=percent_methylation, group=pool)) + 
  geom_boxplot(outlier.size=1) +
  geom_jitter(height=.05, width=.05, size=1) +  
  theme_bw() +
  ylab("% of methylated cytosines (CHH)") +
  ggtitle("CpG methylation whole genome") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + scale_y_continuous(limits=c(0, 14), breaks=seq(0,14,2)) +
  scale_x_discrete(labels=c("Col-3","cmt2-5 SAIL"))
```

![](images/mCHH_long_TE_Col3_cmt25.png)

### Statistical test

``` r
# Check variance
bartlett.test(percent_methylation~pool, data = df_mean_subset)
#qqnorm(df_mean_subset$percent_methylation); qqline(df_mean_subset$percent_methylation)

# Variance and normality OK

# Alt Hypothesis is ara1 mCG > Col-0 mCG 
with(df_mean_subset, t.test(percent_methylation[pool=="Col-3"], percent_methylation[pool=="cmt2-5"], alternative="two.sided", var.equal=TRUE))
```

Bartlett test of homogeneity of variances

data: percent_methylation by pool Bartlett’s K-squared = 2.337, df = 1,
p-value = 0.1263

    Two Sample t-test

data: percent_methylation\[pool == “Col-3”\] and
percent_methylation\[pool == “cmt2-5”\] t = 19.682, df = 4, p-value =
3.93e-05 alternative hypothesis: true difference in means is not equal
to 0 95 percent confidence interval: 8.07735 10.73048 sample estimates:
mean of x mean of y 11.26492 1.86100

## Compare cmt2 and SA with CMT2stop

``` r
df_accessions <- read.table("data/df_accessions_all.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE, na.strings = "")

path_DB <- "F:/NETSCRATCH/methylKit_DB_files/GC_3427_3542_3599_4050_4220_4373_TAIR10"

df_name <- "df_mean_TEs_4kb_161"
title <- "Weighted Methylation Level for long TEs"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)

load_df_wml(path_DB, df_name)

# Merge dataframes of methylation levels and df_accessions to get detailed information about dataset
df_mean <- merge(get(df_name), df_accessions, by = "sample")

# Keep only CHH
df_mean_formated_CHH <- df_mean %>%
  filter(context == "CHH") %>%
  filter(location %in% c("Col-3", "cmt2-5", "SA")) %>%
  filter(!sample %in% c("4220_AE_cmt2-5", "3427_C_Col-0", "3599_A_Col-0", "4220_Z_Col-0", "4220_AG_Col-3"))

# Arrange by name on the plot
df_mean_formated_CHH$location <- factor(df_mean_formated_CHH$location, levels = c("Col-3", "cmt2-5", "SA"), ordered = TRUE)

df_mean_formated_CHH$CMT2 <- as.factor(df_mean_formated_CHH$CMT2)

# Create new variable to separate SA with CMT2stop and CMT2ref on the plot
df_SA_CMT2ref <- df_mean_formated_CHH %>% filter(location == "SA", CMT2 == 0)
df_SA_CMT2ref$pool <- "SA_CMT2ref"

df_SA_CMT2stop <- df_mean_formated_CHH %>% filter(location == "SA", CMT2 == 1)
df_SA_CMT2stop$pool <- "SA_CMT2stop"

df_col <- df_mean_formated_CHH %>% filter(location == "Col-3")
df_col$pool <- "Col-3"

df_cmt2 <- df_mean_formated_CHH %>% filter(location == "cmt2-5")
df_cmt2$pool <- "cmt2-5"

# Rejoin dataframe
df_mean_formated_CHH <- rbind(df_SA_CMT2ref, df_SA_CMT2stop, df_col, df_cmt2)

# Reorder factor
df_mean_formated_CHH$pool <- factor(df_mean_formated_CHH$pool, levels = c("SA_CMT2ref", "SA_CMT2stop", "Col-3", "cmt2-5"))

# Where Cvi-0 is
which(grepl("4220_Y_Cvi-0", df_mean_formated_CHH$sample))
# 56

give.n <- function(x) {
  return(c(y = mean(x), label = length(x)))
}

ggplot(data = df_mean_formated_CHH, aes(x = pool, y = percent_methylation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = df_mean_formated_CHH[-56, ], height = .05, width = .05, size = 0.8) +
  geom_jitter(data = df_mean_formated_CHH[56, ], size = 2, fill = "white", shape = 23, width = 0.1, height = 0) +
  ylab("% of methylated cytosines (CHH)") +
  ggtitle("CHH methylation in long TEs") +
  xlab("") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_summary(fun.data = give.n, geom = "text", position = position_stack(vjust = 0)) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::comma_format())
```

![](images/mCHH_long_TEs_SA_cmt2_mutant.png)

### Statistical difference

``` r
require(onewaytests)
# Compare SA CMT2

df_mean_formated_1 <- df_mean_formated_CHH %>% filter(!location %in% c("Col-3","cmt2-5"))

# Test if variances equals
bartlett.test(percent_methylation~VIM2, data=df_mean_formated_1)
# Variances are equal

with(df_mean_formated_1, t.test(percent_methylation[CMT2==0], percent_methylation[CMT2==1], var.equal=TRUE, alternative="two.sided"))
# <2.2e-16

# Compare cmt2-5 and Col-3
df_mean_formated_2 <- df_mean_formated_CHH %>% filter(location %in% c("Col-3","cmt2-5"))

# Test if variances equals
bartlett.test(percent_methylation~pool, data=df_mean_formated_2)
# Variances are equal

with(df_mean_formated_2, t.test(percent_methylation[pool=="Col-3"], percent_methylation[pool=="cmt2-5"], var.equal=TRUE, alternative="two.sided"))
# 3.93e-05
```

Bartlett test of homogeneity of variances

data: percent_methylation by VIM2 Bartlett’s K-squared = 0.85006, df =
1, p-value = 0.3565

    Two Sample t-test

data: percent_methylation\[CMT2 == 0\] and percent_methylation\[CMT2 ==
1\] t = 33.064, df = 81, p-value \< 2.2e-16 alternative hypothesis: true
difference in means is not equal to 0 95 percent confidence interval:
8.016816 9.043461 sample estimates: mean of x mean of y 11.15977 2.62963

    Bartlett test of homogeneity of variances

data: percent_methylation by pool Bartlett’s K-squared = 2.337, df = 1,
p-value = 0.1263

    Two Sample t-test

data: percent_methylation\[pool == “Col-3”\] and
percent_methylation\[pool == “cmt2-5”\] t = 19.682, df = 4, p-value =
3.93e-05 alternative hypothesis: true difference in means is not equal
to 0 95 percent confidence interval: 8.07735 10.73048 sample estimates:
mean of x mean of y 11.26492 1.86100

# Analysis DNA methylation fbx5 CRISPR lines

## Samples

| Library title | name_R | genotype | Library number | treatment | conversion_efficiency | total reads | final_usable |
|----|----|----|----|----|----|----|----|
| ara1 \#12-3 | fbx5_12_3_rep1 | fbx5_12_3 | 5007_C | 1 | 98.77% | 7,661,742.0 | 34.0% |
| ara1 \#12-3 rep2 | fbx5_12_3_rep2 | fbx5_12_3 | 5083_E | 1 | 99.10% | 7,390,669.00 | 38.1% |
| ara1 \#12-3 rep3 | fbx5_12_3_rep3 | fbx5_12_3 | 5083_F | 1 | 99.20% | 7,459,641.00 | 37.9% |
| ara1 \#3-23 | fbx5_3_23_rep1 | fbx5_3_23 | 5007_B | 1 | 98.75% | 7,597,128.0 | 33.9% |
| ara1 \#3-23 rep2 | fbx5_3_23_rep2 | fbx5_3_23 | 5083_C | 1 | 98.78% | 7,495,590.00 | 39.0% |
| ara1 \#3-23 rep3 | fbx5_3_23_rep3 | fbx5_3_23 | 5083_D | 1 | 98.97% | 7,375,831.00 | 37.6% |
| ara1 \#8-9 | fbx5_8_9_rep1 | fbx5_8_9 | 5007_A | 1 | 98.97% | 7,370,862.0 | 34.5% |
| ara 1 \# 8-9 rep2 | fbx5_8_9_rep2 | fbx5_8_9 | 5083_O | 1 | 98.92% | 5,153,243.00 | 20.2% |
| ara 1 \# 8-9 rep3 | fbx5_8_9_rep3 | fbx5_8_9 | 5083_P | 1 | 99.11% | 5,677,122.00 | 16.7% |
| S7-B5 | S7_B5_rep1 | S7_B5 | 5007_D | 0 | 99.05% | 7,629,428.0 | 33.7% |
| S7-B5 rep2 | S7_B5_rep2 | S7_B5 | 5007_E | 0 | 99.04% | 7,456,458.0 | 33.8% |
| S7-B5 rep3 | S7_B5_rep3 | S7_B5 | 5083_G | 0 | 99.16% | 7,491,817.00 | 39.2% |

## Create methylRawListDB objects

Upload the files in my local PC to make things faster. I just need the
files `*_report_only_chr.txt`. Omit CX files.

``` r
import_bismark_cytosine_report(path_bismark_files, list_DB_paths, list_samples, list_treatments)
```

## Load methylRawListDB objects

Once created, load methylRawListDB objects. The files won’t actually be
loaded but accessed in real time when needed.

``` r
list_methylRawLists_raw <- load_methylRawListDB(list_DB_paths, type="raw", list_samples, list_treatments)
```

## Filter methylRawList raw

Keep only cytosine positions that have a define minimum coverage. This
threshold is usually set at around 5 in most WGBS analyses but since our
samples were sequenced at the minimum depth allowed by the sequencing
facility, we defined a lower threshold (minimum 2). This approach is
valid considering that we look at pattern across large genomic regions.
We assumed we would catch any strong signal if any.

``` r
filter_methylRawList(list_methylRawLists_raw)
```

## Load filtered methylRawListDB objects

``` r
list_methylRawLists <- load_methylRawListDB(list_DB_paths, type="filtered", list_samples, list_treatments)
```

## Subset genomic regions

We want now to analyze methylation patterns is specific genomic regions.
For this, we need to subset our data and generate new DB flat files for
the different regions.

## Subset data

``` r
# Create subset for methylRawList
subset_methylObject(list_methylRawLists, list_DB_paths, bed_genes, "genes", "methylRaw")

subset_methylObject(list_methylRawLists, list_DB_paths, bed_TEs, "TEs", "methylRaw")

subset_methylObject(list_methylRawLists, list_DB_paths, bed_TEs_4kb, "TEs_4kb", "methylRaw")

subset_methylObject(list_methylRawLists, list_DB_paths, bed_TEs_500bp, "TEs_500bp", "methylRaw")
```

## Load methylRaw objects per regions

``` r
# Load methylRawListDB objects 
list_methylRawLists_genes <- load_methylRawListDB(list_DB_paths, type="genes", list_samples, list_treatments)

list_methylRawLists_TEs <- load_methylRawListDB(list_DB_paths, type="TEs", list_samples, list_treatments)

list_methylRawLists_TEs_4kb <- load_methylRawListDB(list_DB_paths, type="TEs_4kb", list_samples, list_treatments)

list_methylRawLists_TEs_500bp <- load_methylRawListDB(list_DB_paths, type="TEs_500bp", list_samples, list_treatments)
```

## Create methlBaseDB objects

Create methylBase objects based on given list_methylRawLists object.
Create DB files if not existing.

``` r
list_methylBases <- merged_methylRawList(list_methylRawLists)
```

## Load methylBaseDB objects

``` r
# Make lists of objects
list_methylBases <- load_methylBaseDB(list_DB_paths, list_samples, list_treatments)
```

## Long TE methylation

### mCG

``` r
df_name <- "df_mean_TEs_4kb"
title <- "Weighted Methylation Level for long TEs"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)

load_df_wml(path_DB, df_name)

df_mean_TEs_4kb <- df_mean_TEs_4kb %>% filter(sample %in% as.vector(paste(list_samples, "_TEs_4kb", sep = "")))

# Create a vector of 3 times the genotype names

# Create variable without the rep
groups <- df_accessions %>%
  separate(genotype, c("genotype", "rep"), "_rep") %>%
  mutate_at("genotype", as.factor) %>%
  pull("genotype")

# add genotype
df_mean_TEs_4kb$genotype <- groups


df_mean_TEs_4kb$sample <- rep(df_accessions$genotype, 3)

# Plot (use get() to pass the string name of the dataframe as a R object)
# For CG
df_fbx5_CG <- df_mean_TEs_4kb %>%
  filter(genotype %in% c("fbx5_8_9", "fbx5_3_23", "fbx5_12_3", "S7_B5")) %>%
  filter(context == "CpG") %>%
  droplevels()

# Put WT sample as control
df_fbx5_CG$genotype <- factor(df_fbx5_CG$genotype, levels = c("S7_B5", "fbx5_8_9", "fbx5_12_3", "fbx5_3_23"), ordered = TRUE)


ggboxplot(df_fbx5_CG, x = "genotype", y = "percent_methylation", add = "jitter", title = "CpG methylation in long TEs") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + scale_x_discrete(labels = c("S7-B5", "fbx5-1", "fbx5-2", "fbx5-3"))
```

![](images/mCG_long_TEs_CRISPR_fbx5.png)

#### Statistical test

``` r
lm1 <- aov(percent_methylation ~ genotype, data=df_fbx5_CG)
summary(glht(lm1, linfct=mcp(genotype="Dunnett"), alternative="two.sided"))
```

Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Dunnett Contrasts

Fit: aov(formula = percent_methylation ~ genotype, data = df_fbx5_CG)

Linear Hypotheses: Estimate Std. Error t value Pr(\>\|t\|)  
fbx5_8_9 - S7_B5 == 0 2.0900 0.4066 5.140 0.00248 ** fbx5_12_3 - S7_B5
== 0 2.2600 0.4066 5.558 0.00152 ** fbx5_3_23 - S7_B5 == 0 1.0833 0.4066
2.664 0.06841 .

### mCHG

``` r
df_name <- "df_mean_TEs_4kb"

title <- "Weighted Methylation Level for long TEs"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)

load_df_wml(path_DB, df_name)

df_mean_TEs_4kb <- df_mean_TEs_4kb %>% filter(sample %in% as.vector(paste(list_samples, "_TEs_4kb", sep = "")))

# Create a vector of 3 times the genotype names

# Create variable without the rep
groups <- df_accessions %>%
  separate(genotype, c("genotype", "rep"), "_rep") %>%
  mutate_at("genotype", as.factor) %>%
  pull("genotype")

# add genotype
df_mean_TEs_4kb$genotype <- groups


df_mean_TEs_4kb$sample <- rep(df_accessions$genotype, 3)

# Plot (use get() to pass the string name of the dataframe as a R object)
# For CG
df_fbx5_CHG <- df_mean_TEs_4kb %>%
  filter(genotype %in% c("fbx5_8_9", "fbx5_3_23", "fbx5_12_3", "S7_B5")) %>%
  filter(context == "CHG") %>%
  droplevels()

# Put WT sample as control
df_fbx5_CHG$genotype <- factor(df_fbx5_CHG$genotype, levels = c("S7_B5", "fbx5_8_9", "fbx5_12_3", "fbx5_3_23"), ordered = TRUE)


ggboxplot(df_fbx5_CHG, x = "genotype", y = "percent_methylation", add = "jitter", title = "CHG methylation in long TEs") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + scale_x_discrete(labels = c("S7-B5", "fbx5-1", "fbx5-2", "fbx5-3"))
```

![](images/mCHG_long_TE_CRISPR_fbx5.png)

#### Statistical test

``` r
lm1 <- aov(percent_methylation ~ genotype, data=df_fbx5_CHG)
summary(glht(lm1, linfct=mcp(genotype="Dunnett"), alternative="two.sided"))
```

Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Dunnett Contrasts

Fit: aov(formula = percent_methylation ~ genotype, data = df_fbx5_CHG)

Linear Hypotheses: Estimate Std. Error t value Pr(\>\|t\|) fbx5_8_9 -
S7_B5 == 0 -0.7767 1.3100 -0.593 0.882 fbx5_12_3 - S7_B5 == 0 0.7433
1.3100 0.567 0.894 fbx5_3_23 - S7_B5 == 0 0.7800 1.3100 0.595 0.881

### mCHH

``` r
df_name <- "df_mean_TEs_4kb"

title <- "Weighted Methylation Level for long TEs"

get_df_wml(list_methylRawLists_TEs_4kb, path_DB, df_name)

load_df_wml(path_DB, df_name)

df_mean_TEs_4kb <- df_mean_TEs_4kb %>% filter(sample %in% as.vector(paste(list_samples, "_TEs_4kb", sep = "")))

# Create a vector of 3 times the genotype names

# Create variable without the rep
groups <- df_accessions %>%
  separate(genotype, c("genotype", "rep"), "_rep") %>%
  mutate_at("genotype", as.factor) %>%
  pull("genotype")

# add genotype
df_mean_TEs_4kb$genotype <- groups


df_mean_TEs_4kb$sample <- rep(df_accessions$genotype, 3)

# Plot (use get() to pass the string name of the dataframe as a R object)
df_fbx5_CHH <- df_mean_TEs_4kb %>%
  filter(genotype %in% c("fbx5_8_9", "fbx5_3_23", "fbx5_12_3", "S7_B5")) %>%
  filter(context == "CHH") %>%
  droplevels()

# Put WT sample as control
df_fbx5_CHH$genotype <- factor(df_fbx5_CHH$genotype, levels = c("S7_B5", "fbx5_8_9", "fbx5_12_3", "fbx5_3_23"), ordered = TRUE)


ggboxplot(df_fbx5_CHH, x = "genotype", y = "percent_methylation", add = "jitter", title = "CHH methylation in long TEs") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + scale_x_discrete(labels = c("S7-B5", "fbx5-1", "fbx5-2", "fbx5-3"))
```

![](images/mCHH_long_TE_CRISPR_fbx5.png)

#### Statistical test

``` r
lm1 <- aov(percent_methylation ~ genotype, data=df_fbx5_CHH)
summary(glht(lm1, linfct=mcp(genotype="Dunnett"), alternative="two.sided"))
```

     Simultaneous Tests for General Linear Hypotheses

    Multiple Comparisons of Means: Dunnett Contrasts


    Fit: aov(formula = percent_methylation ~ genotype, data = df_fbx5_CHH)

    Linear Hypotheses:
                           Estimate Std. Error t value Pr(>|t|)
    fbx5_8_9 - S7_B5 == 0    -1.963      1.056  -1.860    0.223
    fbx5_12_3 - S7_B5 == 0    0.100      1.056   0.095    0.999
    fbx5_3_23 - S7_B5 == 0   -0.250      1.056  -0.237    0.990
    (Adjusted p values reported -- single-step method)

# RNA-seq library preparation

We generated leaf transcriptomes for 97 Santo Antão accessions,
including 11 accessions with three biological replicates. Rosette leaves
from 4-true leaves (20 DAG) were collected between ZT3 and ZT6 and
flash-frozen into liquid nitrogen. Samples were ground in 2 ml Eppendorf
tubes containing one tungsten carbide ball in a TissueLyser II (Qiagen).
An aliquot of about 20 mg of powder was transferred into a 96 well
plate, and total RNA was prepared with the NucleoMag® 96 RNA kit
(Macherey Nagel). Libraries were prepared with the NEBNext Ultra™
Directional RNA Library Prep Kit for Illumina sequencing (New England
Biolabs). Approximately about 7 million reads of 150 bp single-end reads
were generated on the Illumina sequencer HiSeq3000. The adaptors were
trimmed using Cutadapt (parameters -m 20 -q 35) (Martin 2011). Reads
were mapped on TAIR10 reference genome and Araport11 gene annotation
using HISAT2 (v2.2.0) (Kim et al. 2019, 2). The read count was performed
with HTSeq (v0.12.4) (Anders et al. 2015). Differential expression
analysis was performed in R using the DESeq2 package (v1.28.1) (Love et
al. 2014).

137 samples 98 accessions (97 SA + Col-0) 19 samples in 3 replicates (57
samples) + 1 sample in duplicate (S4-B3-14, 4623_AR and 4623_BF)
CMT2stop = 37/98 FBX5stop = 42/98 VIM2del = 45/98 VIM3alt = 82/98

VIM3alt (1) is defined as the SNP at position Chr5:15047549, which
correspond to the most significant SNP found on Chromosome 5 for gbM
GWAS, and located 790 kb downstream of VIM3 (AT5G39550,
Chr5:15837178..15840678).

| Library Name | seqID   | VIM2_allele | CMT2_allele | FBX5_allele | VIM3_allele |
|--------------|---------|-------------|-------------|-------------|-------------|
| 4559_AA      | 27175   | 0           | 1           | 0           | 1           |
| 4559_AB      | 15675   | 1           | 0           | 1           | 1           |
| 4559_AC      | 15675   | 1           | 0           | 1           | 1           |
| 4559_AD      | 15675   | 1           | 0           | 1           | 1           |
| 4559_AE      | 27180   | 0           | 1           | 0           | 1           |
| 4559_AH      | 22634   | 0           | 0           | 0           | 1           |
| 4559_AI      | 22634   | 0           | 0           | 0           | 1           |
| 4559_AJ      | 22634   | 0           | 0           | 0           | 1           |
| 4559_AK      | 2876_X  | 0           | 1           | 0           | 1           |
| 4559_AL      | 2876_X  | 0           | 1           | 0           | 1           |
| 4559_AM      | 2876_X  | 0           | 1           | 0           | 1           |
| 4559_AN      | 2876_J  | 1           | 0           | 1           | 1           |
| 4559_AO      | 2876_J  | 1           | 0           | 1           | 1           |
| 4559_AP      | 2876_J  | 1           | 0           | 1           | 1           |
| 4559_AQ      | 2876_B  | 1           | 0           | 0           | 1           |
| 4559_AR      | 2876_B  | 1           | 0           | 0           | 1           |
| 4559_AS      | 2876_B  | 1           | 0           | 0           | 1           |
| 4559_AT      | 20683   | 0           | 0           | 0           | 1           |
| 4559_AU      | 20683   | 0           | 0           | 0           | 1           |
| 4559_AV      | 20683   | 0           | 0           | 0           | 1           |
| 4559_AW      | 12912   | 1           | 0           | 0           | 1           |
| 4559_AX      | 12912   | 1           | 0           | 0           | 1           |
| 4559_AY      | 12912   | 1           | 0           | 0           | 1           |
| 4559_AZ      | 13173   | 1           | 0           | 1           | 1           |
| 4559_BA      | 13173   | 1           | 0           | 1           | 1           |
| 4559_BB      | 13173   | 1           | 0           | 1           | 1           |
| 4559_BC      | 13581   | 0           | 1           | 0           | 1           |
| 4559_BD      | 13581   | 0           | 1           | 0           | 1           |
| 4559_BE      | 13581   | 0           | 1           | 0           | 1           |
| 4559_BF      | 22624   | 1           | 0           | 0           | 1           |
| 4559_BG      | 22624   | 1           | 0           | 0           | 1           |
| 4559_BH      | 22624   | 1           | 0           | 0           | 1           |
| 4559_BI      | 35519   | 0           | 1           | 0           | 1           |
| 4559_BJ      | 35519   | 0           | 1           | 0           | 1           |
| 4559_BK      | 35519   | 0           | 1           | 0           | 1           |
| 4559_BL      | 27172   | 0           | 1           | 0           | 1           |
| 4559_BM      | 27172   | 0           | 1           | 0           | 1           |
| 4559_BN      | 27172   | 0           | 1           | 0           | 1           |
| 4559_BO      | 13578   | 0           | 0           | 0           | 1           |
| 4559_BP      | 13578   | 0           | 0           | 0           | 1           |
| 4559_BQ      | 13578   | 0           | 0           | 0           | 1           |
| 4559_BR      | 2876_AD | 0           | 1           | 0           | 1           |
| 4559_BS      | 2876_AD | 0           | 1           | 0           | 1           |
| 4559_BT      | 2876_AD | 0           | 1           | 0           | 1           |
| 4559_BV      | 2876_V  | 0           | 0           | 0           | 1           |
| 4559_C       | 6909    | 0           | 0           | 0           | 0           |
| 4559_D       | 4073_M  | 1           | 0           | 1           | 1           |
| 4559_E       | 4073_M  | 1           | 0           | 1           | 1           |
| 4559_F       | 4073_M  | 1           | 0           | 1           | 1           |
| 4559_J       | 15673   | 0           | 0           | 0           | 1           |
| 4559_K       | 15673   | 0           | 0           | 0           | 1           |
| 4559_L       | 15673   | 0           | 0           | 0           | 1           |
| 4559_M       | 15671   | 1           | 0           | 0           | 1           |
| 4559_N       | 15671   | 1           | 0           | 0           | 1           |
| 4559_O       | 15671   | 1           | 0           | 0           | 1           |
| 4559_P       | 2876_AN | 1           | 0           | 0           | 1           |
| 4559_S       | 16293   | 1           | 0           | 0           | 1           |
| 4559_T       | 16293   | 1           | 0           | 0           | 1           |
| 4559_U       | 16293   | 1           | 0           | 0           | 1           |
| 4559_V       | 27174   | 1           | 1           | 1           | 1           |
| 4559_Y       | 27175   | 0           | 1           | 0           | 1           |
| 4559_Z       | 27175   | 0           | 1           | 0           | 1           |
| 4568_A       | 22637   | 1           | 1           | 1           | 1           |
| 4568_B       | 27158   | 1           | 1           | 0           | 1           |
| 4568_C       | 22615   | 0           | 0           | 0           | 0           |
| 4568_D       | 22620   | 0           | 1           | 1           | 1           |
| 4568_E       | 27170   | 0           | 0           | 1           | 0           |
| 4568_F       | 27160   | 1           | 1           | 0           | 1           |
| 4568_G       | 27161   | 0           | 0           | 1           | 0           |
| 4568_H       | 22631   | 0           | 0           | 1           | 0           |
| 4568_I       | 20682   | 0           | 0           | 1           | 0           |
| 4568_J       | 22622   | 1           | 1           | 0           | 1           |
| 4568_K       | 2876_C  | 0           | 0           | 1           | 0           |
| 4568_L       | 2876_M  | 0           | 0           | 1           | 0           |
| 4568_M       | 2876_R  | 0           | 1           | 1           | 1           |
| 4568_N       | 2876_F  | 0           | 1           | 1           | 1           |
| 4568_O       | 2876_AR | 1           | 0           | 1           | 1           |
| 4568_P       | 22630   | 0           | 1           | 0           | 1           |
| 4568_Q       | 22633   | 1           | 0           | 1           | 1           |
| 4568_R       | 13179   | 1           | 0           | 1           | 1           |
| 4568_S       | 13183   | 1           | 0           | 1           | 1           |
| 4568_T       | 16150   | 0           | 0           | 0           | 1           |
| 4568_U       | 22626   | 0           | 0           | 0           | 0           |
| 4623_A       | 2876_AM | 1           | 1           | 0           | 1           |
| 4623_AB      | 2876_S  | 0           | 1           | 0           | 1           |
| 4623_AC      | 2876_G  | 1           | 0           | 1           | 1           |
| 4623_AD      | 2876_AS | 1           | 0           | 1           | 1           |
| 4623_AE      | 20686   | 1           | 0           | 1           | 1           |
| 4623_AF      | 20685   | 1           | 0           | 1           | 1           |
| 4623_AG      | 2876_K  | 1           | 0           | 1           | 1           |
| 4623_AH      | 2876_Q  | 1           | 0           | 1           | 1           |
| 4623_AI      | 2876_L  | 1           | 0           | 1           | 1           |
| 4623_AJ      | 2876_P  | 1           | 0           | 1           | 1           |
| 4623_AK      | 22643   | 0           | 1           | 0           | 1           |
| 4623_AL      | 27163   | 1           | 1           | 0           | 1           |
| 4623_AM      | 2876_AA | 1           | 0           | 1           | 1           |
| 4623_AN      | 13175   | 1           | 0           | 1           | 1           |
| 4623_AO      | 2876_AT | 1           | 0           | 1           | 1           |
| 4623_AP      | 13172   | 1           | 0           | 1           | 1           |
| 4623_AQ      | 2876_AK | 1           | 0           | 0           | 1           |
| 4623_AR      | 27177   | 1           | 0           | 1           | 1           |
| 4623_AS      | 2876_AU | 1           | 0           | 0           | 1           |
| 4623_AT      | 22632   | 0           | 1           | 1           | 1           |
| 4623_AU      | 22618   | 1           | 0           | 0           | 1           |
| 4623_AW      | 13177   | 1           | 0           | 0           | 1           |
| 4623_AX      | 2876_AV | 0           | 1           | 0           | 1           |
| 4623_AY      | 27165   | 0           | 1           | 0           | 1           |
| 4623_AZ      | 15669   | 0           | 1           | 0           | 1           |
| 4623_B       | 22619   | 1           | 1           | 0           | 1           |
| 4623_BA      | 22617   | 0           | 1           | 0           | 1           |
| 4623_BB      | 2876_Z  | 0           | 1           | 0           | 1           |
| 4623_BC      | 27179   | 0           | 1           | 0           | 1           |
| 4623_BD      | 2876_AW | 0           | 1           | 0           | 1           |
| 4623_BE      | 16292   | 0           | 0           | 0           | 0           |
| 4623_BF      | 27177   | 1           | 0           | 1           | 1           |
| 4623_C       | 15672   | 0           | 0           | 0           | 1           |
| 4623_D       | 22642   | 0           | 1           | 0           | 1           |
| 4623_E       | 22623   | 1           | 1           | 0           | 1           |
| 4623_F       | 2876_AO | 0           | 1           | 0           | 0           |
| 4623_G       | 15674   | 0           | 1           | 0           | 1           |
| 4623_H       | 35517   | 1           | 0           | 0           | 1           |
| 4623_I       | 2876_AP | 0           | 0           | 0           | 1           |
| 4623_J       | 27169   | 0           | 0           | 0           | 1           |
| 4623_K       | 27181   | 1           | 0           | 1           | 1           |
| 4623_L       | 16151   | 1           | 0           | 1           | 1           |
| 4623_N       | 27182   | 0           | 1           | 0           | 1           |
| 4623_P       | 2876_D  | 0           | 0           | 1           | 0           |
| 4623_Q       | 27162   | 0           | 0           | 1           | 0           |
| 4623_R       | 2876_W  | 1           | 0           | 0           | 1           |
| 4623_S       | 22616   | 0           | 1           | 0           | 1           |
| 4623_T       | 22641   | 0           | 1           | 0           | 1           |
| 4623_U       | 20681   | 0           | 0           | 1           | 0           |
| 4623_V       | 22013   | 1           | 0           | 1           | 1           |
| 4623_W       | 22636   | 0           | 0           | 0           | 1           |
| 4623_X       | 2876_O  | 0           | 0           | 1           | 0           |
| 4623_Y       | 2876_N  | 0           | 0           | 1           | 0           |
| 4623_Z       | 2876_AB | 0           | 1           | 1           | 1           |

## Read trimming

Remove reads smaller than 20 nt and nucleotides at 5’ and 3’ ends with
with quality below 35.

``` bash
for i in *.fastq.gz; do
    name=$(basename $i | cut -d. -f1 -)
    if [ ! -e trimmed_reads/${name}.fastq* ]; then
          zcat $i | cutadapt -m 20 -q 35 -o trimmed_reads/${name}.fastq -
      fi
done
```

## Mapping

### Get reference genome

``` bash

# Download fasta file https://www.arabidopsis.org/download/file?path=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz

# Replace chromosome names to add prefix Chr
sed -i 's/^>\([1-5]\)/>Chr\1/g' TAIR10.fa
```

### Get gene annotation Araport11

Download from
<https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release>
the file Araport11_GFF3_genes_transposons.201606.gff.gz

Convert it in gtf with gffread (argument
`-T main output will be GTF instead of GFF3`) (version 0.11.4)

``` bash
# Uncompress
gunzip Araport11_GFF3_genes_transposons.201606.gff.gz

# Convert with gffread 
gffread Araport11_GFF3_genes_transposons.201606.gff -T -o Araport11_GFF3_genes_transposons.201606_gffread.gtf

# Change the names of mitochondria and chloroplast to match the fasta reference

cut -f1 Araport11_GFF3_genes_transposons.201606_gffread.gtf | sort - | uniq
Chr1
Chr2
Chr3
Chr4
Chr5
ChrC
ChrM

sed -i 's/ChrC/chloroplast/g' Araport11_GFF3_genes_transposons.201606_gffread.gtf
sed -i 's/ChrM/mitochondria/g' Araport11_GFF3_genes_transposons.201606_gffread.gtf

cut -f1 Araport11_GFF3_genes_transposons.201606_gffread.gtf | sort - | uniq
Chr1
Chr2
Chr3
Chr4
Chr5
chloroplast
mitochondria
```

Note that the protocol
<https://www.nature.com/articles/nprot.2016.095#Tab1> recommends to
create a file for exon and intron before creating the HISAT index.

The two files generated by the python scripts hisat2_extract_exons.py
and hisat2_extract_splice_sites.py are bed files of the locations of
splicing sites and exons and the corresponding strand orientation.

``` bash
# Get python script
wget https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/hisat2_extract_exons.py
wget https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/hisat2_extract_splice_sites.py


# Extract exon and intron coordinates
python hisat2_extract_exons.py ../Araport11_GFF3_genes_transposons.201606_gffread.gtf > ARAPORT11.exon
python hisat2_extract_splice_sites.py ../Araport11_GFF3_genes_transposons.201606_gffread.gtf > ARAPORT11.ss

# Build the index
hisat2-build --ss ARAPORT11.ss --exon ARAPORT11.exon -f TAIR10.fa TAIR10
```

### Mapping

``` bash
for i in *.fastq.gz; do
  name=$(basename $i | cut -d. -f1)
    hisat2 -p 8 --rna-strandness R --dta -x /path/to/index/TAIR10 -U $i | \
    samtools view -bS -F 4 - | samtools sort - -o ${name}.bam
done
```

## Read counting

### Install HTseq

``` bash
# Install HTseq (version 0.12.4)
pip install --user numpy --upgrade
pip install htseq --upgrade
```

### Gene annotation for htseq

### Keep only CDS and exons

``` bash
# Keep only exon and CDS
grep -vw 'transcript' Araport11_GFF3_genes_transposons.201606_gffread.gtf > Araport11_GFF3_genes.gtf
```

### Remove miRNAs genes

htseq-count bugs because the miRNA have no gene_id field (only
transcript_id).

``` bash
grep -v "ath-miR" Araport11_genes.gtf > Araport11_genes_wo_miRNAs.gtf
```

### Counting

``` bash

mkdir read_count_htseq

for i in *bam; do
    name_file=$(basename $i | cut -d. -f1 -)
    if [! -e ./read_count_htseq/${name_file}.count ]; then
      samtools view $i | \
      htseq-count -s reverse -t exon -i gene_id - Araport11_genes_wo_miRNAs.gtf > ./read_count_htseq/${name_file}.count
  fi
done
```

### Count merging

Gather the read counts of different samples in one file using the bash
script [scripts/merge_counts.sh](merge_counts.sh).

``` bash
cd read_count_htseq

# Create list_file_counts.txt and sample_names.txt file, I would use the library code for sample names
for i in *count; do
  echo $i >> list_file_counts.txt
  echo $i | cut -d. -f1 >> sample_names.txt
done

# Merge
bash merge_counts.sh list_file_counts.txt sample_names.txt > cts.txt
```

The `cts.txt` file contains the raw read count with the samples in
columns and the genes in row.

## Analysis in R

### R Libraries

### Load cts and coldata

``` r
# Load coldata with sample description and variables
coldata <- read.table("data/coldata.txt", header=TRUE, row.names=1, check.names = FALSE)
# Convert the variables into factors
col_names <- names(coldata)
coldata[,col_names] <- lapply(coldata[,col_names] , factor)

# Load cts
# cts <- read.table("data/cts.txt", header=TRUE, row.names=1, check.names = FALSE)

# Turn into a Rds object for smaller size (twice smaller)
# saveRDS(cts,"data/cts.Rds")

# Delete cts.txt file to save space on GitHub
# file.remove("data/cts.txt")

# Load Rds object 
cts <-  readRDS("data/cts.Rds")

### Check that sample names match in both files
all(colnames(cts) == rownames(coldata))
```

### PCA

``` r
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~1
)

# Keep only genes with at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Create PCA
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("sample"), returnData = TRUE)

# Create a population variable by spliting sample
pcaData2 <- pcaData %>%
  separate(group, c("population", "second", "third"), "_") %>%
  select(PC1, PC2, population, sample, name)

# Get percentage variation
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData2, aes(PC1, PC2, color = population, label = sample)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

p + geom_text(size = 4) + ggtitle("PCA RNA-seq (136 samples)") + theme(plot.title = element_text(hjust = 0.5))
```

![](images/PCA_136_samples.png)

We can see that samples cluster partially within populations.

## Replicates analysis (20 accessions)

We have in total 19 accessions with 3 replicates and one accession has 2
replicates (S4-B3-14) so 59 samples.

    Cvi-0
    S10-16
    S10-29
    S11-20
    S11-45
    S11-63
    S15-7
    S15-T2-15-41
    S16-T1-15-47
    S17-1
    S18-2
    S3-12
    S3-4
    S4-B2-9
    S5-100
    S7-B2
    S7-B20
    S7-B5
    S7-T1-15-2
    S4-B3-14

``` r
rep_samples <- as.vector(read.table("data/samples_replicates_RNAseq.txt")$V1)

cts_rep <- cts[,rep_samples]

coldata_rep <- coldata[rep_samples,] %>% droplevels()

# Check if order followed
all(colnames(cts_rep) == rownames(coldata_rep))
```

### PCA analysis

Check distribution of replicates in a PCA

``` r
dds <- DESeqDataSetFromMatrix(
  countData = cts_rep,
  colData = coldata_rep,
  design = ~1
)

# Keep only genes with at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Create PCA
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("sample"), returnData = TRUE)

# Get percentage variation
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = sample, label = sample)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

p + geom_text(size = 4) + ggtitle("PCA RNA-seq replicates") + theme(plot.title = element_text(hjust = 0.5))
```

![](images/PCA_replicates.png)

Most biological replicates cluster tightly together, indicating
homogenous conditions across the growth chamber (replicates were
separated randomly across the trays).

## Analysis without replicates (97 accessions)

Remove Col-0 4559_C and keep first replicate for each accessions with
replicates. There is 97 samples left.

``` r
# Remove Col-0
coldata_filtered <- coldata %>%
  distinct(sample, .keep_all = TRUE) %>%
  filter(sample != "Col_0") %>%
  droplevels()

# I end up with 136 lines

# Filter also cts
cts_filtered <- cts %>% dplyr::select(rownames(coldata_filtered))

# Rename to coldata and cts
coldata <- coldata_filtered
cts <- cts_filtered

rm(coldata_filtered, cts_filtered)
```

### DEG analysis by VIM2 allele

``` r
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~VIM2
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

## Generate expression model
# Get the read count and perform 3 steps: estimation of size factors, estimation of dispersion
# and negative binomial GLM fitting and Wald Statistics
dds <- DESeq(dds)

# Create a DESeqResults object
res <- results(dds)

sum(res$padj < 0.05, na.rm = TRUE)

DEG <- as.data.frame(res) %>% rownames_to_column("geneID")

sigDEG <- as.data.frame(res) %>%
  rownames_to_column("geneID") %>%
  filter(padj < 0.05)

write.table(sigDEG, "data/significant_94_DEGs_VIM2.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```

94 DEGs.

#### Volcano plot

``` r
EnhancedVolcano(res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  axisLabSize = 12,
  xlab = bquote(~ Log[2] ~ Fold ~ Change),
  ylab = bquote(~ -Log[10] ~ P ~ Value),
  pCutoff = 0.05,
  FCcutoff = 2,
  col = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
  labSize = 4,
  title = "",
  subtitle = "DEGs VIM2/4del vs VIM2/4ref",
  legendPosition = "bottom",
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)
```

![](images/volcano_plots_VIM2.png)

[AT5G36100](https://www.arabidopsis.org/locus?key=132582) is the most
significant and is upregulated in VIM2/4del. SLDP2 (SEED LIPID DROPLET
PROTEIN2).Protein of unknown function that is found on the surfaces of
lipid droplets and may function to anchor the droplets to the plasma
membrane.

This gene has two mCG DMRs in its downstream region:

![](images/AT5G36100_snapshot.png)

#### PCA analysis

``` r
# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Create PCA
pcaData <- plotPCA(vsd, intgroup = c("sample", "VIM2"), returnData = TRUE)

# Get percentage variation
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = VIM2, label = sample)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

p + geom_text(size = 3) + ggtitle("PCA RNA-seq 97 samples") + theme(plot.title = element_text(hjust = 0.5))
```

![](images/PCA_97_samples.png)

#### VIM2/VIM4 and VIM3 expression

Combine read counts for *VIM2* and *VIM4* as they are too similar to be
able to differentiate their expression.

``` r
# Get data
# Get data
vim2 <- plotCounts(dds, gene = "AT1G66050", intgroup = "VIM2", main = "VIM2", returnData = T)
vim4 <- plotCounts(dds, gene = "AT1G66040", intgroup = "VIM2", main = "VIM2", returnData = T)
vim3 <- plotCounts(dds, gene = "AT5G39550", intgroup = "VIM2", main = "VIM2", returnData = T)

vim24 <- cbind(vim2, vim4)

# Create a new variable which summarizes the read count
vim24$VIM24_count <- vim24[, 1] + vim24[, 3]

# Remove extra VIM2 column
vim24 <- vim24[, -2]

# Highlight Cvi-0
# Get library names
library_Cvi_0 <- coldata %>%
  filter(sample == "Cvi_0") %>%
  rownames()

# Get index corresponding in vim24 that got Col-0 removed
index_cvi <- which(rownames(vim24) %in% library_Cvi_0)


p1 <- ggboxplot(vim24, y = "VIM24_count", x = "VIM2", add = "jitter", title = "VIM2/4 (AT1G66050/AT1G66040") + theme(plot.title = element_text(hjust = 0.5)) + ylab("Normalized read count") + xlab("VIM2 allele") + geom_point(data = vim24[index_cvi, ], color = "red", pch = 18, cex = 3)

p2 <- ggboxplot(vim3, y = "count", x = "VIM2", add = "jitter", title = "VIM3 (AT5G39550)") + ylab("Normalized read count") + xlab("VIM2 allele") + theme(plot.title = element_text(hjust = 0.5)) + geom_point(data = vim3[index_cvi, ], color = "red", pch = 18, cex = 3)

grid.arrange(p1, p2, nrow = 1)
```

![](images/VIM24_VIM3_expression.png)

#### Statistical test expression VIM2/4

``` r
# Count number of genotype
table(vim24$VIM2)

library(onewaytests)

# Test if variances equals
bartlett.test(VIM24_count ~ VIM2, data = vim24)

# Variance homogenous

# Do t-test (two-tailed)
with(vim24, t.test(VIM24_count[VIM2 == "0"], VIM24_count[VIM2 == "1"], var.equal = TRUE))
```

        Bartlett test of homogeneity of variances

    data:  VIM24_count by VIM2
    Bartlett's K-squared = 0.74395, df = 1, p-value = 0.3884


        Two Sample t-test

    data:  VIM24_count[VIM2 == "0"] and VIM24_count[VIM2 == "1"]
    t = -8.3376, df = 95, p-value = 5.822e-13
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -54.70224 -33.66205
    sample estimates:
    mean of x mean of y 
     36.18342  80.36557 

Highly significant different between the two groups VIM2del and VIM2ref.

#### Statistical test expression VIM3

``` r
# Count number of genotype
table(vim3$VIM2)

library(onewaytests)

# Test if variances equals
bartlett.test(count ~ VIM2, data = vim3)

# Variance homogenous

# Do t-test (two-tailed)
with(vim3, t.test(count[VIM2 == "0"], count[VIM2 == "1"], var.equal = TRUE))
```

        Bartlett test of homogeneity of variances

    data:  count by VIM2
    Bartlett's K-squared = 5.1689, df = 1, p-value = 0.02299


        Two Sample t-test

    data:  count[VIM2 == "0"] and count[VIM2 == "1"]
    t = 3.9866, df = 95, p-value = 0.0001316
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
      4.815942 14.370551
    sample estimates:
    mean of x mean of y 
     37.33505  27.74180 

# Analysis TE expression

## R Libraries

``` r
# RNA-seq analysis
library(DESeq2)

# Plotting functions
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(EnhancedVolcano)

# Perform Bayesian shrinkage estimators for effect sizes in GLM models
#library(apeglm)

# Subset and organize dataframes
library(dplyr)
library(tidyverse)
library(magrittr)

# Handle DESeResults to Dataframe conversion
library(tibble)

# Test for significant overlaps
library(GeneOverlap)

# Make Venn diagram
library(RVenn)

# Heat map
library(pheatmap)

# Function
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}
```

To quantify the expression of transposable elements (TEs) in our RNA-seq
dataset, we used the RepEnrich2 pipeline
(<https://github.com/nerettilab/RepEnrich2>)

## TE annotation

`TAIR10_Transposable_Elements.txt` can be downloaded on
<https://www.arabidopsis.org/download/file?path=Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_Transposable_Elements.txt>

The format for RepEnrich1 should be:

Column 1: Chromosome Column 2: Start Column 3: End Column 4: Repeat_name
Column 5: Class (family) Column 6: Family (superfamily)

``` bash
wc -l TAIR10_Transposable_Elements.txt
35090 TAIR10_Transposable_Elements.txt

# Remove TE genes
sed '/AT[1-5]G/d' TAIR10_Transposable_Elements.txt > TAIR10_Transposable_Elements_only_TEs.txt

# Delete header
sed -i '1d' TAIR10_Transposable_Elements_only_TEs.txt

# Create chromosome column
sed -E 's/^AT([0-9]+)TE[0-9]+/Chr\1/' TAIR10_Transposable_Elements_only_TEs.txt | cut -f1 > chr_names.txt 

paste chr_names.txt TAIR10_Transposable_Elements_only_TEs.txt > TAIR10_Transposable_Elements_only_TEs_with_chr.txt

# Rearrange order column
awk '{print $1,$4,$5,$2,$6,$7}'  OFS='\t'  TAIR10_Transposable_Elements_only_TEs_with_chr.txt > TAIR10_Transposable_Elements_only_TEs_with_chr_reordered.txt

wc -l TAIR10_Transposable_Elements_only_TEs_with_chr_reordered.txt
31189 TAIR10_Transposable_Elements_only_TEs_with_chr_reordered.txt

# Replace slashes by underscore to avoid potential parsing issues
sed -i 's/\//_/g' TAIR10_Transposable_Elements_only_TEs_with_chr_reordered.txt

# Rename
mv TAIR10_Transposable_Elements_only_TEs_with_chr_reordered.txt TAIR10_Transposable_Elements_RepEnrich.txt
```

## Reference genome

Get reference genome from
<https://www.arabidopsis.org/download/list?dir=Genes%2FTAIR10_genome_release>
and build bowtie2 indexes

``` bash
# Make a copy of the fasta file
cp ~/TAIR10_Chr_no_Pt_Mt/TAIR10.fa .

# Check chromosome names
grep ">" TAIR10.fa
>Chr1
>Chr2
>Chr3
>Chr4
>Chr5

# Build a bowtie2 indexes
bowtie2-build TAIR10.fa TAIR10_bowtie2
```

## Bowtie2 index for each TE

For each TE, a fasta file and bowtie2 index files (bt and bt2) are
created. The file `repgenomes_key.txt` contains the name of each TEs and
a unique ID (digit). `repnames.bed` Contains the position of each TE.
This step is performed only once.

``` bash

git clone https://github.com/nerettilab/RepEnrich2.git

# Run
python ./RepEnrich2/RepEnrich2_setup.py --is_bed TRUE TAIR10_Transposable_Elements_RepEnrich.txt TAIR10.fa /srv/netscratch/irg/grp_hancock/johan/repenrich2_pipeline/repenrich2_araport11/
```

## STEP I: Mapping RNA-seq data as first step of RepEnrich2

bowtie2 and bedtools also need to be installed. Used trimmed reads from
the RNA-seq (see RNA-seq part).

``` bash
echo "###############################################################################"
echo "STEP I: Mapping RNA-seq data as first step of RepEnrich2"
echo "###############################################################################"
# bowtie2 in /opt/share/software/bin/bowtie2
# version 2.2.8

fastq_dir="/srv/biodata/irg/grp_hancock/NGS_data/GC_4559/raw_fastq/final/trimmed_reads"
bowtie_index="~/TAIR10_Chr_no_Pt_Mt/bowtie2_index/TAIR10_bowtie2"
output_dir="/srv/netscratch/irg/grp_hancock/johan/repenrich2_pipeline/output_4559"

for i in ${fastq_dir}/*.fastq.gz; do
  name_fastq=$(basename $i | cut -d. -f1)
  # Create a directory to store the output bam file
  if [ ! -d ${output_dir}/${name_fastq} ]; then
    mkdir  ${output_dir}/${name_fastq}
    echo "bowtie2 -q -p 16 -x $bowtie_index -U $i | samtools view -bS - > ${output_dir}/${name_fastq}/${name_fastq}.bam"
    bowtie2 -q -p 16 -x $bowtie_index -U $i | samtools view -bS - > ${output_dir}/${name_fastq}/${name_fastq}.bam
  else
    echo "${output_dir}/${name_fastq} already exists"
  fi
done
```

## STEP II: Split uniquely mapped and multimapping reads

The `RepEnrich2_subset.py` script should be run to output discrete files
for uniquely and multi-mapping reads. MAPQ threshold of 30. Retrieve
uniquely mapped read using the command
`samtools view -h -F 4 -bq 30 file.bam`. The q of -bq selects for
uniquely mapped reads (quality \>=30) as mapping quality of multimapping
reads = 0 if I am correct.

Also, the script selects multimapping reads by first retrieving all
mapped reads (`"samtools view -h -F 4 -b file.bam`). Then
`"samtools view -U _multimap_filtered.bam -bq 30 _map.bam`, which
retrieves all reads not matching -q 30 thanks to parameter -U
`-U FILE  output reads not selected by filters to FILE`

It also converts the \_multimap_filtered.bam file into a fastq file
using command
`bedtools bamtofastq -i _multimap_filtered.bam -fq _multimap.fastq`.

``` bash
echo "###############################################################################"
echo "STEP II: Split uniquely mapped and multimapping reads"
echo "###############################################################################"

repenrich2_index="/srv/netscratch/irg/grp_hancock/johan/repenrich2_pipeline/repenrich2_araport11"
script_repenrich2="/home/zicola/SCRIPTS/analysis_RNA-seq_GC4559_4568_4623/analysis_RepEnrich/RepEnrich2"


# Loop over created directories in step I
for dir in ${output_dir}/*/; do
name_sample=$(basename $dir)
  if [ -e ${dir}${name_sample}_multimap.fastq ] && [ -e ${dir}${name_sample}_unique.bam ] && [ -e ${dir}${name_sample}_multimap_filtered.bam ]; then
    echo "${dir}${name_sample}_unique.bam already exists"
    echo "${dir}${name_sample}_multimap.fastq already exists"
    echo "${dir}${name_sample}_multimap_filtered.bam already exists"
  else
    name_sample=$(basename $dir)
    echo "python ${script_repenrich2}/RepEnrich2_subset.py ${dir}${name_sample}.bam 30 ${dir}${name_sample} --pairedend FALSE"
    python ${script_repenrich2}/RepEnrich2_subset.py ${dir}${name_sample}.bam 30 ${dir}${name_sample} --pairedend FALSE
  fi
done
```

## STEP III: Mapping to TEs

Read count for each TE.

``` bash

echo "###############################################################################"
echo "STEP III: Mapping to TEs"
echo "###############################################################################"

te_annotation="/home/zicola/SCRIPTS/analysis_RNA-seq_GC4559_4568_4623/analysis_RepEnrich/araport11_TE_annotation_fixed.txt"

# When the index will be renamed
for dir in ${output_dir}/*/; do
  name_sample=$(basename $dir)
  if [ ! -e ${dir}${name_sample}_fraction_counts.txt ]; then
    echo "python ${script_repenrich2}/RepEnrich2.py $te_annotation ${dir} ${name_sample} ${repenrich2_index} ${dir}${name_sample}_multimap.fastq ${dir}${name_sample}_unique.bam --is_bed TRUE --cpus 16 --pairedend FALSE"
    python ${script_repenrich2}/RepEnrich2.py $te_annotation ${dir} ${name_sample} ${repenrich2_index} ${dir}${name_sample}_multimap.fastq ${dir}${name_sample}_unique.bam --is_bed TRUE --cpus 16 --pairedend FALSE
  else
    echo "${dir}${name_sample}_fraction_counts.txt already exists"
  fi
done
```

## Summarize data in read count matrix

``` bash
# Get script from my Git
wget https://raw.githubusercontent.com/johanzi/RNA-seq_pipeline/master/merge_counts.sh

# Get raw read count
mkdir read_counts

while read i; do
find output_repenrich -name "${i}_fraction_counts.txt" -exec cp {} read_counts \;
done < sample_names.txt

cd read_counts

ls *txt > list_file_read_counts.txt

# Keep only geneID and read count
while read i; do
cut -f1,4 $i > tmp && mv tmp $i
done < list_file_read_counts.txt

ls *txt > list_file_read_counts.txt
bash ../merge_counts.sh list_file_read_counts.txt ../sample_names.txt > cts_TEs.txt
```

## Analysis in R

``` r
cts <- read.table("data/cts_TEs.txt", header=TRUE, row.names=1, check.names = FALSE)

# Load coldata with sample description and variables
coldata <- read.table("data/coldata.txt", header=TRUE, row.names=1, check.names = FALSE)
# Convert the variables into factors
col_names <- names(coldata)
coldata[,col_names] <- lapply(coldata[,col_names] , factor)

cts <- readRDS("data/cts_TEs.Rds")
# Add population in coldata

## Create a population variable by splitting sample
#coldata_pop <- coldata %>% separate(sample, c("population","second","third"), "_") %>% dplyr::select(population)

#coldata$population <- coldata_pop$population
#coldata$population <- as.factor(coldata$population)

#coldata <- coldata %>% dplyr::select(seqID, sample, population, everything())

# saveRDS(coldata, "data/coldata_TEs.Rds")

coldata <- readRDS("data/coldata_TEs.Rds")
```

### Replicates analysis

``` r
rep_samples <- as.vector(read.table("data/samples_replicates_RNAseq.txt")$V1)

cts_rep <- cts[,rep_samples]

coldata_rep <- coldata[rep_samples,]

# Check if order followed
all(colnames(cts_rep) == rownames(coldata_rep))
```

Check distribution of replicates in a PCA

``` r
dds <- DESeqDataSetFromMatrix(
  countData = cts_rep,
  colData = coldata_rep,
  design = ~1
)

# Keep only genes with at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# Variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# Create PCA
pcaData <- DESeq2::plotPCA(vsd, intgroup=c("sample"), returnData=TRUE)

# Get percentage variation
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = sample, label = sample)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

p + geom_text(size = 4) + ggtitle("PCA TE expression replicates") + theme(plot.title = element_text(hjust = 0.5))
```

![](images/PCA_replicates_TE_expression.png)

Most replicates cluster together for the different accessions but it
less clustered than for gene expression. Remove replicates to keep only
one sample per accession. Remove Col-0 from the analysis. We end up with
97 samples.

### Analysis on single accessions and long TEs

Keep only one replicate per accession.

``` r
coldata_filtered <- coldata %>%
  distinct(sample, .keep_all = TRUE) %>%
  filter(sample != "Col_0")

# I end up with 97 lines

# Filter also cts
cts_filtered <- cts %>% dplyr::select(rownames(coldata_filtered))


# Rename to coldata and cts
coldata <- coldata_filtered
cts <- cts_filtered

rm(coldata_filtered, cts_filtered)
```

Focus on long TEs (\> 4kb)

``` r
# TE_classification <- read.delim("data/TAIR10_Transposable_Elements_RepEnrich.txt", header=FALSE)

# Rename
# colnames(TE_classification) <- c("chr","start","end","name","family","superfamily")

# Change variable
# TE_classification <- TE_classification %>% mutate_at(c(4,5,6), as_factor)

# saveRDS(TE_classification, "data/TAIR10_Transposable_Elements_RepEnrich.Rds")

TE_classification <- readRDS("data/TAIR10_Transposable_Elements_RepEnrich.Rds")

# Get TEs longer than 4 kb
TE_classification_long <- TE_classification %>% filter(end - start > 3999)
# 1235 TEs

# Subset cts with only long TEs
cts_long_TEs <- cts[TE_classification_long$name, ]

idx <- match(TE_classification_long$name, rownames(cts))

cts_long_TEs <- cts[idx, ]

# Check if names match
all(rownames(cts_long_TEs) %in% TE_classification_long$name)

saveRDS(cts_long_TEs, "data/cts_long_TEs.Rds")

cts_long_TEs <- readRDS("data/cts_long_TEs.Rds")
```

### Analysis by CMT2 allele

``` r
dds <- DESeqDataSetFromMatrix(
  countData = cts_long_TEs,
  colData = coldata,
  design = ~CMT2
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# Create a DESeqResults object
res <- results(dds)

sum(res$padj < 0.05, na.rm = TRUE)
# 38

DEG <- as.data.frame(res) %>% rownames_to_column("geneID")

sigDEG <- as.data.frame(res) %>%
  rownames_to_column("geneID") %>%
  filter(padj < 0.05)

sum(sigDEG$log2FoldChange < 0)
# 23
sum(sigDEG$log2FoldChange > 0)
# 15

sigDEG_annot <- merge.data.frame(sigDEG, TE_classification, by.x = "geneID", by.y = "name")

sigDEG_annot <- merge.data.frame(sigDEG, TE_classification, by.x = "geneID", by.y = "name")

write.table(sigDEG_annot, "data/significant_38_TEs_CMT2.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```

#### Volcano plots

``` r
require(EnhancedVolcano)

TE_subset <- TE_classification[(TE_classification$name %in% dds@rowRanges@partitioning@NAMES), ]

p1 <- EnhancedVolcano(res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  axisLabSize = 12,
  xlab = bquote(~ Log[2] ~ Fold ~ Change),
  ylab = bquote(~ -Log[10] ~ P ~ Value),
  pCutoff = 0.05,
  FCcutoff = 2,
  col = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
  labSize = 4,
  title = "",
  subtitle = "TE name",
  legendPosition = "bottom",
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)

# Provide now a new label
p2 <- EnhancedVolcano(res,
  lab = TE_subset$family,
  x = "log2FoldChange",
  y = "padj",
  axisLabSize = 12,
  xlab = bquote(~ Log[2] ~ Fold ~ Change),
  ylab = bquote(~ -Log[10] ~ P ~ Value),
  pCutoff = 0.05,
  FCcutoff = 2,
  col = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
  labSize = 4,
  title = "",
  subtitle = "TE family",
  legendPosition = "bottom",
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)

p3 <- EnhancedVolcano(res,
  lab = TE_subset$superfamily,
  x = "log2FoldChange",
  y = "padj",
  axisLabSize = 12,
  xlab = bquote(~ Log[2] ~ Fold ~ Change),
  ylab = bquote(~ -Log[10] ~ P ~ Value),
  pCutoff = 0.05,
  FCcutoff = 2,
  col = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
  labSize = 4,
  title = "",
  subtitle = "TE superfamily",
  legendPosition = "bottom",
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)

grid.arrange(p1, p2, p3, nrow = 1)
```

![](images/volcano_plots_TEs_CMT2.png)

### Analysis by FBX5 allele

``` r
dds <- DESeqDataSetFromMatrix(
  countData = cts_long_TEs,
  colData = coldata,
  design = ~FBX5
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# and negative binomial GLM fitting and Wald Statistics
dds <- DESeq(dds)

# Create a DESeqResults object
res <- results(dds)

sum(res$padj < 0.05, na.rm = TRUE)
# 30

DEG <- as.data.frame(res) %>% rownames_to_column("geneID")

sigDEG <- as.data.frame(res) %>%
  rownames_to_column("geneID") %>%
  filter(padj < 0.05)

sum(sigDEG$log2FoldChange < 0)
# 16
sum(sigDEG$log2FoldChange > 0)
# 14

sigDEG_annot <- merge.data.frame(sigDEG, TE_classification, by.x = "geneID", by.y = "name")

write.table(sigDEG_annot, "data/significant_30_TEs_FBX5.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```

#### Volcano plots

``` r
require(EnhancedVolcano)

TE_subset <- TE_classification[(TE_classification$name %in% dds@rowRanges@partitioning@NAMES), ]

p1 <- EnhancedVolcano(res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  axisLabSize = 12,
  xlab = bquote(~ Log[2] ~ Fold ~ Change),
  ylab = bquote(~ -Log[10] ~ P ~ Value),
  pCutoff = 0.05,
  FCcutoff = 2,
  col = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
  labSize = 4,
  title = "",
  subtitle = "TE name",
  legendPosition = "bottom",
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)

# Provide now a new label
p2 <- EnhancedVolcano(res,
  lab = TE_subset$family,
  x = "log2FoldChange",
  y = "padj",
  axisLabSize = 12,
  xlab = bquote(~ Log[2] ~ Fold ~ Change),
  ylab = bquote(~ -Log[10] ~ P ~ Value),
  pCutoff = 0.05,
  FCcutoff = 2,
  col = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
  labSize = 4,
  title = "",
  subtitle = "TE family",
  legendPosition = "bottom",
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)

p3 <- EnhancedVolcano(res,
  lab = TE_subset$superfamily,
  x = "log2FoldChange",
  y = "padj",
  axisLabSize = 12,
  xlab = bquote(~ Log[2] ~ Fold ~ Change),
  ylab = bquote(~ -Log[10] ~ P ~ Value),
  pCutoff = 0.05,
  FCcutoff = 2,
  col = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
  labSize = 4,
  title = "",
  subtitle = "TE superfamily",
  legendPosition = "bottom",
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.75
)

grid.arrange(p1, p2, p3, nrow = 1)
```

![](images/volcano_plots_TEs_FBX5.png)

### Permutation

Check how many TEs are differentially expressed when permutating alleles
within populations.

``` r
permutate_deseq2 <- function(allele_name) {
  coldata_shuffled <- coldata %>%
    group_by(population) %>%
    mutate(RAND = sample(get(allele_name)))

  dds <- DESeqDataSetFromMatrix(
    countData = cts_long_TEs,
    colData = coldata_shuffled,
    design = ~RAND
  )

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]

  dds <- DESeq(dds, quiet = TRUE)

  # Create a DESeqResults object
  res <- results(dds)

  sigDEG <- as.data.frame(res) %>%
    rownames_to_column("geneID") %>%
    filter(padj < 0.05)

  # nb_TEs_overexpressed <- sum(sigDEG$log2FoldChange>0)

  # return(nb_TEs_overexpressed)
  return(sigDEG)
}
```

``` r
permutation_TE_CMT2 <- vector(mode="list",length=100)
for(i in 1:100){permutation_TE_CMT2[[i]] <- permutate_deseq2("CMT2")}
saveRDS(permutation_TE_CMT2, "data/permutation_TE_CMT2.Rds")

permutation_TE_FBX5 <- vector(mode="list",length=100)
for(i in 1:100){permutation_TE_FBX5[[i]] <- permutate_deseq2("FBX5")}
saveRDS(permutation_TE_FBX5, "data/permutation_TE_FBX5.Rds")
```

``` r
plot_permutation <- function(vector_permutation, observed_value, ylab) {
  observed_value <- as.integer(observed_value)
  vector_permutation <- as.vector(vector_permutation)

  df <- as.data.frame(matrix(nrow = length(vector_permutation), ncol = 1))

  names(df) <- "value"

  df$value <- as.integer(vector_permutation)

  ggplot(df, aes(x = "", y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter()) +
    theme_bw() +
    ylab(ylab) +
    xlab("") +
    geom_hline(yintercept = observed_value, col = "red") +
    theme(
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    )
}
```

``` r
permutation_CMT2 <- readRDS("data/permutation_TE_CMT2.Rds")

sig_TE_CMT2 <- read.table("data/significant_38_TEs_CMT2.txt", header = T)
sum(sig_TE_CMT2$log2FoldChange > 2)
# 9

nb_up_CMT2 <- as.vector(lapply(permutation_CMT2, function(x) sum(x$log2FoldChange > 2)))

p1 <- plot_permutation(nb_up_CMT2, 9, ylab = "# TEs log2FC > 2") + ggtitle("Permutation CMT2") + theme(plot.title = element_text(hjust = 0.5))

permutation_FBX5 <- readRDS("data/permutation_TE_FBX5.Rds")

sig_TE_FBX5 <- read.table("data/significant_30_TEs_FBX5.txt", header = T)
sum(sig_TE_FBX5$log2FoldChange < -2)
# 3

nb_down_FBX5 <- as.vector(lapply(permutation_FBX5, function(x) sum(x$log2FoldChange < -2)))

p2 <- plot_permutation(nb_down_FBX5, 3, ylab = "# TEs log2FC < -2") + ggtitle("Permutation FBX5") + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2, nrow = 1)
```

![](images/permutation_TEs.png)

Observed values (horizontal red lines) are above the permutation
distribution.

# TE transposition using read coverage as proxy

Map reads to the TAIR10 TE annotation allowing two mismatches (default
bowtie) and look at the reads per millions for each accessions

``` bash
# Raw reads
/srv/netscratch/irg/grp_hancock/raw_CVI_MAD_fastq

 ~/bin/seqkit stats 4073_L_R1.fastq.gz
file                format  type    num_seqs        sum_len  min_len  avg_len  max_len
4073_L_R1.fastq.gz  FASTQ   DNA   10,566,430  1,590,859,487       15    150.6      151

# Data have been trimmed already

# TE annotation (31189)

~/bin/seqkit stats /home/zicola/TAIR10_annotations/TAIR10_TE.fas
file           format  type  num_seqs     sum_len  min_len  avg_len  max_len
TAIR10_TE.fas  FASTA   DNA     31,189  23,315,940       10    747.6   31,019

/home/zicola/TAIR10_annotations/TAIR10_Transposable_Elements_EvaFormat.bed

# Get only overlapping
cd /srv/netscratch/irg/grp_hancock/johan/TE_mapping

cp /home/zicola/TAIR10_annotations/TAIR10_TE.fas .

bowtie-build -f TAIR10_TE.fas TAIR10_TE

index="/srv/netscratch/irg/grp_hancock/johan/TE_mapping/TAIR10_TE"

while read i; do
if [ -e /srv/netscratch/irg/grp_hancock/raw_CVI_MAD_fastq/${i}_R1.fastq.gz ]; then
  echo "$i present"
else
  echo "$i ABSENT!!!!!!!!!!!!!!!!!!"
fi
done < list_accessions.txt

# OK, all present

path_fastq="/srv/netscratch/irg/grp_hancock/raw_CVI_MAD_fastq"
index="/srv/netscratch/irg/grp_hancock/johan/TE_mapping/TAIR10_TE"

while read i; do
  bowtie --threads 6 $index -S -1 ${path_fastq}/${i}_R1.fastq.gz -2 ${path_fastq}/${i}_R2.fastq.gz  |  samtools sort | samtools view -bS -F4 -o mapped_TEs/${i}.bam
done < test

# reads processed: 22195073
# reads with at least one reported alignment: 58166 (0.26%)
# reads that failed to align: 22136907 (99.74%)


# Many error message "Warning: Exhausted best-first chunk memory for read ..."
# Seem to be due to limited memoru
# https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/13897-bowtie-memory-warning
# Recommend using option --chunkmbs 200

# http://dell-head.mpipz.mpg.de/ganglia/?c=DC_Cluster&h=clusternode2.cluster.local&m=load_one&r=hour&s=by%20name&hc=4&mc=2
# Indicates there is plenty of mem on dell-node-2 (up to 280 Gb)

while read i; do
  if [ -e mapped_TEs/${i}.bam ]; then
    echo "mapped_TEs/${i}.bam already exists"
  else
    bowtie --threads 6 --chunkmbs 200  $index -S -1 ${path_fastq}/${i}_R1.fastq.gz -2 ${path_fastq}/${i}_R2.fastq.gz  |  samtools sort | samtools     view -bS -F4 -o mapped_TEs/${i}.bam
  fi
done < list_accessions.txt
```

Summarize mapping statistics

``` bash

cd /srv/netscratch/irg/grp_hancock/johan/TE_mapping/mapped_TEs

for i in *bam; do
  reads=$(samtools flagstat $i | head -n1 | cut -d' ' -f1)
  echo -e "${i%%.*}\t$reads" >> summary_TE_mapping.txt
done

sort -k1 summary_TE_mapping.txt > summary_TE_mapping.sorted.txt

# Get total reads

cd /srv/netscratch/irg/grp_hancock/johan/TE_mapping
cd /srv/netscratch/irg/grp_hancock/raw_CVI_MAD_fastq

while read i; do
  seqkit stats -T ${i}_R1.fastq.gz >> summary.txt
done < list_accessions.txt

grep -v "file" summary.txt | cut -f1, 4 > summary_read_count.txt

sed -i 's/_R1.fastq.gz//g' summary_read_count.txt

sort -k1 summary_read_count.txt > summary_read_count.sorted.txt

paste summary_read_count.sorted.txt summary_TE_mapping.sorted.txt | cut -f1,2,4 > summary_read_count_and_TEs.sorted.txt
```

Combine total read with TE-mapped read count

``` bash
cut -f2 summary_TE_mapping.txt > TE_readcount
paste summary_read_count.txt TE_readcount > summary_read_count_total_and_TEs.txt
```

``` r
allele_id <- read.delim("data/tepid/allele_id.txt")

df_TE_count <- read.delim("data/summary_read_count_and_TEs.sorted.txt", header=FALSE)

names(df_TE_count) <- c("seqID", "total_reads","TE_reads")
df_TE_count$seqID <- as.factor(df_TE_count$seqID)

df_TE_count$ratio_TE <- df_TE_count$TE_reads/df_TE_count$total_reads

df_TE_count_allele <-  merge(df_TE_count, allele_id, by="seqID")

df_TE_count_allele$FBX5 <- as.factor(df_TE_count_allele$FBX5)
df_TE_count_allele$CMT2 <- as.factor(df_TE_count_allele$CMT2)

# Plot nb of reads
TE_plot_per_allele <- function(df, type, group){
  
  require(ggplot2)
  require(tidyr)
  
  # Remove potential rows with NA for genotype (one case for FBX5)
  df_clean <- df %>% drop_na(group)
  
  give.n <- function(x){
    return(c(y = mean(x), label = length(x)))
  }
  
  ggplot(data=df_clean, aes_string(x=group, y=type, group=group))+
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(data=df_clean[-190, ], size=1, width=0.2) +
    geom_point(data=df_clean[190,], size=2, fill="red", shape=23, width=0.2) +
    ggtitle(paste(type," for ",group, sep="")) + 
    theme_bw() + ylab("% reads mapping to TEs") + xlab("") +
    stat_summary(fun.data = give.n, geom = "text") +
    theme(axis.text.x = element_text(color="black"), 
      axis.text.y = element_text(color="black"),
      axis.ticks = element_line(color = "black")) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      scale_y_continuous(labels = scales::percent_format())
}
```

### CMT2

``` r
TE_plot_per_allele(df_TE_count_allele, "ratio_TE","CMT2")
```

![](images/TE_mapping_WGS_CMT2.png)

``` r
# Libraries for Welch test
require(onewaytests)

# Normality
lm <- aov(ratio_TE~CMT2, data=df_TE_count_allele)
lm.stdres <- rstandard(lm)
shapiro.test(lm.stdres)
# Nope

# Homoscedasticity
bartlett.test(ratio_TE~CMT2, df_TE_count_allele)
# Yep

wilcox.test(df_TE_count_allele$ratio_TE[df_TE_count_allele$CMT2==0], df_TE_count_allele$ratio_TE[df_TE_count_allele$CMT2==1], alternative = "two.sided")
```

        Shapiro-Wilk normality test

    data:  lm.stdres
    W = 0.88789, p-value = 9.953e-11


        Bartlett test of homogeneity of variances

    data:  ratio_TE by CMT2
    Bartlett's K-squared = 2.1333, df = 1, p-value = 0.1441


        Wilcoxon rank sum test with continuity correction

    data:  df_TE_count_allele$ratio_TE[df_TE_count_allele$CMT2 == 0] and df_TE_count_allele$ratio_TE[df_TE_count_allele$CMT2 == 1]
    W = 2894, p-value = 0.004321
    alternative hypothesis: true location shift is not equal to 0

### FBX5

``` r
TE_plot_per_allele(df_TE_count_allele, "ratio_TE","FBX5")
```

![](images/TE_mapping_WGS_FBX5.png)

``` r
# Libraries for Welch test
require(onewaytests)

# Normality
lm <- aov(ratio_TE~FBX5, data=df_TE_count_allele)
lm.stdres <- rstandard(lm)
shapiro.test(lm.stdres)
# Nope

# Homoscedasticity
bartlett.test(ratio_TE~FBX5, df_TE_count_allele)
# Yep

wilcox.test(df_TE_count_allele$ratio_TE[df_TE_count_allele$FBX5==0], df_TE_count_allele$ratio_TE[df_TE_count_allele$FBX5==1], alternative = "two.sided")
```

        Shapiro-Wilk normality test

    data:  lm.stdres
    W = 0.8609, p-value = 3.513e-12


        Bartlett test of homogeneity of variances

    data:  ratio_TE by FBX5
    Bartlett's K-squared = 0.063316, df = 1, p-value = 0.8013


        Wilcoxon rank sum test with continuity correction

    data:  df_TE_count_allele$ratio_TE[df_TE_count_allele$FBX5 == 0] and df_TE_count_allele$ratio_TE[df_TE_count_allele$FBX5 == 1]
    W = 4142, p-value = 0.6254
    alternative hypothesis: true location shift is not equal to 0

# Marginal genealogical tree with RELATE

Install the RELATE software from <https://myersgroup.github.io/relate/>
(see paper <https://www.nature.com/articles/s41588-019-0484-x>).

## VIM2

``` bash
RELATE="/path_to_RELATE/Relate"

bp=24586731

## Relate
${RELATE}/bin/Relate \
--mode All \
--haps chr1_SantoAntao_haploid.haps \
--sample chr1_SantoAntao_haploid.sample \
--map chr1_recmap.map \
--memory 20 --coal SA_relate_popsize.coal -m 2.2131e-09 \
-o  chr1_relate

# SampleBranchLengths (--format a)
${RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
-i chr1_relate \
-o chr1_relate_resample200 \
-m 2.2131e-09 \
--coal SA_relate_popsize.coal \
--format a \
--num_samples 200 \
--first_bp ${bp} \
--last_bp ${bp} \
--seed 1 

# SampleBranchLengths (CLUES)
${RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
-i chr1_relate \
-o chr1_relate_resample200_CLUES \
-m 2.2131e-09 \
--coal SA_relate_popsize.coal \
--format b \
--num_samples 200 \
--first_bp ${bp} \
--last_bp ${bp} \
--seed 1 

## TreeViewMutation
${RELATE}/scripts/TreeView/TreeViewMutation.sh \
--haps chr1_SantoAntao_haploid.haps \
--sample chr1_SantoAntao_haploid.sample \
--anc chr1_relate.anc \
--mut chr1_relate.mut \
--poplabels SantoAntao_poplabels.txt \
--bp_of_interest ${bp} \
--years_per_gen 1 \
-o treeview_mutation_VIM2

## TreeViewSample
${RELATE}/scripts/TreeView/TreeViewSample.sh \
--haps chr1_SantoAntao_haploid.haps \
--sample chr1_SantoAntao_haploid.sample \
--anc chr1_relate_resample200.anc \
--mut chr1_relate_resample200.mut \
--dist chr1_relate_resample200.dist \
--poplabels SantoAntao_poplabels.txt \
--bp_of_interest ${bp} \
--years_per_gen 1 \
-o treeview_sample_VIM2
```

## CMT2

``` bash
RELATE="/path_to_RELATE/Relate"

bp=10707974

## Relate
${RELATE}/bin/Relate \
--mode All \
--haps chr4_SantoAntao_haploid.haps \
--sample chr4_SantoAntao_haploid.sample \
--map chr4_recmap.map \
--memory 20 --coal SA_relate_popsize.coal -m 2.142e-09 \
-o  chr4_relate

# SampleBranchLengths (--format a)
${RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
-i chr4_relate \
-o chr4_relate_resample200 \
-m 2.142e-09 \
--coal SA_relate_popsize.coal \
--format a \
--num_samples 200 \
--first_bp ${bp} \
--last_bp ${bp} \
--seed 1 

# SampleBranchLengths (CLUES)
${RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
-i chr4_relate \
-o chr4_relate_resample200_CLUES \
-m 2.142e-09 \
--coal SA_relate_popsize.coal \
--format b \
--num_samples 200 \
--first_bp ${bp} \
--last_bp ${bp} \
--seed 1 

## TreeViewMutation
${RELATE}/scripts/TreeView/TreeViewMutation.sh \
--haps chr4_SantoAntao_haploid.haps \
--sample chr4_SantoAntao_haploid.sample \
--anc chr4_relate.anc \
--mut chr4_relate.mut \
--poplabels SantoAntao_poplabels.txt \
--bp_of_interest ${bp} \
--years_per_gen 1 \
-o treeview_mutation_CMT2

## TreeViewSample
${RELATE}/scripts/TreeView/TreeViewSample.sh \
--haps chr4_SantoAntao_haploid.haps \
--sample chr4_SantoAntao_haploid.sample \
--anc chr4_relate_resample200.anc \
--mut chr4_relate_resample200.mut \
--dist chr4_relate_resample200.dist \
--poplabels SantoAntao_poplabels.txt \
--bp_of_interest ${bp} \
--years_per_gen 1 \
-o treeview_sample_CMT2
```

## FBX5

``` bash
RELATE="/path_to_RELATE/Relate"

bp=18513626

## Relate
${RELATE}/bin/Relate \
--mode All \
--haps chr2_SantoAntao_haploid.haps \
--sample chr2_SantoAntao_haploid.sample \
--map chr2_recmap.map \
--memory 20 --coal SA_relate_popsize.coal -m 2.142e-09 \
-o  chr2_relate

# SampleBranchLengths (--format a)
${RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
-i chr2_relate \
-o chr2_relate_resample200 \
-m 2.142e-09 \
--coal SA_relate_popsize.coal \
--format a \
--num_samples 200 \
--first_bp ${bp} \
--last_bp ${bp} \
--seed 1 

# SampleBranchLengths (CLUES)
${RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
-i chr2_relate \
-o chr2_relate_resample200_CLUES \
-m 2.142e-09 \
--coal SA_relate_popsize.coal \
--format b \
--num_samples 200 \
--first_bp ${bp} \
--last_bp ${bp} \
--seed 1 

## TreeViewMutation
${RELATE}/scripts/TreeView/TreeViewMutation.sh \
--haps chr2_SantoAntao_haploid.haps \
--sample chr2_SantoAntao_haploid.sample \
--anc chr2_relate.anc \
--mut chr2_relate.mut \
--poplabels SantoAntao_poplabels.txt \
--bp_of_interest ${bp} \
--years_per_gen 1 \
-o treeview_mutation_FBX5

## TreeViewSample
${RELATE}/scripts/TreeView/TreeViewSample.sh \
--haps chr2_SantoAntao_haploid.haps \
--sample chr2_SantoAntao_haploid.sample \
--anc chr2_relate_resample200.anc \
--mut chr2_relate_resample200.mut \
--dist chr2_relate_resample200.dist \
--poplabels SantoAntao_poplabels.txt \
--bp_of_interest ${bp} \
--years_per_gen 1 \
-o treeview_sample_FBX5
```

# Inference of selection coefficient

Install software Clues on GitHub
<https://github.com/standard-aaron/clues> (see paper
<https://doi.org/10.1371/journal.pgen.1008384>).

## VIM2

``` bash
clues_inference="/path_to_clues/clues"
timeBins="/path_to_timebins/timebins_1epoch.txt"

## Inference of selection and allele frequency trajectory for VIM2 (chr1:24586731)
python ${clues_inference}/inference.py \
--times chr1_relate_resample200_CLUES \
--popFreq 0.468 \
--tCutoff 5000 \
--timeBins ${timeBins} \
--coal SA_relate_popsize.coal \
--sMax 1 \
--df 100 \
--dom 0 \
--out chr1_VIM2_inference
```

## CMT2

``` bash
clues_inference="/path_to_clues/clues"
timeBins="/path_to_timebins/timebins_1epoch.txt"

## Inference of selection and allele frequency trajectory for CMT2 (Chr4:10707974)
python ${clues_inference}/inference.py \
--times chr4_relate_resample200_CLUES \
--popFreq 0.33 \
--tCutoff 5000 \
--timeBins ${timeBins} \
--coal SA_relate_popsize.coal \
--sMax 1 \
--df 100 \
--dom 0 \
--out chr4_CMT2_inference
```

## FBX5

``` bash
clues_inference="/path_to_clues/clues"
timeBins="/path_to_timebins/timebins_1epoch.txt"

## Inference of selection and allele frequency trajectory for FBX5 (chr2:18513626)
python ${clues_inference}/inference.py \
--times chr2_relate_resample200_CLUES \
--popFreq 0.33 \
--tCutoff 5000 \
--timeBins ${timeBins} \
--coal SA_relate_popsize.coal \
--sMax 1 \
--df 100 \
--dom 0 \
--out chr2_FBX5_inference
```

# Selective sweep analysis (figure 7)

Scripts from Ahmed Elfarargi.

``` bash
# 1. Selective sweep analysis (RAiSD)
  # bgzip and tabix
  bgzip -c subset_189_SA_accessions_biallelic_DP3_GQ25.recode.vcf > subset_189_SA_accessions_biallelic_DP3_GQ25.recode.vcf.gz
  tabix -p vcf subset_189_SA_accessions_biallelic_DP3_GQ25.recode.vcf.gz
  # Split by Chr
  for i in {1,2,3,4,5}; do bcftools view -r Chr${i} subset_189_SA_accessions_biallelic_DP3_GQ25.recode.vcf.gz > subset_189_SA_accessions_biallelic_DP3_GQ25_Ch${i}.vcf; done
  # Run RAiSD
  for i in {1..5}; do RAiSD -n Santo_RAiSD_w50_chr${i} -I subset_189_SA_accessions_biallelic_DP3_GQ25_Ch${i}.vcf -y 1 -A 0.905 -M 0 -w 50 -a 123 -D -R -P -O -s -f; done

# 2. Pairwise Fst analysis using vcftools
  vcftools \
  --vcf subset_189_SA_accessions_biallelic_DP3_GQ25.recode.vcf \
  --weir-fst-pop picoespong.txt \
  --weir-fst-pop covafig.txt \
  --fst-window-size 10000 \
  --fst-window-step 1000 \
  --out picoespong_vs_covafig_w10k_s1k
```

``` r
# libraries
library(tidyverse)
library(sf)
library(maptiles)
library(tidyterra)
library(terra)
library(ggrepel)
library(ggspatial)
library(RColorBrewer)
library(ggtext)
library(glue)
library(ggpubr) 
library(cowplot)

# ==============================================================================
# load data 
# ==============================================================================
file_coords <- "data/figure7/santo_coords_inds.csv"
file_pheno  <- "data/figure7/santo_pheno_data.csv"

# order for variants
target_order <- c("VIM2/VIM4", "CMT2", "FBX5")

# RAiSD Files (Panel C)
files_raisd <- c(
  "Chr1" = "data/figure7/Santo_RAiSD_w50_chr1.Chr1",
  "Chr2" = "data/figure7/Santo_RAiSD_w50_chr2.Chr2",
  "Chr3" = "data/figure7/Santo_RAiSD_w50_chr3.Chr3",
  "Chr4" = "data/figure7/Santo_RAiSD_w50_chr4.Chr4",
  "Chr5" = "data/figure7/Santo_RAiSD_w50_chr5.Chr5"
)

# Fst File (Panel D)
file_fst <- "data/figure7/picoespong_vs_covafig_w10k_s1k.windowed.weir.fst"

# ==============================================================================
# Panel A: Geographical map
# ==============================================================================
data_a <- read.table("data/figure7/Santo_coords_groups.csv", header = TRUE, sep = ",") %>%
  rename(Group = Region, Cluster = Group, n = N, lon = Longitude, lat = Latitude) %>%
  arrange(Cluster, Group) %>%
  mutate(Legend_Label = factor(paste0(Cluster, ": ", Group), levels = unique(paste0(Cluster, ": ", Group))))

data_a$Group <- factor(data_a$Group, levels = unique(data_a$Group))

points_sf <- st_as_sf(data_a, coords = c("lon", "lat"), crs = 4326)

base_cols <- c("#0075DC", "yellow4", "#C20088", "#2BCE48")
final_colors <- c()
clusters <- unique(data_a$Cluster)

for(i in seq_along(clusters)) {
  grps <- data_a %>% filter(Cluster == clusters[i]) %>% pull(Legend_Label) %>% unique() %>% sort()
  pal <- colorRampPalette(c("white", base_cols[i], "black"))(length(grps) + 4)[3:(length(grps) + 2)]
  names(pal) <- grps
  final_colors <- c(final_colors, pal)
}

bbox <- st_bbox(c(xmin = -25.1, xmax = -25.01, ymin = 17.09, ymax = 17.125), crs = st_crs(4326))
tiles <- get_tiles(st_as_sfc(bbox), provider = "Esri.WorldImagery", zoom = 15, crop = TRUE)

p_a <- ggplot() +
  geom_spatraster_rgb(data = tiles, alpha = 0.6) +
  geom_sf(data = points_sf, aes(color = Legend_Label, size = n), show.legend = "point") +
  geom_label_repel(data = data_a, aes(lon, lat, label = Group, fill = Legend_Label),
                   fontface = "bold", size = 4, box.padding = 0.5, color = "white", alpha = 0.8,
                   min.segment.length = 0, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = final_colors) +
  scale_fill_manual(values = final_colors) +
  scale_size_area(max_size = 10, name = NULL) +
  guides(color = "none", fill = "none", size = guide_legend(title.position = "top")) +
  annotate("text", x=-25.085, y=17.12, label="Espongeiro", color="black", size = 6, fontface="bold") +
  annotate("text", x=-25.060, y=17.10, label="Cova", color="black", size = 6, fontface="bold") +
  annotate("text", x=-25.040, y=17.098, label="Figueira", color="black", size = 6, fontface="bold") +
  annotate("text", x=-25.02, y=17.1075, label="Pico", color="black", size = 6, fontface="bold") +
  annotation_scale(location = "br", width_hint = 0.2, style = "ticks", line_col = "black", text_col = "black") +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(line_col = "black", text_col = "black"), height = unit(1, "cm")) +
  coord_sf(xlim = c(-25.1, -25.01), ylim = c(17.09, 17.125), expand = FALSE) +
  labs(tag = "a") +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), 
    axis.title = element_blank(),
    plot.tag = element_text(size = 20, face = "bold"),
    legend.position = c(0.999, 0.999), 
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "black", linewidth = 0.5),
    legend.key = element_blank(),
    legend.title = element_text(face = "bold", size = 10, hjust = 0.5),
    legend.text = element_text(size = 10)
  )

# ==============================================================================
# Panel B: Derived allele frequency
# ==============================================================================
cols_c <- c("Cova"="#0075DC", "Espongeiro"="yellow4", "Figueira"="#C20088", "Pico"="#2BCE48")
cols_g <- c("CMT2"="navyblue", "FBX5"="coral2", "VIM2/VIM4"="palegreen4")
ord_c <- c("Espongeiro", "Cova", "Figueira", "Pico")

df_b <- read.csv(file_coords) %>%
  pivot_longer(c(CMT2, FBX5, VIM2_4_DEL), names_to="Gene", values_to="Allele") %>%
  mutate(Gene = ifelse(Gene == "VIM2_4_DEL", "VIM2/VIM4", Gene)) %>%
  # CHANGED: Force Factor Order for Panel B
  mutate(Gene = factor(Gene, levels = target_order)) %>%
  group_by(Cluster, Group) %>%
  mutate(n = n_distinct(SeqID)) %>%
  group_by(Cluster, Group, n, Gene) %>%
  summarise(Freq = sum(Allele == "Derived", na.rm=TRUE) / n(), .groups="drop") %>%
  mutate(
    xlb = paste0(Group, "\n(n=", n, ")"),
    c_lbl = glue("<span style='color:{cols_c[Cluster]};'>{Cluster}</span>"),
    num = as.numeric(str_extract(Group, "\\d+"))
  ) %>%
  arrange(factor(Cluster, levels=ord_c), num) %>%
  mutate(
    Cluster = factor(Cluster, levels=ord_c),
    xlb = factor(xlb, levels=unique(xlb)),
    c_lbl = factor(c_lbl, levels=unique(c_lbl))
  )

p_b <- ggplot(df_b, aes(xlb, Freq, fill=Gene)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  scale_fill_manual(values=cols_g) +
  scale_y_continuous(labels=scales::percent, limits=c(0, 1.05), expand=c(0,0)) +
  facet_grid(~c_lbl, scales="free_x", space="free_x", switch="x") +
  labs(x=NULL, y="Derived allele frequency", fill=NULL, tag = "b") + # CHANGED to lowercase
  theme_bw(base_size=14) +
  theme(
    strip.text.x.bottom=element_markdown(face="bold", size=12, margin=margin(t=5, b=5)),
    strip.background=element_rect(fill="white", color="black"),
    strip.placement="outside",
    axis.text.x=element_text(angle=0, hjust=0.5, size=10, color="black"),
    axis.title.y=element_text(face="bold", margin=margin(r=10)),
    legend.position="top",
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    plot.tag = element_text(size = 20, face = "bold")
  )


# ==============================================================================
# Panel C: RAiSD histogram
# ==============================================================================
genes_c <- tibble(
  name = c("CMT2", "FBX5", "VIM2/VIM4"),
  chrom = c("Chr4", "Chr2", "Chr1"),
  pos = c(10420088, 18513626, 24586731)
)
cols_c_plot <- c("CMT2" = "navyblue", "FBX5" = "coral2", "VIM2/VIM4" = "palegreen4")

df_c <- map_dfr(files_raisd, read_tsv, col_names = c("p", "s", "e", "v", "sf", "ld", "stat"), col_types = cols(.default = "d"), skip = 1, .id = "chr")
scores_c <- genes_c %>% rowwise() %>%
  mutate(stat = df_c %>% filter(chr == chrom, s <= pos, e >= pos) %>% pull(stat) %>% max(na.rm = TRUE),
         # CHANGED: Force Factor Order for Legend in C
         name = factor(name, levels = target_order)) 
         
cut_c <- quantile(df_c$stat, 0.95, na.rm = TRUE)
mx_c <- max(ggplot_build(ggplot(df_c, aes(stat)) + geom_histogram(bins = 50))$data[[1]]$count)

p_c <- ggplot(df_c, aes(stat)) +
  geom_histogram(bins = 50, fill = "grey50", color = "black", alpha = 0.8) +
  geom_vline(xintercept = cut_c, linetype = 2, color = "black", linewidth = 1.2) +
  geom_vline(data = scores_c, aes(xintercept = stat, color = name), linetype = 1, linewidth = 1.2) +
  scale_color_manual(values = cols_c_plot, name = NULL) +
  labs(x = expression(bold("RAiSD " ~ mu ~ " statistic")), y = expression(bold(Count)), tag = "c") + # CHANGED to lowercase
  coord_cartesian(ylim = c(0, mx_c * 1.05), clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 1.5, color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dotted"),
    legend.position = "none",
    plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"),
    plot.tag = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(color="black", size=12, hjust=1, margin = margin(r = 2)),
    axis.title.y = element_text(face="bold", margin = margin(r = 5)),
    axis.ticks.length = unit(0.15, "cm")
  )


# ==============================================================================
# Panel D: Fst histogram
# ==============================================================================
genes_d <- tibble(
  name = c("CMT2", "FBX5", "VIM2/VIM4"),
  chr = c("Chr4", "Chr2", "Chr1"),
  pos = c(10420088, 18513626, 24586731)
)

df_d <- read_tsv(file_fst, show_col_types = FALSE) %>%
  mutate(FST = as.numeric(WEIGHTED_FST), chr_n = as.numeric(gsub("Chr", "", CHROM))) %>%
  filter(FST >= 0, !is.na(chr_n), !is.na(BIN_START))

cands_d <- genes_d %>% rowwise() %>%
  mutate(FST = df_d %>% filter(CHROM == chr, BIN_START <= pos, BIN_END >= pos) %>% pull(FST) %>% max(na.rm = TRUE),
         # CHANGED: Force Factor Order for Legend in D
         name = factor(name, levels = target_order))

cut95_d <- quantile(df_d$FST, 0.95, na.rm = TRUE)
ymax_d <- max(ggplot_build(ggplot(df_d, aes(FST)) + geom_histogram(bins = 50))$data[[1]]$count) * 1.05

p_d <- ggplot(df_d, aes(FST)) +
  geom_histogram(bins = 50, fill = "grey50", color = "black", alpha = 0.8) +
  geom_vline(xintercept = cut95_d, linetype = 2, linewidth = 1.2) +
  geom_vline(data = cands_d, aes(xintercept = FST, color = name), linetype = 1, linewidth = 1.2) +
  scale_color_manual(values = cols_c_plot, name = NULL) +
  labs(x = expression(bold(~F[ST])), y = expression(bold(Count)), tag = "d") + # CHANGED to lowercase
  coord_cartesian(ylim = c(0, ymax_d), clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dotted"),
    legend.position = "none",
    plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"),
    plot.tag = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(color="black", size=12, hjust=1, margin = margin(r = 2)),
    axis.title.y = element_text(face="bold", margin = margin(r = 5)),
    axis.ticks.length = unit(0.15, "cm")
  )


# ==============================================================================
# Panel E: Phenotypic Effects
# ==============================================================================
vars_e <- c("CMT2", "FBX5", "VIM2_4_DEL")
cols_e <- c("Ancestral" = "#606060", "CMT2" = "navyblue", "FBX5" = "coral2", "VIM2/VIM4" = "palegreen4")

df_e <- read.csv(file_pheno) %>%
  select(SeqID, FRI, all_of(vars_e), FloweringDays, delta13C_31_Avg) %>%
  pivot_longer(cols = all_of(vars_e), names_to = "Variant", values_to = "Genotype") %>%
  filter(!is.na(Genotype), !is.na(FRI)) %>%
  mutate(
    Genotype = factor(Genotype, levels = c("Ancestral", "Derived")),
    FRI = factor(FRI, levels = c("Ancestral", "Derived")),
    Var_Lab = ifelse(Variant == "VIM2_4_DEL", "VIM2/VIM4", Variant),
    # CHANGED: Force Factor Order for Panel E axis
    Var_Lab = factor(Var_Lab, levels = target_order),
    X_Lab = glue("{Var_Lab}<sub>{ifelse(Genotype=='Ancestral', 'Anc', 'Der')}</sub>"),
    Col_Grp = ifelse(Genotype == "Ancestral", "Ancestral", as.character(Var_Lab))
  )

x_ord_e <- df_e %>% arrange(Var_Lab, Genotype) %>% distinct(X_Lab) %>% pull(X_Lab)
df_e$X_Lab <- factor(df_e$X_Lab, levels = x_ord_e)

stats_ft <- df_e %>% filter(!is.na(FloweringDays)) %>% group_by(Variant, Var_Lab) %>%
  summarise(
    p = summary(lm(FloweringDays ~ Genotype + FRI))$coefficients["GenotypeDerived", 4],
    y = max(FloweringDays) + diff(range(FloweringDays)) * 0.15,
    group1 = X_Lab[Genotype == "Ancestral"][1],
    group2 = X_Lab[Genotype == "Derived"][1],
    .groups = "drop"
  ) %>% mutate(label = scales::pvalue(p, accuracy = 0.0001, add_p = TRUE))

stats_wue <- df_e %>% filter(!is.na(delta13C_31_Avg)) %>% group_by(Variant, Var_Lab) %>%
  summarise(
    p = summary(lm(delta13C_31_Avg ~ Genotype))$coefficients["GenotypeDerived", 4],
    y = max(delta13C_31_Avg) + diff(range(delta13C_31_Avg)) * 0.15,
    group1 = X_Lab[Genotype == "Ancestral"][1],
    group2 = X_Lab[Genotype == "Derived"][1],
    .groups = "drop"
  ) %>% mutate(label = scales::pvalue(p, accuracy = 0.0001, add_p = TRUE))

n_ft <- df_e %>% filter(!is.na(FloweringDays)) %>% count(X_Lab) %>% mutate(lab = paste0("n=", n))
n_wue <- df_e %>% filter(!is.na(delta13C_31_Avg)) %>% count(X_Lab) %>% mutate(lab = paste0("n=", n))

base_theme_e <- theme_classic(base_size = 9) + theme(
  legend.position = "none",
  axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 5)),
  axis.text.y = element_text(size = 10, color = "black"),
  axis.line = element_line(linewidth = 1.2),
  axis.ticks = element_line(linewidth = 1.2),
  axis.ticks.length = unit(0.2, "cm"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
)

p_e1 <- ggplot(df_e %>% filter(!is.na(FloweringDays)), aes(X_Lab, FloweringDays, fill = Col_Grp)) +
  geom_vline(xintercept = seq(2.5, length(x_ord_e) - 1.5, 2), color = "gray60", linewidth = 0.8) +
  geom_jitter(aes(color = Col_Grp), width = 0.2, shape = 16, size = 2, alpha = 0.6) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", color = "black", size = 0.8, fatten = 2.5) +
  stat_pvalue_manual(stats_ft, label = "label", y.position = "y", tip.length = 0.01, size = 3) +
  geom_text(data = n_ft, aes(X_Lab, -Inf, label = lab), vjust = -1.5, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = cols_e) + scale_color_manual(values = cols_e) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(y = "Days to flowering", x = NULL, tag = "e") + # CHANGED to lowercase
  base_theme_e +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        plot.margin = margin(t = 5.5, r = 5.5, b = -6, l = 5.5),
        plot.tag = element_text(size = 20, face = "bold"))

p_e2 <- ggplot(df_e %>% filter(!is.na(delta13C_31_Avg)), aes(X_Lab, delta13C_31_Avg, fill = Col_Grp)) +
  geom_vline(xintercept = seq(2.5, length(x_ord_e) - 1.5, 2), color = "gray60", linewidth = 0.8) +
  geom_jitter(aes(color = Col_Grp), width = 0.2, shape = 16, size = 2, alpha = 0.6) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", color = "black", size = 0.8, fatten = 2.5) +
  stat_pvalue_manual(stats_wue, label = "label", y.position = "y", tip.length = 0.01, size = 3) +
  geom_text(data = n_wue, aes(X_Lab, -Inf, label = lab), vjust = -1.5, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = cols_e) + scale_color_manual(values = cols_e) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(y = expression(bold(paste("WUE (", delta^13, "C)"))), x = NULL) +
  base_theme_e +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1, size = 9, face = "bold", color = "black"), 
        plot.margin = margin(t = -6, r = 5.5, b = 5.5, l = 5.5))

# models and stats
for (var in unique(df_e$Var_Lab)) {
  sub_df <- df_e %>% filter(Var_Lab == var)
  
  # 1. flowering time (controlled for FRIK232X)
  df_ft <- sub_df %>% filter(!is.na(FloweringDays))
  if(nrow(df_ft) > 0) {
    mod_ft <- summary(lm(FloweringDays ~ Genotype + FRI, data = df_ft))
    beta_ft <- mod_ft$coefficients["GenotypeDerived", "Estimate"]
    se_ft <- mod_ft$coefficients["GenotypeDerived", "Std. Error"]
    p_ft <- mod_ft$coefficients["GenotypeDerived", "Pr(>|t|)"]
    
    p_ft_format <- scales::pvalue(p_ft, accuracy = 0.001, add_p = TRUE)
  }
  
  # 2. Water Use Efficiency
  df_wue <- sub_df %>% filter(!is.na(delta13C_31_Avg))
  if(nrow(df_wue) > 0) {
    mod_wue <- summary(lm(delta13C_31_Avg ~ Genotype, data = df_wue))
    beta_wue <- mod_wue$coefficients["GenotypeDerived", "Estimate"]
    se_wue <- mod_wue$coefficients["GenotypeDerived", "Std. Error"]
    p_wue <- mod_wue$coefficients["GenotypeDerived", "Pr(>|t|)"]
    
    p_wue_format <- scales::pvalue(p_wue, accuracy = 0.001, add_p = TRUE)
  }
  
  # print stats
  cat(paste0("\n--- ", var, " ---\n"))
  cat(sprintf("Flowering Time: (beta = %.2f, SE = %.2f, %s)\n", beta_ft, se_ft, p_ft_format))
  cat(sprintf("WUE (delta13C): (beta = %.2f, SE = %.2f, %s)\n", beta_wue, se_wue, p_wue_format))
}

# ==============================================================================
# Final figure with all panels
# ==============================================================================
col_right <- plot_grid(p_e1, p_e2, ncol = 1, align = "v", rel_heights = c(1, 1.2))
col_left <- plot_grid(p_c, p_d, ncol = 1, align = "v")
bottom_row <- plot_grid(col_left, col_right, ncol = 2, rel_widths = c(1, 1))
final_plot <- plot_grid(p_a, p_b, bottom_row, ncol = 1, rel_heights = c(1.5, 1, 2.5))

# Save
ggsave("Figure7.pdf", final_plot, width = 14, height = 16)
```

# Authors

- **Johan Zicola** - [johanzi](https://github.com/johanzi)
- **Emmanuel Tergemina** -
  [EmmanuelTergemina](https://github.com/EmmanuelTergemina)
- **Ahmed F. Elfarargi** -
  [AhmedElfarargi](https://github.com/AhmedElfarargi)

# License

This project is licensed under the MIT License - see the
[LICENSE](LICENSE) file for details
