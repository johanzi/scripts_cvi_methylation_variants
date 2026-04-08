#!/bin/bash
#===============================================================================#
#		  FILE: run_GWASkmers_methylome_SA.sh
#		 USAGE: ---
#  DESCRIPTION: ---
#	   OPTIONS: ---
# REQUIREMENTS: ---
#		  BUGS: ---
#		 NOTES: ---
#	    AUTHOR: Emmanuel Tergemina
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#	   VERSION: 1.0
#	   CREATED: 20220824
#	  REVISION: ---
#===============================================================================

kmerGWAS='/srv/netscratch/irg/grp_hancock/Manu/software/kmersGWASv0.2'

files=(/srv/netscratch/irg/grp_hancock/johan/kmerGWAS/phenotype/*.tsv)


# Multiple jobs
# index=$LSB_JOBINDEX
# var=$(expr $index - 1)
# echo `basename ${files[$var]} .tsv`
# # echo ${files[$var]}
# python2.7 $kmerGWAS/kmers_gwas.py --pheno ${files[$var]} --kmers_table kmers_table -l 31 -p 8 --outdir output_dir/`basename ${files[$var]} .tsv`


# Single jobs
# $1 = first argument
# python2.7 $kmerGWAS/kmers_gwas.py --pheno phenotype/$1.tsv --kmers_table kmers_table -l 31 -p 8 --outdir output_dir/$1

#
python2.7 $kmerGWAS/kmers_gwas.py --pheno phenotype/$1.tsv --kmers_table kmers_table -l 31 -p 8 --maf 0.1 --mac 8 --outdir output_dir/$1