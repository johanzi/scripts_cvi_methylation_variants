#!/bin/bash
#===============================================================================#
#		  FILE: combine_kmers_methylome_SA.sh
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

KMC='/srv/netscratch/irg/grp_hancock/Manu/software/kmersGWASv0.2/external_programs/kmc_v3'
BIN='/srv/netscratch/irg/grp_hancock/Manu/software/kmersGWASv0.2/bin/'
DIR='/srv/netscratch/irg/grp_hancock/johan/kmerGWAS/'



# ${BIN}list_kmers_found_in_multiple_samples -l ${DIR}kmers_list_paths.txt -k 31 --mac 5 -p 0.2 -o ${DIR}kmers_to_use
# ${BIN}build_kmers_table -l ${DIR}kmers_list_paths.txt -k 31 -a kmers_to_use -o ${DIR}kmers_table
# ${BIN}emma_kinship_kmers -t ${DIR}kmers_table -k 31 --maf 0.05 > ${DIR}kmers_table_MAF10.kinship



${BIN}list_kmers_found_in_multiple_samples -l ${DIR}kmers_list_paths.txt -k 31 --mac 6 -p 0.2 -o ${DIR}kmers_to_use
${BIN}build_kmers_table -l ${DIR}kmers_list_paths.txt -k 31 -a kmers_to_use -o ${DIR}kmers_table
${BIN}emma_kinship_kmers -t ${DIR}kmers_table -k 31 --maf 0.1 > ${DIR}kmers_table.kinship
