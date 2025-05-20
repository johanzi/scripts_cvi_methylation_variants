# First argument => name of rda object
# Second argument => name of context
# Third argument => region
# Fourth argument => accessions_order file (should not have header)
# Fifth argument => 'sample' name (as given in rda file) and 'seq_ID' (as in VCF file) for each sample (header needed)

# The script can be run from another directory as long as dependent scripts are located in same directory as caller script (i.e. this script)

ABSPATH=$(readlink -f $0)
ABSDIR=$(dirname $ABSPATH)

generate_phenotype_gemma(){
	if [[ "$#" != 5 ]]; then
		echo "Arguments are missing"
		exit 0
	fi
	
	rda_object=$1
	name_df=$(basename ${1%.rda})
	context=$2
	region=$3
	gemma_file=${context}_${region}.tsv
	accessions_order=${4}
	accessions_info=${5}

	# Make readable dataframe for 1 methylation context
	echo -e "Rscript ${ABSDIR}/subset_df_methylation.R $context $rda_object ${name_df}.txt\n"
	Rscript ${ABSDIR}/subset_df_methylation.R $context $rda_object ${name_df}.txt
	
	# Generate a file with only methylation data, in the order specified in the VCF file used for GWAS
	echo -e "Rscript ${ABSDIR}/create_gemma_phenotype.R $accessions_order $accessions_info ${name_df}.txt $gemma_file \n"
	Rscript ${ABSDIR}/create_gemma_phenotype.R $accessions_order $accessions_info ${name_df}.txt $gemma_file
	
	# Delete the temporary ${name_df}.txt file
	rm -v ${name_df}.txt
}

# Example:
# generate_phenotype_gemma df_mean_intergenic_regions.rda CHH intergenic_regions order_accessions.txt
generate_phenotype_gemma $1 $2 $3 $4 $5

