#!/bin/bash

# AUTHOR: Johan Zicola
# DATE: 2018-01-26

# DESCRIPTION: This script performs GWAS using as input a vcf file and a phenotype file
# The phenotype file should contain as many values as individuals in the vcf file.
# The order of the accessions in the vcf file should be the same than the order of the 
# phenotype in the phenotype file

# USAGE:
# run_gwas_gemma.sh <phenotype_file.tsv> <vcf_file.vcf> [covariate_file.txt]

# Output is a .assoc.clean.txt file which can be loaded in R to look at significant SNPs.
# All generated files are located where <phenotype_file.tsv> is and the subdirectory "output"

# Path to python script assoc2qqman.py (assumes the script is in the same
# directory as this script)

# Note that gemma output files will be generated in the directory
# in which the user launches the script (not forcibly the directory containing the script or the input files)
# We call this directory $current_path
current_path=$(pwd)

# Get directory of running bash script to identify where the python script assoc2qqman.py is 
path_script=$(dirname $0)
path_script=$(cd $path_script && pwd)
# Assign path of the python script assoc2qqman.py to $assoc2qqman
assoc2qqman="${path_script}/assoc2qqman.py"

# Softwares required:
#	python2.7
#	gemma 0.94 (tested)
#	vcftools 0.1.14 (tested)
#	bcftools 1.2 (tested)
#	PLINK v1.90b3.40 (tested)
#
####################################################################################

# Help command if no arguments or -h argument is given 
if [ "$1" == "-h" -o "$#" -eq 0 ] ; then
	echo -e "\n" 
	echo "Usage: `basename $0` <phenotype_file.tsv> <vcf_file.vcf> [covariate_file.txt] [-h]"
	echo -e "\n" 
	echo "Description: Provide as main argument a phenotype file" 
	echo "generated in R (should have a .tsv extension) and a vcf files containing "
	echo "the accessions described in the phenotype file (should have a .vcf extension)."
	echo "Optional argument "covariate_file.txt" can be provided in third position"
	echo "if covariate analysis is to be performed."
	echo "`basename $0` generates several files derived from vcftools, p-link, and gemma"
	echo "The file with ".assoc.clean.txt" as suffix contains the GWAS results and can be "
	echo "uploaded in R for visualization (qqman library)"
	echo -e "\n"
exit 0
fi

# Test if softwares are available
test_program(){
	command -v $1 >/dev/null 2>&1 || { echo >&2 " $1 required  but it's not installed. Aborting"; exit 1; }
}

test_program "python"
test_program "gemma"
test_program "vcftools"
test_program "bcftools"
test_program "plink" 

# Time script
start=$(date +%s)

# Test if 2 arguments were provided
if [ "$#" -lt 2	]; then
	echo "Argument(s) missing"
	exit 1
fi

# Phenotype file is generated by R and has a tsv extension (can have several phenotypes
# tab-separated in different columns. It should be in unix format and should not 
# contai empty line
phenotype_file=$1
if [ ! -e $phenotype_file ]; then
	echo "File $phenotype_file does not exist"
    exit 1
elif [[ $phenotype_file != *.tsv ]]; then
	echo "Provide a phenotype file with .tsv extension"
	exit 1
fi 

# Define location of vcf file to analyze
vcf_file=$2
if [ ! -e $vcf_file ]; then
	echo "File $vcf_file does not exist"
     	exit 0
elif [[ $vcf_file != *.vcf ]] && [[ $vcf_file != *.vcf.gz ]]; then
	echo "Provide a vcf file with .vcf or vcf.gz extension"
	exit 0
fi

# Check if phenotype file contains as many lines as there are samples in the VCF file.
# If not, this will yield in segmentation fault when making the kinship matrix

nb_samples=$(bcftools query -l $vcf_file | wc -l)
nb_phenotypes=$(cat $phenotype_file | wc -l)
if [[ $nb_samples != $nb_phenotypes ]]; then
	echo "Number of samples and number of phenotypes are not equal! verify files"
	exit 0
fi

# Check if covariate file is provided
if [ "$#" -eq 3 ];	then
	covariate_file=$3
	# Check if file exists
	if [ -e $covariate_file ]; then		
		# Check if number of lines equal number of samples
		line_covariate_file=$(cat $covariate_file | wc -l)
		if [[ $line_covariate_file != $nb_phenotypes ]]; then
			echo "Number of lines in covariate file and number of phenotypes are not equal! verify files"
			exit 0
		else
			# Get prefix of the covariate file to add as suffix in final output name
			covariate_prefix="${covariate_file%.*}"
		fi
	else
			echo "Covariate file $covariate_file does not exist"
			exit 0
	fi
fi

	
##################################################################################

# Directory containing the phenotype file and which will contain final GWAS results
# in the output/ directory
# Get prefix of phenotype file
dir_file=$(pwd $phenotype_file)

# Get prefix from phenotype name (assume the phenotype file has a .tsv extension)
# Prefix for GWAS results. Based on phenotype file name
prefix_gwas=$(basename -s .tsv $phenotype_file)

# Prefix to generate plink files from input VCF file
# Per default will be name of the VCF file before the first point
path_prefix_vcf=$(echo $vcf_file | cut -d'.' -f1)

# In case the file of the VCF file contains a path, get suffix of the vcf file
# for output name for gemma
prefix_vcf=$(basename $path_prefix_vcf)

echo -e "\n###################### INPUT VARIABLES ###################################\n"

echo "phenotype file: $phenotype_file"
echo "vcf file: $vcf_file"
echo "dir_file (directory where output/ will be): $dir_file"
echo "current_path: $current_path"
echo "prefix plink output files: $prefix_vcf"
echo "prefix gemma output files: $prefix_gwas"

echo -e "\n###################### CONVERT VCF TO PLINK FORMAT #######################\n"


# VCF into bed file => make .ped and .map files
# These files must be made only once, then only the fam file should be modified for the 
# tested phenotype
# Check if plink files already exist
# Also check if input vcf file is compressed or not

echo -e "Generate ped and map files:"
if [ -e ${path_prefix_vcf}.ped ] && [ -e ${path_prefix_vcf}.map ]; then
	echo -e "${path_prefix_vcf}.ped and ${path_prefix_vcf}.map already exists. Go to next step."
else
	if [[ $vcf_file == *.vcf ]]; then
		printf "vcftools --vcf $vcf_file --plink --out ${path_prefix_vcf}"
		vcftools --vcf $vcf_file --plink --out ${path_prefix_vcf}
	elif [[ $vcf_file == *.vcf.gz ]]; then
		printf "vcftools --gzvcf $vcf_file --plink --out ${path_prefix_vcf}"
		vcftools --gzvcf $vcf_file --plink --out ${path_prefix_vcf} 
	fi
fi

echo -e "\nGenerate bed, bim, and fam files:"
# Make bed files: 3 files are created => .bed, .bim, .fam
if [ -e ${path_prefix_vcf}.bed ] && [ -e ${path_prefix_vcf}.bim ] && [ -e ${path_prefix_vcf}.fam ]; then
	echo "File ${path_prefix_vcf}.bed, ${path_prefix_vcf}.bim, ${path_prefix_vcf}.fam already exist. Go to next step."
else
	printf "plink --file ${path_prefix_vcf} --make-bed --out ${path_prefix_vcf} --noweb" 
	plink --file ${path_prefix_vcf} --make-bed --out ${path_prefix_vcf} --noweb
fi


echo -e "\nKeep only 5 first column from .fam file (remove 6th column of -9)"

echo "cut -d' ' -f1,2,3,4,5 ${path_prefix_vcf}.fam > ${path_prefix_vcf}_modified.fam"
cut -d' ' -f1,2,3,4,5 ${path_prefix_vcf}.fam > ${path_prefix_vcf}_modified.fam


echo -e "\nPaste phenotype data to fam file and reformat it:"
# Paste to fam file
echo "paste -d ' ' ${path_prefix_vcf}_modified.fam $phenotype_file > ${path_prefix_vcf}.fam"
paste -d ' ' ${path_prefix_vcf}_modified.fam $phenotype_file > ${path_prefix_vcf}.fam


echo -e "\n############################# RUN GEMMA ##################################\n"
# Run Gemma
# Per default, gemma put the results in an "output" directory located in working directory ($current_path)
# There is apparently no way to change this (adding fullpath in -o  /srv/biodata/dep_coupland/grp_hancock/johan/GWAS/rDNA_copy_number_MOR)
# does not work: 
# error writing file: ./output//srv/biodata/dep_coupland/grp_hancock/johan/GWAS/rDNA_copy_number_MOR.sXX.txt
# Instead, just transfer the generated files into the ${dir_file}/output directory at the end

# Estimate relatedness matrix from genotypes (n x n individuals)
# Generate relatedness matrix based on centered genotypes (sXX)
# centered matrix preferred in general, accounts better for population structure
# If standardized genotype matrix is needed, change to -gk 2 (sXX)
# standardized matrix preferred if SNPs with lower MAF have larger effects 
# Note that this file remains the same whatever phenotype is studied as it considers 
# only the VCF file.
# Add flag (-maf 0.05) to set it at 5% for instance).

echo -e "Generate relatedness matrix:"
if [ -e ${current_path}/output/${prefix_vcf}.sXX.txt ]; then
	echo -e "${current_path}/output/${prefix_vcf}.sXX.txt file already exists. Go to next step."
else
	echo -e "\ngemma -bfile ${path_prefix_vcf} -gk 2 ${prefix_vcf}"
	gemma -bfile ${path_prefix_vcf} -gk 2 -o $prefix_vcf
fi

## If needed, the relatedness matrix can be transformed into eigen values and eigen vectors
## Generates 3 files: log, eigen values (1 column of na elements) and eigen vectors  (na x na matrix)
## Use of eigen transformation allows quicker analysis (if samples > 10 000)
# gemma -bfile ${prefix} -k ${current_path}/output/${prefix}.sXX.txt -eigen -o ${prefix}

## Association Tests with Univariate Linear Mixed Models
# Use lmm 2 to performs likelihood ratio test
# Note that minor allele frequency (MAF( threshold is set at 1% per default. 
# Add flag (-maf 0.05) to set it at 5% for instance).
# prefix.log.txt contains PVE estimate and its standard error in the null linear mixed model.
# assoc.txt file contains the results
# If covariate file is provided as third argument, perform covariate analysis

echo -e "\nPerform the association test:"

# Rename variable $prefix_gwas if covariate file was provided so
# the output can be differentiated from same analysis without covariate file
if [ $covariate_file ]; then
	prefix_gwas="${prefix_gwas}_${covariate_prefix}"
fi

if [ -e ${current_path}/output/${prefix_gwas}.assoc.txt ]; then
	echo "${current_path}/output/${prefix_gwas}.assoc.txt already exists. Go to next step"
else
	if [ $covariate_file ]; then	
		echo -e "\ngemma -bfile ${path_prefix_vcf} -k ${current_path}/output/${prefix_vcf}.sXX.txt -lmm 2 -maf 0.05 -o ${prefix_gwas} -c $covariate_file"
		gemma -bfile ${path_prefix_vcf} -k ${current_path}/output/${prefix_vcf}.sXX.txt -lmm 2 -maf 0.05 -o ${prefix_gwas} -c $covariate_file
	else
		echo -e "\ngemma -bfile ${path_prefix_vcf} -k ${current_path}/output/${prefix_vcf}.sXX.txt -lmm 2 -maf 0.05 -o ${prefix_gwas}"
		gemma -bfile ${path_prefix_vcf} -k ${current_path}/output/${prefix_vcf}.sXX.txt -lmm 2 -maf 0.05 -o ${prefix_gwas}
	fi
fi

# # Association Tests with Multivariate Linear Mixed Models
# # several phenotypes can be given (for instance columns 1,2,3 of the phenotype file). Less than 10
# # phenotypes are recommended
# gemma -bfile ${path_prefix_vcf} -k ${current_path}/output/${prefix}.sXX.txt -lmm 2 -n 1 2 3 -o ${prefix}

## Bayesian Sparse Linear Mixed Model
## Use a standard linear BSLMM (-bslmm 1)
## Does not require a relatedness matrix (calculates it internally)
## Generates 5 output files: log, hyp.txt (estimated hyper-parameters), param.txt (posterior mean 
## estimates for the effect size parameters), prefix.bv contains 
# gemma -bfile ${path_prefix_vcf} -bslmm 1 -o ${prefix}

# To analyze the output of bslmm in R, one needs to plot the gamma value multiplied by beta value from the # param.txt (see more details on https://visoca.github.io/popgenomworkshop-gwas_gemma/)
# see R script gemma_param.R


echo -e "\n########################## POLISHING GEMMA OUTPUT ############################\n"

# Polish file for R
echo -e "Reformat assoc.txt file to be compatible with manhattan library in R:"

echo "python $assoc2qqman ${current_path}/output/${prefix_gwas}.assoc.txt > ${current_path}/output/${prefix_gwas}.assoc.clean.txt"
python $assoc2qqman ${current_path}/output/${prefix_gwas}.assoc.txt > ${current_path}/output/${prefix_gwas}.assoc.clean.txt

# Move the output data into dir_file

# Check if output directory is present in dir_file (case when script is $dir_file == $path_script
# If not, create it and move output/ in it

# Transfer output files if necessary (in case $current_path != $dir_path
if [[ $current_path != $dir_file ]]; then
	echo "Transfer output files in ${dir_file}/output"
	echo "mkdir ${dir_file}/output"
	mkdir ${dir_file}/output
	echo "mv ${current_path}/output/* ${dir_file}/output/"
	mv ${current_path}/output/* ${dir_file}/output/
	echo "rm -r ${current_path}/output"
	rm -r ${current_path}/output
fi

end=$(date +%s)
runtime=$((end-start))

echo -e "\n############################ GENERATE LOG FILE  #################################\n"

# Create a log output file and integrate gemma log file

# Check if log from gemma was successfully generated
if [ -e ${current_path}/output/${prefix_gwas}.log.txt ]; then
	echo -e "\nLog generated in ${current_path}/output/${prefix_gwas}.log"
	echo -e "LOG FILE\n" > ${current_path}/output/temp.txt
	echo "Command: run_gwas_gemma.sh $phenotype_file $vcf_file" >> ${current_path}/output/temp.txt
	echo "Phenotype file: $phenotype_file" >> ${current_path}/output/temp.txt
	echo "VCF file: $vcf_file" >> ${current_path}/output/temp.txt
	if [ $covariate_file ]; then
		echo "Covariate file: $covariate_file" >> ${current_path}/output/temp.txt
	fi
	echo "Output file: ${dir_file}/output/${prefix_gwas}.assoc.clean.txt" >> ${current_path}/output/temp.txt
	echo -e "Run finished on $(date)" >> ${current_path}/output/temp.txt
	echo -e "Total time of the run: $runtime seconds \n" >> ${current_path}/output/temp.txt
	echo "Log output from GEMMA:" >>  ${current_path}/output/temp.txt
	cat ${current_path}/output/temp.txt ${current_path}/output/${prefix_gwas}.log.txt > ${current_path}/output/${prefix_gwas}.log
	rm ${current_path}/output/temp.txt ${current_path}/output/${prefix_gwas}.log.txt
else
	echo -e "Gemma did not generate a log, a problem occured in the pipeline!"
fi


echo -e "\n############################ END OF THE PIPELINE  #################################\n"

