# Functions for bismark/methylKit analysis pipeline

# Written: 2018-07-10
# Author: Johan Zicola


# Libraries ---------------------------------------------------------------

## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("methylKit")
# biocLite("GenomicRanges")
# biocLite("genomation")

# New method with R3.6
#BiocManager::install("methylKit")
#BiocManager::install("GenomicRanges")
#BiocManager::install("genomation")

library(methylKit)
library(genomation)
library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(data.table)
library(knitr)
library(qqman)
library(stringr)
library(dplyr)
library(onewaytests)
#library(gridExtra)

# Check help on https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#12_high-throughput_bisulfite_sequencing


# General functions -------------------------------------------------------

# Get weighted methylation levels mean value by parsing the output of getMethylationStats function
# The mean value represents the weighted methylation level (Schultz 2012)
# It represents the total of methylated cytosines (C) divided by the total of calls (C+T)
get_wml <- function(object){
  summary <- capture.output(getMethylationStats(object, both.strands=FALSE, plot = FALSE))
  # Parse and get the mean element of the function
  mean <- as.list(strsplit(summary[4], "\\s+")[[1]][[5]])
  mean <- as.double(mean)
}


# Make a dataframe with name of the accessions, methylation context and weighted methylation level
# Requires a list of methylRawList (can be DB or not)
# Requires function getAverageMethylation
get_df_wml <- function(list_methylRawLists, path_DB, df_name){
  
  # Get path name to object of dataframe
  df_rda <- paste(path_DB, paste(df_name, ".rda", sep=""), sep="/")
  
  # Check if R object exists, if not, create it
  if (file.exists(df_rda)){
    warning(sprintf('The file "%s" already exist. Erase it first if needed to be regenerated.', df_rda))
  } else {

    # How many contexts are studied
    nb_contexts <- length(list_methylRawLists)
    
    # Make a empty vector of contexts
    list_contexts <- c()
    
    # I take the context from the first sample (Data[[1]]) but it does not matter which
    for (i in 1:nb_contexts){
      list_contexts[i] <- list_methylRawLists[[i]]@.Data[[1]]@context
    }
    
    # Get nb of samples for the first methylation context (should be the same for each context)
    nb_samples <- length(list_methylRawLists[[1]])
    
    # Make a vectors of samples
    list_samples <- c()
    
    # Retrieve the sample names
    for (i in 1:nb_samples){
      list_samples[i] <- list_methylRawLists[[1]]@.Data[[i]]@sample.id
    }
    
    # Get methylation level
    list_means <- c()
    index=1
    
    for (i in 1:nb_contexts){
      for (j in 1:nb_samples){
        list_means[index] <- get_wml(list_methylRawLists[[i]][[j]])
        index = index + 1
      }
    }
    
    # Make a dataframe of averages (all elements should be vectors)
    # Repeat vectors list_samples and list_contexts to have as many rows as for methylation
    # means.
    samples <- rep(list_samples, nb_contexts)
    contexts <- rep(list_contexts, each=nb_samples)
    df_mean <- data.frame(samples, contexts, list_means)
    
    
    # Rename headers of dataframe
    colnames(df_mean) <- c("sample","context","percent_methylation")
    
    # Save the dataframe (df_name should still be a string in this step)
    df_rda <- paste(path_DB, paste(df_name, ".rda", sep=""), sep="/")
    
    # Assign dataframe to the string variable df_name (becomes a df)
    assign(df_name, df_mean)
  
    # Use list= to avoid saving the string value instead of the df R object
    save(list=df_name, file=df_rda)
    #load(file=df_rda, envir= .GlobalEnv)
  }
}


# Load dataframe with weighted methylation level or created it if not existing
load_df_wml <- function(path_DB, df_name){
  # Get path name to object of dataframe
  df_rda <- paste(path_DB, paste(df_name, ".rda", sep=""), sep="/")

    # Check if R object exists, if not, create it
  if (file.exists(df_rda)){
      # Load dataframe as R object
      load(file=df_rda, envir= .GlobalEnv)
  } else {
      stop(sprintf('The file "%s" does not exist. Generate it first with get_df_wml', df_rda))
  }
}


reformat_df_methylation <- function(df){
  
  # Get plot without CX
  df <- df[(df$context %in% c("CpG","CHG","CHH")),]
  
  # Get accession name w/o library code
  sample <- lapply(as.vector(df$sample), function(x) strsplit(x, "...._.+_")[[1]][[2]])
  sample <- unlist(sample)
  df$sample <- sample
  df$sample <- as.factor(df$sample)
  return(df)
}


# Plot weighted methylation levels
ggplot_all <- function(df, title="Weighted Methylation Level"){
  ggplot(data = df, aes(x=sample, y=percent_methylation)) + 
    geom_bar(aes(x=sample, y=percent_methylation, fill=context), position="dodge", stat="identity") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


# Plot getCoverageStats in batch
getCoverageStats_batch <- function(list_methylRawLists){
  
  # How many contexts are studied
  nb_contexts <- length(list_methylRawLists)
  
  # Get nb of samples for the first methylation context (should be the same for each context)
  nb_samples <- length(list_methylRawLists[[1]])
  
  for (i in 1:nb_contexts){
    for (j in 1:nb_samples){
      getCoverageStats(list_methylRawLists[[i]][[j]], plot=TRUE, both.strands=FALSE)
    }
  }
}


# Plot getCoverageStats in batch and displays values instead of an histogram
# Choose context to output
getCoverageStats_values_batch <- function(list_methylRawLists, context){
  
  # Context studied (1,2,3)
  i <- context
  
  # Get nb of samples for the first methylation context (should be the same for each context)
  nb_samples <- length(list_methylRawLists[[1]])
  
  cat("Sample\tContext\tMin\t1st_Qu\tMedian\tMean\t3rd_Qu\tMax") 
  cat("\n")
  
  for (j in 1:nb_samples){
    # Name sample
    name <- list_methylRawLists[[i]]@.Data[[j]]@sample.id
    # Name context
    context <- list_methylRawLists[[i]]@.Data[[j]]@context
    # Get summary stat for a subelement
    summary <- capture.output(getCoverageStats(list_methylRawLists[[i]][[j]]))
    # Extract mean
    #mean <- unlist(as.list(strsplit(summary[[4]],"\\s+")[[1]][[4]]))
    values <- as.list(strsplit(summary[[4]],"\\s+")[[1]][c(2,3,4,5,6,7)])
    cat(name, context, unlist(values), sep="\t")
    #cat(name, context, mean, sep="\t")
    cat("\n")
  }
}



# Display plots with all samples but highlighted CMT2 allele (defined as variable group)
# Provide possibility to order the data by ascending level of methylation
ggplot_all_allele <- function(df, context, group=NULL, order_data=FALSE){
  # Get subset of data with wanted methylation context
  df_context <- df[(df$context == context),]
  
  # Order the methylation data in ascending order if optional argument 'order_data' is TRUE
  if((isTRUE(order_data))){
    df_context$sample <- factor(df_context$sample, levels = df_context$sample[order(df_context$percent_methylation)])
  }
  
  df_context <- as.data.frame(df_context)
  
  # Check if group is given
  if((is.null(group))){
    ggplot(data=df_context) +
      geom_bar(aes(x=sample, y=percent_methylation), stat="identity") +
      ggtitle(paste(context," methylation",sep="")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else {
  ggplot(data=df_context) +
    geom_bar(aes_string(x="sample", y="percent_methylation", fill=group),
             position="dodge", stat="identity") +
    ggtitle(paste(context," methylation",sep="")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
}



# Display boxplot of the methylation in chosen context for the 2 alleles of CMT2
# Display number of accession for each allele

ggplot_per_context <- function(df, context, group){
  
  give.n <- function(x){
    return(c(y = mean(x), label = length(x)))
  }
  
  df_context <- df[(df$context == context),]

  ggplot(data=df_context, aes_string(x=group, y="percent_methylation", group=group))+
    ggtitle(paste(context, " methylation", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_boxplot() +
    stat_summary(fun.data = give.n, geom = "text")
  
}


# Plot getMethylationStats in batch
getMethylationStats_batch <- function(list_methylRawLists){
  
  # How many contexts are studied
  nb_contexts <- length(list_methylRawLists)
  
  # Get nb of samples for the first methylation context (should be the same for each context)
  nb_samples <- length(list_methylRawLists[[1]])
  
  # Loop over the list_methylRawLists
  for (i in 1:nb_contexts){
    for (j in 1:nb_samples){
      getMethylationStats(list_methylRawLists[[i]][[j]],plot=TRUE,both.strands=FALSE)
    }
  }
}


# Get list of files according to the type precised (either raw or filtered)
# Upgrade by allowing the loading of subset regions using type as pattern
get_list_files <- function(path_DB, list_samples, type){
  if (type == "raw"){
    pattern=".txt.bgz"
  } else {
    pattern=paste("_", type, ".txt.bgz", sep="")
  }
  list_files <- lapply(list_samples, function(x) x = paste(path_DB,"/", x, pattern, sep=""))
}


# Make tiles for DMR analysis
tileMethylCounts_batch <- function(list_methylBases, win_size=300, step_size=100, suffix=NULL, list_dbdir=NULL, cov.bases=0){
  
  # Check if list_dbdir is provided
  if(missing(list_dbdir)){
    stop("Argument list_dbdir not provided")
  }
  
  # Suffix indicates windows and step sizes, additional suffix if needed is optional
  suffix <- paste("tiled", suffix, win_size, step_size, sep="_")
  
  # Loop over the list_methylBases object
  for (i in 1:3){
    methylBase_object <- DB_to_RAM_conversion(list_methylBases[[i]], type="methylBase")
    tileMethylCounts(methylBase_object, 
                     win.size=win_size, step.size=step_size, 
                     save.db=TRUE, dbdir=list_dbdir[[i]], suffix=suffix, cov.bases = cov.bases)
    rm(methylBase_object)
  }
}

# For instance, a file called methylBase_Col-0_vs_Cvi-0.txt.bgz would be named
# methylBase_Col-0_vs_Cvi-0_tiled_300_100.txt.bgz per default.


# Identify DMRs
get_DMRs <- function(list_methylBases){
  
  list_diff <- lapply(list_methylBases, calculateDiffMeth, overdispersion="MN",test="Chisq",mc.cores=1, save.db = TRUE)
  
  list_hypo <- lapply(list_diff, getMethylDiff, difference=25, qvalue=0.01, type="hypo", save.db=TRUE)
  
  list_hyper <- lapply(list_diff, getMethylDiff, difference=25, qvalue=0.01, type="hyper", save.db=TRUE)
  
  list_DMRs <- list(list_diff, list_hypo, list_hyper)
  
  names(list_DMRs) <- c("all_DMRs", "hypo_DMRs", "hyper_DMRs")
}


# Corrected version of .setMethylDBNames. Will be fixed in the new methylKit version
setMethylDBNames_new <- function (df, methylDBclass = c("methylRawDB", "methylBaseDB", 
                                "methylDiffDB")) 
{
  if (nrow(df) == 0) 
    return(df)
  if (missing(methylDBclass)) {
    if (length(df) == 7 & unique(sapply(df, class)[5:7]) == 
        "integer") {
      setnames(x = df, old = names(df), new = c("chr", 
                                                "start", "end", "strand", "coverage", "numCs", 
                                                "numTs"))
    }
    else if (length(df) == 7 & unique(sapply(df, class)[5:7]) == 
             "numeric") {
      setnames(x = df, old = names(df), new = c("chr", 
                                                "start", "end", "strand", "pvalue", "qvalue", 
                                                "meth.diff"))
    }
    else if (length(df) > 7) {
      setnames(x = df, old = names(df)[1:4], new = c("chr", 
                                                     "start", "end", "strand"))
      numsamples = (length(df) - 4)/3
      coverage.ind = seq(5, by = 3, length.out = numsamples)
      numCs.ind = coverage.ind + 1
      numTs.ind = coverage.ind + 2
      setnames(df, names(df)[coverage.ind], paste(c("coverage"), 
                                                  1:numsamples, sep = ""))
      setnames(df, names(df)[numCs.ind], paste(c("numCs"), 
                                               1:numsamples, sep = ""))
      setnames(df, names(df)[numTs.ind], paste(c("numTs"), 
                                               1:numsamples, sep = ""))
    }
  }
  else {
    if (methylDBclass == "methylRawDB") {
      setnames(x = df, old = names(df), new = c("chr", 
                                                "start", "end", "strand", "coverage", "numCs", 
                                                "numTs"))
    }
    else if (methylDBclass == "methylBaseDB") {
      setnames(x = df, old = names(df)[1:4], new = c("chr", 
                                                     "start", "end", "strand"))
      numsamples = (length(df) - 4)/3
      coverage.ind = seq(5, by = 3, length.out = numsamples)
      numCs.ind = coverage.ind + 1
      numTs.ind = coverage.ind + 2
      setnames(df, names(df)[coverage.ind], paste(c("coverage"), 
                                                  1:numsamples, sep = ""))
      setnames(df, names(df)[numCs.ind], paste(c("numCs"), 
                                               1:numsamples, sep = ""))
      setnames(df, names(df)[numTs.ind], paste(c("numTs"), 
                                               1:numsamples, sep = ""))
    }
    else if (methylDBclass == "methylDiffDB") {
      setnames(x = df, old = names(df), new = c("chr", 
                                                "start", "end", "strand", "pvalue", "qvalue", 
                                                "meth.diff"))
    }
  }
}



# Convert a methylKitDB object to methylKit object (put data from DB files in RAM)
DB_to_RAM_conversion <- function(x, type=c("methylRaw","methylBase","methylDiff")) {
  
  df <- methylKit:::fread.gzipped(x@dbpath,
                                  stringsAsFactors = FALSE,
                                  data.table = FALSE,
                                  skipDecompress = FALSE)
  #methylKit:::.setMethylDBNames(df,paste(type,"DB",sep=""))
  setMethylDBNames_new(df, methylDBclass=paste(type,"DB",sep=""))
  
  if (type == "methylBase"){
    new(type,
        df,
        sample.ids = x@sample.ids,
        assembly = x@assembly,
        context = x@context,
        treatment = x@treatment,
        coverage.index = x@coverage.index,
        numCs.index = x@numCs.index,
        numTs.index = x@numTs.index,
        destranded = x@destranded,
        resolution = x@resolution
      )
  }  else if (type == "methylRaw") {
    new(type,
        df,
        sample.id = x@sample.id,
        assembly = x@assembly,
        context = x@context,
        resolution = x@resolution
      )
    } else if (type == "methylDiff") {
      new(type,
          df,
          sample.ids = x@sample.ids,
          assembly = x@assembly,
          context = x@context,
          treatment = x@treatment,
          destranded = x@destranded,
          resolution = x@resolution
        )
    }
}



# GWAS code -------------------------------

GWAS_run <- function(output_gemma, threshold_pvalue="0", highlighted_SNP=""){
  
  # Highlighted_SNP allows to display in green the SNP of interested on the Manahattan plot
  # It can be 1 SNP (e.g. highlighted_SNP="Chr4_10420088") or several SNPs, passed as a vector
  # (e.g. highlighted_SNP=c("Chr4_10420088","Chr5_112000"). No SNP highlighted by default
  
  # Import GEMMA output file
  gwas.results <- read.delim(path.file, sep="\t")
  
  # Plot QQ plot (need to precise the package as lattice has a similar function
  #qqman::qq(gwas.results$P, main=file.name)
  
  # One can select SNPs above the Bonferroni corrected p-value threshold
  # by using the argument "bonferroni"
  if(threshold_pvalue == "bonferroni"){
    # Calculate Bonferroni threshold with risk 5%
    ## Get total number of SNPs
    nb_snps <- dim(gwas.results)[[1]]
    
    ## Calculate Bonferroni corrected P-value threshold
    bonferroni_threshold <- 0.05/nb_snps
    
    threshold_pvalue <- bonferroni_threshold
  } else {
    # In case the variable was entered as string and is not "bonferroni"
    # convert to numeric. Set to 0 by default if user does not want any threshold
    threshold_pvalue <-  as.numeric(threshold_pvalue)
  }
  
  # Get positions of the chromosome with SNPs having a -log(P) > 5
  gwas_significant <- subset(gwas.results, P < threshold_pvalue)
  
  # Default p-value threshold line commonly used in GWAS -> -log10(5e-8) => red line. 
  # Set genomewideline to False has it makes little sense for Arabidopsis genome
  
  # suggestive line = Bonferroni corrected P-value threshold => blue line
  
  # Plot manhattan plot
  manhattan(gwas.results, highlight=highlighted_SNP, main=file.name, suggestiveline = -log10(threshold_pvalue), genomewideline = FALSE)
  
  #Check if dataframe is not empty (no SNPs above threshold value
  if(dim(gwas_significant)[[1]] != 0){ 
    return(gwas_significant)
  }
}



# Import Bismark data -----------------------------------------------------


# Get proper order of the accession based on library name
# For this, upload the list of files to be upload from Bismark coverage files,
# extract their library name and use it to merge to the df_accessions and keep 
# the ordering of the bismark files
# Note that the name of the fastq files should be the same than the name of the library
order_df_accessions <- function(df){
  # I use here CpG coverage report file but CHG or CHH could work as well
  list_files_CpG <- list.files(path=path_bismark_files, pattern="CpG_report_only_chr.txt")
  
  # Order elements of the list so that it runs the same on Linux and Windows
  list_files_CpG <- as.data.frame(stringr::str_sort(list_files_CpG))
  names(list_files_CpG)[1] <- "library"
  list_files_CpG$library <- as.character(list_files_CpG$library)
  
  # Get only the libray part of the string
  list_files_CpG$library <- unlist(lapply(list_files_CpG$library, function(x) strsplit(x, "_bismark")[[1]][[1]]))
  list_files_CpG$library <- as.factor(list_files_CpG$library )
  
  # Get new df_accessions dataframe but this time with the order defined in list_files_CpG
  df_accessions <- merge(list_files_CpG, df_accessions, by="library", sort=FALSE)
}


# This function can import all contexts
# Use cytosine report files from Bismark. Faster than using bam files but
# requires to generate the cytosine reports beforehand

# The file will be gunziped if compressed so that the function won't work if the script was already run as 
# there will be in the directory compressed (txt.gz) and uncompressed (txt) files, both being recognized
# by the list.file() function. Per default, the run_all_bismark.sh script generates uncompressed txt file
# so this function should work if the user does not compress the file => add a check code to avoid the bug

import_bismark_cytosine_report <- function(path_bismark_files, list_DB_paths, list_samples, list_treatments){
  
  # Check if compressed files are present in the directory
  compressed_files <- list.files(path=path_bismark_files, pattern="_report_only_chr.txt.gz", full.names=TRUE)
  
  if(length(compressed_files) != 0){
    stop("Compressed *report_only_chr.txt.gz files are present in the directory. Uncompress first the files and relaunch the analysis.")
  }

  # Get names files of the directory to get all cytosine report files per context
  # Order elements of files with str_sort (which is default in Windows R but not in Linux R)
  list_files_CpG <- list.files(path=path_bismark_files, pattern="CpG_report_only_chr.txt", full.names=TRUE)
  list_files_CpG <- stringr::str_sort(list_files_CpG)
  list_files_CpG <- as.list(list_files_CpG)
  
  list_files_CHG <- list.files(path=path_bismark_files, pattern="CHG_report_only_chr.txt", full.names=TRUE)
  list_files_CHG <- stringr::str_sort(list_files_CHG)
  list_files_CHG <- as.list(list_files_CHG)
  
  list_files_CHH <- list.files(path=path_bismark_files, pattern="CHH_report_only_chr.txt", full.names=TRUE)
  list_files_CHH <- stringr::str_sort(list_files_CHH)
  list_files_CHH <- as.list(list_files_CHH)
  
  
  list_files_all <- list(list_files_CpG, list_files_CHG, list_files_CHH)
  
  variance <- var(as.numeric(lapply(list_files_all, function(x) length(x))))
  
  if (variance != 0){
    stop("The number of *report_only_chr.txt is not equal for each context")
  } else if (length(list_samples) != length(list_files_all[[1]])){
    stop("The number of samples in list_samples is different than the number of *report_only_chr.txt files")
  }
  
  # Create methylRawListDB objects for each context
  all_CpG <- methRead(list_files_CpG, pipeline = "bismarkCytosineReport",
                      sample.id = list_samples,
                      assembly="TAIR10",
                      context = "CpG",
                      treatment = list_treatments,
                      mincov = 0,
                      dbtype = "tabix",
                      header = FALSE,
                      dbdir=list_DB_paths[[1]])
  
  all_CHG <- methRead(list_files_CHG, pipeline = "bismarkCytosineReport",
                      sample.id = list_samples,
                      assembly="TAIR10",
                      context = "CHG",
                      treatment = list_treatments,
                      mincov = 0,
                      dbtype = "tabix",
                      header = FALSE,
                      dbdir=list_DB_paths[[2]])
  

  all_CHH <- methRead(list_files_CHH, pipeline = "bismarkCytosineReport",
                      sample.id = list_samples,
                      assembly="TAIR10",
                      context = "CHH",
                      treatment = list_treatments,
                      mincov = 0,
                      dbtype = "tabix",
                      header = FALSE,
                      dbdir=list_DB_paths[[3]])
  
  list_methylRawList_raw <- list(all_CpG, all_CHG, all_CHH)

  # Delete all uncompressed files
  list_files_uncompressed <- lapply(list_DB_paths, function(i) list.files(path=i, pattern=".txt$", full.names=TRUE))
  
  # Delete files if there is one (list.files returns character(0) if no result from the search)
  lapply(list_files_uncompressed, function(i) if(!identical(i, character(0))) {file.remove(i)})
}


# Filter dataset ----------------------------------------------------------

# Filter samples based on read coverage
# The code below filters a methylRawList and discards bases that have coverage below 2X and also discards the bases that have more than 99.9th percentile of coverage in each sample.
# Since the objects are methylRawListDB objects, the output of filtering will also be generated as external compressed files in the directory of the parent file

# Return a list of methylRawListDB objects filtered according to wanted parameters
filter_methylRawList <- function(list_methylRawList, mincov=2){
  list_methylRawList_filtered <- lapply(list_methylRawList, 
                                        filterByCoverage, 
                                        lo.count=mincov,
                                        lo.perc=NULL,
                                        hi.count=NULL,
                                        hi.perc=99.9)
}


# Merge dataset ----------------------------------------------------------

# Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. 
# Use destrand=TRUE for CpG as it recovers about 7X more sites than if FALSE is used (personal observation).
# This provides better coverage, but only advised when looking at CpG methylation 
# (for CpH methylation this will cause wrong results in subsequent analyses).
# Note that if the methylRawList object is DB, the unite function will generate 
# a methylBaseDB object per default and the flat file will be in the same directory 
# than the methylRawListDB file. The unite function will generate new destrand files before making
# the final methylBase file.

merged_methylRawList <- function(list_methylRawLists, suffix=NULL){

  merged_CpG <- methylKit::unite(list_methylRawLists[[1]], destrand=TRUE, save.db=TRUE, suffix=suffix)
  
  merged_CHG <- methylKit::unite(list_methylRawLists[[2]], destrand=FALSE, save.db=TRUE, suffix=suffix)
  
  merged_CHH <- methylKit::unite(list_methylRawLists[[3]], destrand=FALSE, save.db=TRUE, suffix=suffix)
  
  # Make lists of objects
  list_methylBases <- list(merged_CpG, merged_CHG, merged_CHH)
}


# Pool methylBase -----------------------------------------------------------

# If the samples needed to be pool by treatment, pool() can be used. The function requires
# a sample vector containing characters (follows the order of the treatment).

pool_methylBase <- function(list_methylBases, list_samples, suffix=NULL){
  
  # Check length of list
  if (length(list_methylBases) < 4){
    warning("Less than 4 methylation contexts were provided")
  } else {
    pool_CpG <- pool(list_methylBases[[1]], list_samples, save.db=TRUE, suffix=suffix)
    
    pool_CHG <- pool(list_methylBases[[2]], list_samples, save.db=TRUE, suffix=suffix)
    
    pool_CHH <- pool(list_methylBases[[3]], list_samples, save.db=TRUE, suffix=suffix)

    # Make lists of objects
    list_methylBases <- list(pool_CpG, pool_CHG, pool_CHH)
  }
}

# Load datasets -----------------------------------------------------------

# Load methylRawListDB objects. Load per default raw data (if no type precised)
# Can load whatever file as long as proper type is explicited
load_methylRawListDB <- function(list_DB_paths, type="raw", list_samples, list_treatments, mincov){
  
  # Get the list of list of files and use function get_list_files to do so
  # Note that type is per default raw (no suffix after accession name) but can also be any 
  # suffix as for example 'TEs' (Til-2_TEs.txt.bgz) or 'genes' (Til-2_genes.txt.bgz)
  list_list_files <- lapply(list_DB_paths, get_list_files, list_samples=list_samples, type=type)
  
  # Checkup to be sure the number of samples are matching the number of files
  for (i in length(list_list_files)){
    if (length(list_list_files[[i]]) < length(list_samples)){
      stop("Files are missing in the database directory. Check if the number of elements in list_samples matches with the number of DB files in database directory")
    }
  }
  
  # Check length of list
  if (length(list_list_files) < 3){
    warning("Less than 4 methylation contexts were provided")
  } else {
    # Upload objects into their respective databases and generate methylRawListDB objects
    all_CpG <- methRead(list_list_files[[1]],
                        sample.id = list_samples,
                        assembly="TAIR10",
                        context = "CpG",
                        treatment = list_treatments,
                        mincov = 0,
                        dbtype = "tabix",
                        dbdir=list_DB_paths[[1]])
    
    all_CHG <- methRead(list_list_files[[2]],
                        sample.id = list_samples,
                        assembly="TAIR10",
                        context = "CHG",
                        treatment = list_treatments,
                        mincov = 0,
                        dbtype = "tabix",
                        dbdir=list_DB_paths[[2]])
    
    all_CHH <- methRead(list_list_files[[3]],
                        sample.id = list_samples,
                        assembly="TAIR10",
                        context = "CHH",
                        treatment = list_treatments,
                        mincov = 0,
                        dbtype = "tabix",
                        dbdir=list_DB_paths[[3]])

    # Make lists of objects
    list_methylRawList <- list(all_CpG, all_CHG, all_CHH)
  }
}



# Load methylBaseDB objects using a prefix to identify the files to load.

load_methylBaseDB <- function(list_DB_paths, list_samples, list_treatments, suffix=NULL){
  
  # The function to upload methylBaseDB objects, no exported functions are available but a function exists
  # https://groups.google.com/forum/#!searchin/methylkit_discussion/methylBaseDB%7Csort:date/methylkit_discussion/xH6JTmmdv_8/E1F-nRHbBQAJ
  # methylKit:::readMethylBaseDB
  
  # list stuck with methylKit:::readMethylBaseDB => invalid object for slot "sample.ids" in class 
  # "methylBaseDB" got class "list", should be or extend class "character".
  if (is.list(list_samples)){
    list_samples <- unlist(as.vector(list_samples))
  }
  
  # Get path to each flat file. If suffix is added, file names are adjusted
  # Note that for tile file, add properly the whole suffix (e.g. 'tiled_300_100')
  path_CpG <- paste(list_DB_paths[[1]], paste("methylBase_", suffix, ".txt.bgz", sep=""), sep="/")
  path_CHG <- paste(list_DB_paths[[2]], paste("methylBase_", suffix, ".txt.bgz", sep=""), sep="/")
  path_CHH <- paste(list_DB_paths[[3]], paste("methylBase_", suffix, ".txt.bgz", sep=""), sep="/")

  # Check if all files to be load exists before trying to load
  if (file.exists(path_CpG) && file.exists(path_CHG) && file.exists(path_CHH)){
    
    # Load data
    merged_CpG <- methylKit:::readMethylBaseDB(dbpath = path_CpG,
                                               dbtype = "tabix", sample.ids = list_samples,
                                               assembly = "TAIR10", context = "CpG",
                                               resolution = "base", treatment = list_treatments,
                                               destranded = TRUE)
    
    merged_CHG <- methylKit:::readMethylBaseDB(dbpath = path_CHG,
                                               dbtype = "tabix", sample.ids = list_samples,
                                               assembly = "TAIR10", context = "CHG",
                                               resolution = "base", treatment = list_treatments,
                                               destranded = FALSE)
    
    merged_CHH <- methylKit:::readMethylBaseDB(dbpath = path_CHH,
                                               dbtype = "tabix", sample.ids = list_samples,
                                               assembly = "TAIR10", context = "CHH",
                                               resolution = "base", treatment = list_treatments,
                                               destranded = FALSE)
    
    # Make lists of objects
    list_methylBases <- list(merged_CpG, merged_CHG, merged_CHH)
  } else {
    warning("One of the methylBase object is not existing or suffix is incorrect, check and create object if needed")
  }
  
}


load_methylDiffDB <- function(list_DB_paths, list_samples, list_treatments, suffix=NULL){
  
  # The function to upload methylBaseDB objects, no exported functions are available but a function exists
  # https://groups.google.com/forum/#!searchin/methylkit_discussion/methylBaseDB%7Csort:date/methylkit_discussion/xH6JTmmdv_8/E1F-nRHbBQAJ
  # methylKit:::readMethylBaseDB
  
  # list stuck with methylKit:::readMethylBaseDB => invalid object for slot "sample.ids" in class 
  # "methylBaseDB" got class "list", should be or extend class "character".
  if (is.list(list_samples)){
    list_samples <- unlist(as.vector(list_samples))
  }
  
  # Get path to each flat file. If suffix is added, file names are adjusted
  # Note that for tile file, add properly the whole suffix (e.g. 'tiled_300_100')
  path_CpG <- paste(list_DB_paths[[1]], paste("methylDiff_", suffix, ".txt.bgz", sep=""), sep="/")
  path_CHG <- paste(list_DB_paths[[2]], paste("methylDiff_", suffix, ".txt.bgz", sep=""), sep="/")
  path_CHH <- paste(list_DB_paths[[3]], paste("methylDiff_", suffix, ".txt.bgz", sep=""), sep="/")

  # Load data
  if (file.exists(path_CpG) && file.exists(path_CHG) && file.exists(path_CHH)){
    merged_CpG <- methylKit:::readMethylDiffDB(dbpath = path_CpG,
                                               dbtype = "tabix", sample.ids = list_samples,
                                               assembly = "TAIR10", context = "CpG",
                                               resolution = "base", treatment = list_treatments,
                                               destranded = FALSE)
    
    merged_CHG <- methylKit:::readMethylDiffDB(dbpath = path_CHG,
                                               dbtype = "tabix", sample.ids = list_samples,
                                               assembly = "TAIR10", context = "CHG",
                                               resolution = "base", treatment = list_treatments,
                                               destranded = FALSE)
    
    merged_CHH <- methylKit:::readMethylDiffDB(dbpath = path_CHH,
                                               dbtype = "tabix", sample.ids = list_samples,
                                               assembly = "TAIR10", context = "CHH",
                                               resolution = "base", treatment = list_treatments,
                                               destranded = FALSE)
    
    # Make lists of objects
    list_methylDiffs <- list(merged_CpG, merged_CHG, merged_CHH)
  } else {
    warning("One of the methylDiff object is not existing or suffix is incorrect, check and create object if needed")
  }
}


# Subset data by regions ---------------------------------------------------


# Test to retrieve subset using individual methyl object instead of lists. This is because
# it becomes too heavy to upload the objects in lists when many. Let's then performed the operation 
# on individual methylRawDB object
# With the new version of methylKit (v1.11.1), the creation of DB files let the uncompressed file
# in the directory. I need to delete them through this function
subset_methylObject <- function(list_methylRawObject, list_DB_paths, 
                                path_bed, region, type=c("methylRaw","methylBase","methylDiff")){
  
  # Import bed files (create a Grange object)
  coordinate_regions <- import.bed(path_bed)

  # Subset list and create DB
  # To do overlap, I need to convert the methylRawListDB in methylRawList
  subset_function <- function(methylObject, type, name_DB_file){
    # Check if DB object already exists. If exists, do nothing
    if (!file.exists(name_DB_file)){
      subset_methylObject_RAM <- DB_to_RAM_conversion(methylObject, type)
      regionCounts(subset_methylObject_RAM, coordinate_regions, save.db = TRUE, 
                 dbdir = list_DB_paths[[i]], suffix = region)
      rm(subset_methylObject_RAM)
      gc()
    }
  }

  # Check what type of list of methylObjects is given and subset data accordingly
  if (type == "methylRaw"){
    for (i in seq_along(list_methylRawObject)){
      for (j in seq_along(list_methylRawObject[[i]])){
        methylObject <- unlist(list_methylRawObject[[i]][[j]])
        name_accession <- methylObject@sample.id
        name_DB_file <- paste(list_DB_paths[[i]], "/", name_accession, "_", region, ".txt.bgz", sep="")
        name_DB_file_unziped <- paste(list_DB_paths[[i]], "/", name_accession, "_", region, ".txt", sep="")
        subset_function(methylObject, type, name_DB_file)
        rm(methylObject)
        if (file.exists(name_DB_file_unziped)){
          file.remove(name_DB_file_unziped)
          gc()
        }
      }
    }  
  } else if (type == "methylBase"){
      for (i in seq_along(list_methylRawObject)){
        methylObject <- unlist(list_methylRawObject[[i]])
        name_DB_file <- paste(list_DB_paths[[i]], "/methylBase_", region, ".txt.bgz", sep="")
        name_DB_file_unziped <- paste(list_DB_paths[[i]], "/methylBase_", region, ".txt", sep="")
        subset_function(methylObject, type, name_DB_file)
        if (file.exists(name_DB_file_unziped)){
          file.remove(name_DB_file_unziped)
          gc()
        }
      }
  }
  
}

# Gene annotation ---------------------------------------------------------

# Note that the bed files should contain 10, 11, and 12th columns from BED12 format.
# To do this, check tutorial https://groups.google.com/forum/#!msg/methylkit_discussion/WC6ZM8GY0kc/rvtqrEBIBAAJ
# Basically, one needs to download GTF file from ENSEMBL, convert it to GenePred, then to BED
# Scripts to do so are gtfToGenePred and genePredToBed (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/ )
annotateWithGeneParts_batch <- function(list_objects, path_to_bed){

  gene_obj <- readTranscriptFeatures(path_to_bed)

  list_granges <- lapply(list_objects, "GRanges")
  
  list_annotations <- lapply(list_granges, annotateWithGeneParts, gene_obj)
  
}


# DMR functions ---------------------------------------------

# Function to merge DMRs (should be in methylDiff format) and get an average 
# methylation value for overlapping DMRs. The output dataframe is in BED format.

merge_DMRs <- function(methylDiff_object){
  require("valr")
  require("tidyr")
  require("GenomicRanges")
  
  # Convert methylDiff into GRange
  gr <- as(methylDiff_object, "GRanges")
  
  # Convert GR to dataframe in bed format
  # Rename seqnames to chrom and convert to char otherwise error from valr  
  # Note Grange is 1-based and bed is 0-based so one needs to subtract 1 to start position
  df <- data.frame(chrom=as.character(seqnames(gr)),
                   start=as.integer(start(gr)-1),
                   end=as.integer(end(gr)),
                   name=c(rep(".", length(gr))),
                   score=as.double(gr$meth.diff))
  
  # Merge first DMRs and keep average methyl.diff value.
  df_merged <- df %>% valr::bed_merge(score=mean(score)) 
  
  # Now add name to each DMR
  len <- dim(df_merged)[[1]]
  name_DMRs <- paste("DMR", seq(1, len, 1), sep="")
  
  df_merged$name <- name_DMRs
  
  # Reorder columns
  df_merged_ordered <-   df_merged %>% 
    dplyr::select(chrom, start, end, name, score)
  
  return(df_merged_ordered)
}

# GO enrichment analysis functions ----


#goterms <- Term(GOTERM)
#term2name <- data.frame("GOID"=names(goterms),"term"=goterms )
#saveRDS(term2name, "data/term2name.rds")

# Objects needed for GO analysis
GO_analysis_data <-  readRDS("data/GO_analysis_data.rds")
term2name <- readRDS("data/term2name.rds")


# Function to perform ego
# Function to generate ego analysis
ego_analysis <- function(vector_geneID){
  
  require(tidyverse)
  require(clusterProfiler)
  require(DOSE)
  require(enrichplot)
  require(xlsx)
  require(readxl)
  require(AnnotationDbi)
  require(GO.db)
  
  # Check if all objects are loaded
  if(!exists("GO_analysis_data")){stop("GO_analysis_data object not loaded")}
  if(!exists("term2name")){stop("term2name object not loaded")}
  
  ego_BP <-enricher(
    gene=vector_geneID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=GO_analysis_data,
    TERM2NAME = term2name$BP
  )
  
  ego_CC <-enricher(
    gene=vector_geneID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=GO_analysis_data,
    TERM2NAME = term2name$CC
  )
  
  
  ego_MF <-enricher(
    gene=vector_geneID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=GO_analysis_data,
    TERM2NAME = term2name$MF
  )
  
  
  ego_list <- list(ego_BP = ego_BP, ego_CC = ego_CC, ego_MF = ego_MF)
  
  # Create variable Enrichment factor (Count / number of gene in the given GO)
  ego_list2 <- lapply(ego_list, function(x) mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))))
  
  # Create variable Fold enrichment (GeneRatio / BgRatio)
  ego_final <- lapply(ego_list2, function(x) mutate(x, FoldEnrich = parse_ratio(GeneRatio) / parse_ratio(BgRatio)))
  
  return(ego_final)
  
}

