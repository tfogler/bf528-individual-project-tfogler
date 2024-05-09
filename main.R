# main.R

#Imports

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install(version = "3.18")
# # Bioconductor version '3.19' requires R version '4.4'; 
# # use `BiocManager::install(version = '3.18')` with R version '4.3'
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
# BiocManager::install(c("DESeq2", "fgsea"))

libs <- c("tidyverse", "BiocManager", "DESeq2", "fgsea")

for (package in libs) {
  suppressPackageStartupMessages(require(package, 
                                         quietly = T, 
                                         character.only = T))
  require(package, character.only = T)
}


#' Load n' trim
#'
#' @param filename full file path as a string of the counts file 
#'
#' @return A _data frame_ with gene names as row names. A tibble will **not** work 
#' with the differential expression packages.
#' 
#' @details As always, we need to load our data and start to shape it into the 
#' form we need for our analysis. Selects only the columns named "gene" we need.
#'
#' @examples counts_df <- load_n_trim("/path/to/counts/verse_counts.tsv")
load_n_trim <- function(filename, ..., row.names) {
  f <- read.delim(filename, row.names=row.names)
  f <- f[,c("CTL_rep1", "CTL_rep2", "CTL_rep3", "KO_rep1", "KO_rep2", "KO_rep3")]
  return(f)
}

#' Perform a DESeq2 analysis of rna seq data
#'
#' @param count_dataframe The data frame of gene names and counts.
#' @param coldata The coldata variable describing the experiment, a dataframe.
#' @param count_filter An arbitrary number of genes each row should contain or 
#' be excluded. DESeq2 suggests 10, but this could be customized while running. 
#' An integer.
#' @param condition_name A string identifying the comparison we are making. It 
#' follows the format "condition_[]_vs_[]". If I wanted to compare day4 and day7 
#' it would be "condition_day4_vs_day7".
#'
#' @return A dataframe of DESeq results. It has a header describing the 
#' condition, and 6 columns with genes as row names. 
#' @details This function is based on the DESeq2 User's Guide. These links describe 
#' the inputs and process we are working with. The output we are looking for comes 
#' from the DESeq2::results() function.
#' https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input
#' https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
#'
#' @examples run_deseq(counts_df, coldata, 10, "condition_day4_vs_day7")
run_deseq <- function(count_dataframe, coldata, count_filter, condition_name) {
  # build DEseq dataset/se object from counts_df and coldata
  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(count_dataframe),
    colData = coldata,
    design = ~ condition
  )
  #filter out counts below filter threshold
  keep <- rowSums(counts(dds)) >= count_filter
  dds <- dds[keep,]
  
  # relevel factors in colDatat of dds
  # set "CTL" to baseline
  dds$condition <- factor(dds$condition) %>% relevel(
    ref = str_extract(condition_name, '[^_]+$')
  )
  # run differential analysis
  dds <- DESeq(dds, quiet=T)
  
  results(dds)
}