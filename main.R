# main.R

#Imports

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
# Bioconductor version '3.19' requires R version '4.4'; 
# use `BiocManager::install(version = '3.16')` with R version '4.2'
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

library(tidyverse)
library(DESeq2)

#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x 1+m) tibble with a 'gene' column followed by
#' sample names as column names.
#'
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#'
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(filename){

  return (read_tsv(filename))
}


