# main.R

#Imports

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install(version = "3.18")
# # Bioconductor version '3.19' requires R version '4.4'; 
# # use `BiocManager::install(version = '3.18')` with R version '4.3'
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
# BiocManager::install(c("DESeq2", "fgsea"))

libs <- c("tidyverse", "BiocManager", "DESeq2", "fgsea", "biomaRt", "RColorBrewer")

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

#' Read Data
#' Remove all zero variance genes
#' 
#' @param filename full file path as a string of the counts file 
#'
#' @examples counts_df <- read_data("/path/to/counts/verse_counts.tsv")
read_data <- function(filename, ...) {
  f <- read_tsv(filename)
  var <- apply(f[-1],
                         MARGIN = 1,
                         FUN = var,
                         na.rm = F)
  f <- f[var > 0,]
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
  
  return(dds)
}


#' Perform and plot PCA using processed data.
#'
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs.
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(dds, title="") {
  # get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  # Perform PCA
  pca_result <- prcomp(t(norm_counts), scale = TRUE)
  # Principal components (loadings)
  principal_components <- pca_result$rotation
  # Data in principal component space (scores)
  scores <- pca_result$x
  
  pca_df <- as.data.frame(pca_result$x)
  pca_df$Sample <- rownames(pca_df)  # Add sample names as a column
  
  ggplot(pca_df,
        mapping=aes(`PC1`, `PC2`,
                    colour = `Sample`,
                    label = `Sample`
        )
  ) +
    scale_color_brewer(palette = 'Set2') +
    geom_point(size=6) +
    labs(title = title,
         x = "Principal Component 1",
         y = "Principal Component 2") %>% 
  return(pc_plot)
}

#' Acquire gene list
#' 
acquire_gene_list <- function(df) {
  mart <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- df$genes
  df$id <- NA
  G_list <- biomaRt::getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                            "entrezgene", "description"),values=genes,mart= mart)
  return(G_list)
}


#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  # make ranked list for fgsea
  # start with reading id2gene list
  id2gene <- read_delim(file=id2gene_path, delim='\t', col_names=c('ensembl_id', 'gene_symbol'))
  # assign a joined tibble with mgi genes ids to labeled_res & sort by log2FC descending
  full_join(labeled_results, id2gene, by=join_by(gene == ensembl_id)) %>%
    arrange(desc(log2FoldChange)) %>%
    filter(!is.na(log2FoldChange)) -> labeled_results
  # subset out the l2FC col and add names from mgi symbols b4 returning
  labeled_results$log2FoldChange %>% `names<-`(labeled_results$gene_symbol) %>%
    return()
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- gmtPathways(gmt_file_path)
  fgsea(pathways, rnk_list, min_size, max_size) %>%
    return()
}

#' Function to display fgsea results as tibble sorted by padj
#' (ascending)
#' 
#' 
#' @examples display_fgsea(fgsea_results)
display_fgsea <- function(fgsea_results) {
  fgsea_results %>%
    arrange(padj) %>%
    head(n=10L) %>%
    return()
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  top_fgsea_results <- fgsea_results %>%
    arrange(NES) %>%
    dplyr::slice(c(1:num_paths, (n()-num_paths):n()))
  `levels<-`(top_fgsea_results$NES, top_fgsea_results$NES)
  ggplot(top_fgsea_results, mapping=aes(`NES`, substr(`pathway`, 10, 19))) +
    geom_col(aes(fill = `NES` > 0)) +
    # scale_color_manual(values=c("#0000ff", "red")) +
    labs(title="fgsea results for Hallmark MSigDB gene pathways") %>%
    return()
}


# Draw Volcano Plot in ggplot2
#'
#' @param dataf The loaded data frame.
#' @param x_name The column name to plot on the x-axis
#' @param y_name The column name to plot on the y-axis
#' @param slider A negative integer value representing the magnitude of
#' p-adjusted values to color. Most of our data will be between -1 and -300.
#' @param color1 One of the colors for the points.
#' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
#'
#' @return A ggplot object of a volcano plot
#'
#' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
#' 
volcano_plot <-
  function(dataf, x_name, y_name, alpha) {
    # Filter out NA rows from selected columns
    dataf <- dataf[complete.cases(dataf[, c(x_name, y_name)]), ]
    
    # get top 10 significant DE genes
    top10_genes <- dataf[base::order(dataf[[y_name]]), ][1:10, ]
    
    # GGplot2 Volcano Plot
    volcano <- ggplot2::ggplot(dataf, mapping = aes(x = .data[[x_name]],
                                                   y = -log10(.data[[y_name]]),
                                                   label = row.names(dataf),
                                                   color = (.data[[y_name]]) < alpha)
    ) + ggplot2::geom_text(
      data = top10_genes,
      aes(label = row.names(top10_genes)),
      vjust = -1,
      hjust = 1.2,
      check_overlap = TRUE
    ) + ggplot2::geom_point(
    ) + ggplot2::theme(legend.position = "bottom"
    ) + scale_color_brewer(palette = "Set2"
    ) + ggplot2::labs(title="Plot of DE Results"
    ) #+ ggplot2::coord_fixed(ratio = 0.045)
    
    volcano <- update_labels(volcano, list(x = x_name, y = y_name))
    
    return(volcano)
  }