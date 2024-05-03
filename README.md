# bf528-individual-project-tfogler
bf528-individual-project-tfogler created by GitHub Classroom

Main author: tfogler

# RNASeq
## Methods
First, mRNAseq samples were obtained corresponding to the wild-type/control (CTL) and knockout (KO) conditions, 3 CTL and 3 KO. To ensure data quality, we performed initial quality control using FastQC v0.12.1.
Next, we aligned the reads to the gencode human primary assembly genome (GRCh38, release 45) using STAR v2.7.11b with default parameters.
From these alignments, we generated gene counts using VERSE v0.1.5 and the gencode 45 primary assembly GTF.
After combining the counts from all 6 samples, we filtered out genes with zero counts. Then, we performed normalization and differential expression analysis using DESeq2 v1.42.1, comparing the two conditions.
Finally, we performed a functional gene set enrichment analysis using FGSEA v1.28.0 to find functions associated with differentially expressed genes.

## Questions to Address

### Must address the following
Briefly remark on the quality of the sequencing reads and the alignment statistics, make sure to specifically mention the following:
1. Are there any concerning aspects of the quality control of your sequencing reads?
    a. Are there any concerning aspects of the quality control related to alignment?
    b. Based on all of your quality control, will you exclude any samples from further analysis?
2. After generating your counts matrix, perform either a PCA or produce a sample-to-sample distance plot as described in the DESeq2 vignette.
    a. Briefly remark on the plot and what it indicates to you in terms of the experiment
3. After performing DE analysis, choose an appropriate FDR threshold to subset your DE results.
    a. How many genes are significant at your chosen statistical threshold?
4. After performing FGSEA (GSEA) using a ranked list of all genes in the experiment and performing gene set enrichment using your list of statistically significant DE genes, please answer the following questions:
    a, How similar are the results from these two analyses? Are there any notable differences?
    b. Do you expect there to be any differences? If so, why?
    c. What do the results imply about potentional biological functions of the factor of interest?

#### Deliverables

1. Produce either a sample-to-sample distance plot from the counts matrix or a PCA biplot

2. A CSV containing all of the results from your DE analysis

3. A histogram showing the distribution of log2FoldChanges of your DE genes

4. A volcano plot that distinguishes between significant and non-significant DE genes as well as labels up- and downregulated genes

    * Label the top ten most significant genes with their associated gene name / gene symbol.

5. Perform a GSEA (FGSEA) analysis on all the genes discovered in the experiment

    * You may choose a ranking metric you feel is appropriate
    * You may use the C2 canonical pathways database or another of your choice

6. Create a single table / figure that reports the most interesting results
    * Use your list of DE genes at your chosen statistical threshold and perform a basic enrichment analysis using a tool of your choice
    * You may use DAVID, enrichR, or other well-validated gene enrichment tools
    * Create a single table / figure that reports the most interesting results

Sample Report: report.Rmd
