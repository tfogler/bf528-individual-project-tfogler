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

Sample Report
