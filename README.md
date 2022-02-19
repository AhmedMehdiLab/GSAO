# GSAO
<!-- badges: start -->
[![R-CMD-check](https://github.com/AhmedMehdiLab/GSAO/workflows/R-CMD-check/badge.svg)](https://github.com/AhmedMehdiLab/GSAO/actions)
<!-- badges: end -->

Gene Set Annotation Overlap

This is a data analysis package allowing the user with a list of genes of interest to find enriched annotations of gene sets from MSigDB C7.

## Installation
To install this package, run:

``` r
# install.packages("remotes")
remotes::install_github("AhmedMehdiLab/GSAO")
```

## Usage
To use the package, run:

``` r
library(GSAO)
genes <- "GENE1 GENE2 GENE3"       # a gene list separated by spaces or commas
input <- process_input_text(genes) 
results <- compute(input)
```

Statistically enriched annotations are stored in a `tidyverse` `tibble`, and can be viewed with:

``` r
results$stats
```

Roxygen documentation is available for all functions.

## Notes
You will be asked to give the path for MSigDB XML file. Please drag the XML file to R console to get the path or alternatively write the path:

Example:
~/Downloads/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1.xml

To run it for gene list, please provide the gene names or gene names;

``` r
input <- process_input_text("IL1B	EGR3	EGR2	CCR1	PTGS2	SLC25A29	RAB12	CXCL1	FOSB	SFRS15	SGK	FRAT2	TRIB1	NGLY1	THOC7	FPRL1	SFRS2IP	TAX1BP1	FLJ42008	GPR109B	TREM1	MAFB	BTBD14A")

results <- compute(input)
head(results$stats)
head(results$matches)
```

For single-cell RNAseq, please provide object processed through Seurat pipeline (normalization, scaling, clustering etc performed. With also a label “grp” indicating the group information.

``` r
input_Seurat_object <- readRDS(‘Your Path to Seurat RDS File’)
```

with input_Seurat_object$seurat_clusters containing the seurat clusters to perform differential expression analyses needed in downstream analyses.
