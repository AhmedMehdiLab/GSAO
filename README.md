# GSAO
Gene Set Annotation Overlap

This is a data analysis package allowing the user with a list of genes of interest to find enriched annotations of gene sets from MSigDB C7.

To run:

```
library(GSAO)
genes <- "GENE1 GENE2 GENE3"       # a gene list separated by spaces or commas
input <- process_input_text(genes) 
results <- compute(input)
```

Statistically enriched annotations are stored in a `tidyverse` `tibble`, and can be viewed with:

```
results$stats
```

Roxygen documentation is available for all functions.
