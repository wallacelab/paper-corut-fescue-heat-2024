# Purpose: Import GO IDs and GO terms from topGO output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load libraries
library(tidyverse)
library(stringr)
library(topGO)
library(clusterProfiler)

# Load GO IDs list
term_to_gene <- read.csv(snakemake@input[["go_info"]], sep="\t")
print(head(term_to_gene))

goterms <- as.data.frame(Term(GOTERM))

goterms <- tibble::rownames_to_column(goterms, "GOID") # Add GOID column

term_to_name <- goterms[goterms$GOID %in% term_to_gene$GOID,] # Filter out GO terms that are not in the GO IDs list
print(head(term_to_name))

write_rds(term_to_gene, path = snakemake@output[["term_to_gene"]], compress = "none") # Save GO IDs list

write_rds(term_to_name, path = snakemake@output[["term_to_name"]], compress = "none") # Save GO term list