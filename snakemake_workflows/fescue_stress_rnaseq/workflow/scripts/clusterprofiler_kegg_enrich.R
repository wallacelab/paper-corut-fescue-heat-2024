# Purpose: Run KEGG enrichment analysis using clusterProfiler
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(ggsci)
library(R.utils)

R.utils::setOption("clusterProfiler.download.method","wget")

# Source the functions script
source(snakemake@params[['funcs_script']])

# Load the dds object
dds <- readRDS(snakemake@input[['dds_subset']])

# Extract analysis name
analysis_name <- gsub("dds_(.*).rds", "\\1", basename(snakemake@input[['dds_subset']]))

tcp_to_gene_df <- read.table(snakemake@input[["tcp_to_gene"]], sep = "\t", header = TRUE)
cat("tcp_to_gene_df\n")
print(head(tcp_to_gene_df))

# Set thresholds
padj_threshold <- snakemake@params[["padj_threshold"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]
p_cutoff <- snakemake@params[["p_cutoff"]]
q_cutoff <- snakemake@params[["q_cutoff"]]
kegg_organism <- snakemake@params[["kegg_organism"]]

# Define directories
plot_dir <- snakemake@output[["plots_dir"]]
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}


# For each result_name, load the corresponding res and res.shrink
result_names <- resultsNames(dds)[!resultsNames(dds) == "Intercept"]

for (res_name in result_names) {
    res_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds")
    res <- readRDS(res_file_path)
    
    print(res)
    
    res_shrink_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds")
    res.shrink <- readRDS(res_shrink_file_path)
    
    print(res.shrink)

    # Create the directory if doesn't exist
    if (!dir.exists(snakemake@output[["tables_dir"]])) {
        dir.create(snakemake@output[["tables_dir"]], recursive = TRUE, showWarnings = FALSE)
    }
    
    ##### KEGG Enrichment Analysis #####
    KEGG.enrich.res <- extract_and_kegg_enrich(tcp_to_gene_df, res, dds, 
                                                    enrichment_type = "pathway", org = 'ath',
                                                    p_cutoff=p_cutoff, q_cutoff=q_cutoff,  
                                                    padj_threshold = padj_threshold, lfc_threshold = lfc_threshold)
    print(KEGG.enrich.res)
    KEGG.module.enrich.res  <- extract_and_kegg_enrich(tcp_to_gene_df, res, dds, 
                                                            enrichment_type = "module", org = 'ath', 
                                                            p_cutoff=p_cutoff, q_cutoff=q_cutoff, 
                                                            padj_threshold = padj_threshold, lfc_threshold = lfc_threshold)
    print(KEGG.module.enrich.res)   

    write.csv(KEGG.enrich.res, 
            paste0(snakemake@output[["tables_dir"]], "/KEGG_enrichment.", analysis_name, ".", res_name, ".csv"), 
            row.names = TRUE)
    write.csv(KEGG.module.enrich.res,
            paste0(snakemake@output[["tables_dir"]], "/KEGG_module_enrichment.", analysis_name, ".", res_name, ".csv"), 
            row.names = TRUE)

    KEGG.enrich.res.shrink <- extract_and_kegg_enrich(tcp_to_gene_df, res.shrink, dds, 
                                                    enrichment_type = "pathway", org = 'ath',
                                                    p_cutoff=p_cutoff, q_cutoff=q_cutoff,  
                                                    padj_threshold = padj_threshold, lfc_threshold = lfc_threshold)
    print(KEGG.enrich.res.shrink)
    
    KEGG.module.enrich.res.shrink  <- extract_and_kegg_enrich(tcp_to_gene_df, res.shrink, dds, 
                                                            enrichment_type = "module", org = 'ath', 
                                                            p_cutoff=p_cutoff, q_cutoff=q_cutoff, 
                                                            padj_threshold = padj_threshold, lfc_threshold = lfc_threshold)
    print(KEGG.module.enrich.res.shrink)
    write.csv(KEGG.enrich.res.shrink,
            paste0(snakemake@output[["tables_dir"]], "/KEGG_enrichment.", analysis_name, ".", res_name, ".shrink.csv"), 
            row.names = TRUE)
    write.csv(KEGG.module.enrich.res.shrink,
            paste0(snakemake@output[["tables_dir"]], "/KEGG_module_enrichment.", analysis_name, ".", res_name, ".shrink.csv"), 
            row.names = TRUE)
    
    ## Bubble Plot
    kegg_bubble <- generate_kegg_enrich_bubble_plot(KEGG.enrich.res, n_category=20)
    ggsave(paste0(plot_dir, "/", "kegg_bubble.", analysis_name, ".", res_name, ".png"), kegg_bubble, width=12, height=10, dpi=350)

    kegg_module_bubble <- generate_kegg_enrich_bubble_plot(KEGG.module.enrich.res, n_category=20)
    ggsave(paste0(plot_dir, "/", "kegg_module_bubble.", analysis_name, ".", res_name, ".png"), kegg_module_bubble, width=12, height=10, dpi=350)

    kegg_bubble_shrink <- generate_kegg_enrich_bubble_plot(KEGG.enrich.res.shrink, n_category=20)
    ggsave(paste0(plot_dir, "/", "kegg_bubble.shrink.", analysis_name, ".", res_name, ".png"), kegg_bubble_shrink, width=12, height=10, dpi=350)

    kegg_module_bubble_shrink <- generate_kegg_enrich_bubble_plot(KEGG.module.enrich.res.shrink, n_category=20)
    ggsave(paste0(plot_dir, "/", "kegg_module_bubble.shrink.", analysis_name, ".", res_name, ".png"), kegg_module_bubble_shrink, width=12, height=10, dpi=350)
            
}

# Close logging
sink(type="message")
sink()
