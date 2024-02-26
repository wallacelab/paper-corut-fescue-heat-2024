# Purpose: Plot the p-value histogram for each DESeq2 analysis
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(reshape2)
library(EnhancedVolcano)

# Source the functions script
source(snakemake@params[['funcs_script']])

# Load the dds object
dds <- readRDS(snakemake@input[['dds_subset']])

# Extract analysis name
analysis_name <- gsub("dds_(.*).rds", "\\1", basename(snakemake@input[['dds_subset']]))

padj_threshold <- snakemake@params[["padj_threshold"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]
p_cutoff <- snakemake@params[["p_cutoff"]]
q_cutoff <- snakemake@params[["q_cutoff"]]

# Define directories
plot_dir <- paste0(snakemake@output[["plots_dir"]])
if (!dir.exists(paste0(plot_dir))) {
    dir.create(paste0(plot_dir), recursive = TRUE, showWarnings = FALSE)
}

# Analysis
result_names <- resultsNames(dds)[!resultsNames(dds) == "Intercept"]

for (res_name in result_names) {
    # Load the results
    res_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds")
    res <- readRDS(res_file_path)

    print(res)
    
    # Load the shrinked results
    res_shrink_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds")
    res.shrink <- readRDS(res_shrink_file_path)

    print(res.shrink)

    # P-value Histogram - Shrinkage
    png(paste0(plot_dir, "/", "pval_hist.shrink.", analysis_name, ".", res_name, ".png") , width=14, height=10)
    hist(res.shrink$pvalue[res.shrink$baseMean > 1], breaks = 0:20/20, 
         col = "grey50", border = "white", xlab = "", ylab = "", 
         main = paste("P-value Histogram (Lfc Shrinked) -", analysis_name, res_name))
    dev.off()

    # P-value Histogram - Unshrinked
    png(paste0(plot_dir, "/", "pval_hist.", analysis_name, ".", res_name, ".png") , width=14, height=10)
    hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, 
         col = "grey50", border = "white", xlab = "", ylab = "", 
         main = paste("P-value Histogram -", analysis_name, res_name))
    dev.off()
}

# Close the log file
sink(type="message")
sink()