# Purpose: Plot the log2FC bar plot for the top 20 transcripts
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(ggplot2)
library(tidyverse)
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
    
    # Top 20 log2FC 
    top_transcripts_lfc_bar_plot <- plot_top_lfc_bar(res, top_n=20)
    ggsave(paste0(plot_dir, "/", "lfc_bar_plot.top_transcripts.", analysis_name, ".", res_name, ".png"), 
            top_transcripts_lfc_bar_plot, width=14, height=12, units="in", dpi=350)

    # Top 20 log2FC - Shrinkage
    top_transcripts_lfc_bar_plot_shrink <- plot_top_lfc_bar(res.shrink, top_n=20)
    ggsave(paste0(plot_dir, "/", "lfc_bar_plot.top_transcripts.shrink.", analysis_name, ".", res_name, ".png"), 
            top_transcripts_lfc_bar_plot_shrink, width=14, height=12, units="in", dpi=350)
}

# Close the log file
sink()
sink(type="message")