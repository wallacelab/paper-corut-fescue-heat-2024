# Purpose: Plot the expression heatmap for the DEG analysis
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)

# Source the functions script
source(snakemake@params[['funcs_script']])

# Load the dds and vsd objects
dds <- readRDS(snakemake@input[['dds_subset']])
vsd <- readRDS(snakemake@input[['vsd_subset']])

# Extract analysis name
analysis_name <- snakemake@params$analysis$name
coef <- snakemake@params$analysis$coef_of_interest

# Set thresholds
padj_threshold <- snakemake@params[["padj_threshold"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]
top_n <- snakemake@params[["top_n"]]

# Define directories
plot_dir <- snakemake@output[["plots_dir"]]
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}

# For each result_name, load the corresponding res and res.shrink
result_names <- resultsNames(dds)[!resultsNames(dds) == "Intercept"]

result_names <- result_names[grep(coef, result_names)]

for (res_name in result_names) {
    # Extract the result for this specific comparison
    res_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds")
    print(paste("Loading", res_file_path))
    res <- readRDS(res_file_path)

    res_shrink_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds")
    print(paste("Loading", res_shrink_file_path))  # Corrected to res_shrink_file_path for accurate message
    res.shrink <- readRDS(res_shrink_file_path)

    # Get the time point as string based on res_name
    # If res_name contains "HR2", the time point is "HR2"
    # If res_name contains "HR12", the time point is "HR12"
    # If res_name contains "D1", the time point is "D1"
    # If res_name contains "D3", the time point is "D3"
    # If res_name contains "W1", the time point is "W1"
    if (grepl("HR2", res_name)) {
        time_point <- "HR2"
    } else if (grepl("HR12", res_name)) {
        time_point <- "HR12"
    } else if (grepl("D1", res_name)) {
        time_point <- "D1"
    } else if (grepl("D3", res_name)) {
        time_point <- "D3"
    } else if (grepl("W1", res_name)) {
        time_point <- "W1"
    } else {
        stop("Could not determine time point for res_name", res_name)
    }

    print(paste("Plotting", res_name, "with time point", time_point))

    png(paste0(plot_dir, "/deg.exp_heatmap.", analysis_name, ".", res_name, ".png"), width=16, height=10, units = "in", res=350)
    exp_heatmap <- plot_exp_heatmap(dds, vsd, res, top_n=top_n, time_point=time_point, alternate_color=FALSE)
    dev.off()

    png(paste0(plot_dir, "/deg.exp_heatmap.shrink.", analysis_name, ".", res_name, ".png"), width=16, height=10, units = "in", res=350)
    exp_heatmap_shrink <- plot_exp_heatmap(dds, vsd, res.shrink, top_n=top_n, time_point=time_point, alternate_color=FALSE)
    dev.off()
}

sink()
sink(type="message")
