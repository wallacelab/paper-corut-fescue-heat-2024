# Purpose: Plot the DEG bar plot for each analysis
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggsci)

# Source the functions script
source(snakemake@params[['funcs_script']])

# Load the dds object
dds <- readRDS(snakemake@input[['dds_subset']])

# Extract analysis name
analysis_name <- snakemake@params$analysis$name
coef <- snakemake@params$analysis$coef_of_interest

# Set thresholds
padj_threshold <- snakemake@params[["padj_threshold"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]

# Define directories
plot_dir <- snakemake@output[["plots_dir"]]
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}

# For each result_name, load the corresponding res and res.shrink
result_names <- resultsNames(dds)[!resultsNames(dds) == "Intercept"]

result_names <- result_names[grep(coef, result_names)]

results_list <- list()
results_list_shrink <- list()

for (res_name in result_names) {
    res_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds")
    print(paste("Loading", res_file_path))
    res <- readRDS(res_file_path)

    res_shrink_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds")
    print(paste("Loading", res_shrink_file_path))  # Corrected to res_shrink_file_path for accurate message
    res.shrink <- readRDS(res_shrink_file_path)

    # Determine the position for insertion based on the specified conditions
    # The initial position is after the last element currently in the list
    insert_pos <- length(results_list) + 1

    # Adjust the position based on the specified ordering rules
    if (grepl("HR2", res_name) && !grepl("HR12", res_name)) {
        insert_pos <- 1
    } else if (grepl("HR12", res_name)) {
        # If there is already an "HR2" element, "HR12" should follow it
        if (any(grepl("HR2", names(results_list)))) {
            insert_pos <- which(grepl("HR2", names(results_list))) + 1
        } else {
            insert_pos <- 1
        }
    } else if (grepl("D1", res_name)) {
        # If there is already an "HR12" element, "D1" should follow it
        if (any(grepl("HR12", names(results_list)))) {
            insert_pos <- which(grepl("HR12", names(results_list))) + 1
        } else if (any(grepl("HR2", names(results_list)))) {
            insert_pos <- which(grepl("HR2", names(results_list))) + 1
        } else {
            insert_pos <- 1
        }
    } else if (grepl("D3", res_name)) {
        # Similar logic for "D3", it should follow "D1" if it exists
        if (any(grepl("D1", names(results_list)))) {
            insert_pos <- which(grepl("D1", names(results_list))) + 1
        } else if (any(grepl("HR12", names(results_list)))) {
            insert_pos <- which(grepl("HR12", names(results_list))) + 1
        } else if (any(grepl("HR2", names(results_list)))) {
            insert_pos <- which(grepl("HR2", names(results_list))) + 1
        } else {
            insert_pos <- 1
        }
    } else if (grepl("W1", res_name)) {
        # "W1" should always be last, so its position is correct by default
    }

    # Now insert the elements at the determined position
    results_list <- append(results_list, setNames(list(res), res_name), insert_pos - 1)
    results_list_shrink <- append(results_list_shrink, setNames(list(res.shrink), res_name), insert_pos - 1)

    print(names(results_list))
    print(names(results_list_shrink))
}

# Call the plotting function
deg_bar_plot <- plot_deg_bar(results_list, padj_cutoff=padj_threshold, lfc_cutoff=lfc_threshold)
deg_bar_plot_shrink <- plot_deg_bar(results_list_shrink, padj_cutoff=padj_threshold, lfc_cutoff=lfc_threshold)

# Save the plot
ggsave(paste0(plot_dir, "/", analysis_name, ".deg_bar_plot.png"), deg_bar_plot, width=12, height=10, dpi=350)
ggsave(paste0(plot_dir, "/", analysis_name, ".deg_bar_plot.shrink.png"), deg_bar_plot_shrink, width=12, height=10, dpi=350)

# Close the sink and the file connection
sink(type = "message")
sink()
close(log)