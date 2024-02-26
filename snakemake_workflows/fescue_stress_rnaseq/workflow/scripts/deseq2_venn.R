# Purpose: Generate venn diagram of significant DE genes from DESeq2 results
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(VennDetail)

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

signif_de_genes_list <- list()
signif_de_genes_list_shrink <- list()

for (res_name in result_names) {
    res_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds")
    print(paste("Loading", res_file_path))
    res <- readRDS(res_file_path)

    res_shrink_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds")
    print(paste("Loading", res_shrink_file_path))  # Corrected to res_shrink_file_path for accurate message
    res_shrink <- readRDS(res_shrink_file_path)

    # Get significant results
    res_shrink_signif <- res_shrink[res_shrink$padj < padj_threshold & !is.na(res_shrink$padj) & abs(res_shrink$log2FoldChange) >= lfc_threshold, ]

    # Get the significant genes
    res_shrink_signif_genes <- rownames(res_shrink_signif)

    # Change list name to include only the timepoint
    # Time point is the string after "time_point"
    time_point <- gsub(".*time_point", "", res_name)

    # Determine the position for insertion based on the specified conditions
    # The initial position is after the last element currently in the list
    insert_pos <- length(signif_de_genes_list_shrink) + 1

    # Adjust the position based on the specified ordering rules
    if (grepl("HR2", time_point) && !grepl("HR12", time_point)) {
        insert_pos <- 1
    } else if (grepl("HR12", time_point)) {
        # If there is already an "HR2" element, "HR12" should follow it
        if (any(grepl("HR2", names(signif_de_genes_list_shrink)))) {
            insert_pos <- which(grepl("HR2", names(signif_de_genes_list_shrink))) + 1
        } else {
            insert_pos <- 1
        }
    } else if (grepl("D1", time_point)) {
        # If there is already an "HR12" element, "D1" should follow it
        if (any(grepl("HR12", names(signif_de_genes_list_shrink)))) {
            insert_pos <- which(grepl("HR12", names(signif_de_genes_list_shrink))) + 1
        } else if (any(grepl("HR2", names(signif_de_genes_list_shrink)))) {
            insert_pos <- which(grepl("HR2", names(signif_de_genes_list_shrink))) + 1
        } else {
            insert_pos <- 1
        }
    } else if (grepl("D3", time_point)) {
        # Similar logic for "D3", it should follow "D1" if it exists
        if (any(grepl("D1", names(signif_de_genes_list_shrink)))) {
            insert_pos <- which(grepl("D1", names(signif_de_genes_list_shrink))) + 1
        } else if (any(grepl("HR12", names(signif_de_genes_list_shrink)))) {
            insert_pos <- which(grepl("HR12", names(signif_de_genes_list_shrink))) + 1
        } else if (any(grepl("HR2", names(signif_de_genes_list_shrink)))) {
            insert_pos <- which(grepl("HR2", names(signif_de_genes_list_shrink))) + 1
        } else {
            insert_pos <- 1
        }
    } else if (grepl("W1", time_point)) {
        # "W1" should always be last, so its position is correct by default
    }

    # Now insert the elements at the determined position
    signif_de_genes_list_shrink <- append(signif_de_genes_list_shrink, setNames(list(res_shrink_signif_genes), time_point), insert_pos - 1)


    print(res_name)
    print(time_point)
    # signif_de_genes_list_shrink[[time_point]] <- res_shrink_signif_genes
}

venn_details <- venndetail(signif_de_genes_list_shrink)

# Plot the venn diagram and save
png(paste0(plot_dir, "/", analysis_name, ".venn_plot.png"), width = 10, height = 8, units = "in", res = 350)
plot(venn_details, type = "venn", cat.cex = 2)
dev.off()