# Load libraries
# Purpose: Generate Venn diagrams for significant genes from DESeq2 results for different comparisons and time points
library(tidyverse)
library(DESeq2)
library(VennDetail)

# Open log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

# Function to create file pattern
create_pattern <- function(base_pattern, time_point) {
  return(paste0(base_pattern, time_point, ".shrink.rds"))
}

# Function to generate names for list elements based on the pattern and file names
generate_names_for_list <- function(file_names, base_pattern, time_point) {
  modified_pattern <- gsub("\\.\\*", "|", base_pattern)
  names_for_list <- gsub(paste0(modified_pattern, time_point, "|.shrink.rds"), "", basename(file_names))
  names_for_list <- paste0(names_for_list, time_point)
  print(names_for_list)
  return(names_for_list)
}

# Function to plot Venn diagram
plot_venn_diagram <- function(sig_genes, plot_dir, comparison, time_point) {
  ven <- venndetail(sig_genes)
  png(paste0(plot_dir, "/", comparison, "/", comparison, ".genotype_comp.", time_point, ".venn_plot.png"), width = 10, height = 8,  units = "in", res = 350)
  plot(ven, type = "venn", cat.cex = 2)
  dev.off()
}

# Function to process files for a given time point
process_time_point <- function(time_point, base_pattern, dir, padj_threshold, lfc_threshold, comparison, plot_dir) {
  # Get the files for the time point
  res_shrink_relative <- list.files(dir, 
                            pattern = create_pattern(base_pattern, time_point),
                            recursive = TRUE)
  
  # Get the full path of the files
  res_shrink <- file.path(dir, res_shrink_relative)
  print(res_shrink)
  # Read the files
  res_shrink_list <- lapply(res_shrink, readRDS)

  #
  names_for_list <- generate_names_for_list(res_shrink, base_pattern, time_point)
  
  # Ensure that the length of names matches the length of the list
  if(length(res_shrink_list) == length(names_for_list)) {
    # Set names for the list elements
    names(res_shrink_list) <- names_for_list
  } else {
    stop("Error: The number of elements and the number of names do not match.")
  }

  # Print the list
  print(res_shrink_list)

  # Create a list to store significant genes
  sig_genes <- list()

  # for each element in res_shrink_list
  for (i in 1:length(res_shrink_list)) {
      # Get the name of the element
      name <- names(res_shrink_list[i])
      print(name)

      # Get the deseq2 results for that element
      res <- res_shrink_list[[i]]
      print(res)

      # Get significant results
      res_shrink_signif <- res[res$padj < padj_threshold & !is.na(res$padj) & abs(res$log2FoldChange) >= lfc_threshold, ]
      
      # Get the significant genes
      res_shrink_signif_genes <- rownames(res_shrink_signif)
      
      # Add the significant genes to the sig_genes list
      sig_genes[[name]] <- res_shrink_signif_genes
  }

  # Plot Venn diagram
  plot_venn_diagram(sig_genes, plot_dir, comparison, time_point)

}

# Main function
run_analysis <- function(dds_dir_list, results_dir_list, time_points, base_pattern, padj_threshold, lfc_threshold, comparison, plot_dir) {
  dir <- dirname(results_dir_list[1])
  print(dir)

  # Create plot directory if it does not exist
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Create plot subdirectory if it does not exist
  if (!dir.exists(paste0(plot_dir, "/", comparison))) {
    dir.create(paste0(plot_dir, "/", comparison), recursive = TRUE, showWarnings = FALSE)
  }

  
  for (time_point in time_points) { # for each time point
    print(time_point)
    # Process files for that time point
    process_time_point(time_point, base_pattern, dir, padj_threshold, lfc_threshold, comparison, plot_dir)
  }
}

# Set parameters (could be read from a config file or snakemake parameters)
time_points <- c("HR2", "HR12", "D1", "D3", "W1")
padj_threshold <- snakemake@params[["padj_threshold"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]
plot_dir <- snakemake@output[["plots_dir"]]

# Run analysis
# res.heat_inf_vs_uninf
base_pattern_heat_inf_vs_uninf <- "res.heat_inf_vs_uninf_.*treatmentinfected.time_point"
run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir, time_points, base_pattern_heat_inf_vs_uninf, 
            padj_threshold, lfc_threshold, 
            comparison= "heat_inf_vs_uninf", 
            plot_dir)

# res.heat_uninf_vs_inf
base_pattern_infected_heat_vs_normal <- "res.infected_heat_vs_normal_.*conditionHS.time_point"
run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir, time_points, base_pattern_infected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "infected_heat_vs_normal", 
            plot_dir)

#res.normal_inf_vs_uninf_
base_pattern_normal_inf_vs_uninf <- "res.normal_inf_vs_uninf_.*treatmentinfected.time_point"
run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir, time_points, base_pattern_normal_inf_vs_uninf, 
            padj_threshold, lfc_threshold, 
            comparison= "normal_inf_vs_uninf", 
            plot_dir)

# res.uninfected_heat_vs_normal
base_pattern_uninfected_heat_vs_normal <- "res.uninfected_heat_vs_normal_.*conditionHS.time_point"
run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir, time_points, base_pattern_uninfected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "uninfected_heat_vs_normal", 
            plot_dir)

# Close log file
sink()
sink(type = "message")

