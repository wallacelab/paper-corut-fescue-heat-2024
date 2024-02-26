# Purpose: Run CEMiTool co-expression analysis

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(CEMiTool)

# Read inputs
dds <- readRDS(snakemake@input$dds_subset)
gmt_in <- CEMiTool::read_gmt(snakemake@input$gmt)
int_df <- read.csv(snakemake@input$interactions, stringsAsFactors = FALSE)

# Get cemitool parameters from Snakemake's params
verbose <- snakemake@params$verbose
plot <- snakemake@params$plot
apply_vst <- snakemake@params$apply_vst
filter <- snakemake@params$filter
count_filter <- snakemake@params$count_filter
min_counts <- snakemake@params$min_counts
min_samples <- snakemake@params$min_samples
network_type <- snakemake@params$network_type
cor_method <- snakemake@params$cor_method

# Get analysis name from Snakemake's params
analysis_name <- snakemake@params$analysis$name

# Get output directory from Snakemake's params
cemitool_outdir <- snakemake@output$cemitool_outdir

# cemitool function
run_cemitool <- function(dds, gmt, int_df, directory, verbose=TRUE, plot=TRUE, apply_vst=FALSE, filter=TRUE, network_type="signed", class = "time_point", cor_method = "spearman", count_filter=TRUE, min_counts=5, min_samples=4) {
      
      # Get counts
      counts <- as.data.frame(counts(dds, normalized=TRUE))

      print("Number of genes before filtering:")
      print(nrow(counts))

      # Filter out genes with low counts
      if (count_filter) {
        counts <- counts[rowSums(counts > min_counts) >= min_samples, ]
      }

      # Check number of genes after filtering
      print("Number of genes after filtering:")
      print(nrow(counts))

      # Get metadata
      dds.meta <- as.data.frame(colData(dds))
      dds.meta <- tibble::rownames_to_column(dds.meta, "SampleName")
      colnames(dds.meta)[colnames(dds.meta) == class] <- "Class"
      
      # Set class levels
      if (class == "time_point") {
        dds.meta$Class <- factor(dds.meta$Class, levels = c("HR0", "HR2", "HR12", "D1", "D3", "W1"))
      } else if (class == "treatment") {
        dds.meta$Class <- factor(dds.meta$Class, levels = c("uninfected", "infected"))
      } else if (class == "condition") {
        dds.meta$Class <- factor(dds.meta$Class, levels = c("NC", "HS"))
      } else {
        dds.meta$Class <- factor(dds.meta$Class)
      }   

      # Run cemitool
      cem <- cemitool(expr=counts, annot=dds.meta, 
                      gmt=gmt, interactions=int_df, 
                      network_type=network_type, 
                      filter = filter, 
                      apply_vst = apply_vst, 
                      verbose=verbose, 
                      plot=plot, 
                      cor_method = "spearman")
      
      # create report as html document
      generate_report(cem, directory = paste0(directory, "/report"), output_format="html_document")
      
      # create diagnosti report as html document
      diagnostic_report(cem, directory = paste0(directory, "/diagnostic"), output_format="html_document")
      
      # write analysis results into files
      write_files(cem, directory=paste0(directory, "/tables"))
      
      # save all plots
      save_plots(cem, "all", directory=paste0(directory, "/plots"))
      
      # return cem object
      return(cem)
}

print("Running cemitool...")

print("Analysis name:")
print(analysis_name)

# Based on analysis name, set the class parameter
if (grepl("inf_vs_uninf", analysis_name)) {
  class <- "treatment"
} else if (grepl("heat_vs_normal", analysis_name)) {
  class <- "condition"
} else {
  class <- "time_point"
}

# Call the function with the class parameter from Snakemake's params
cem <- run_cemitool(dds, gmt_in, int_df, directory= paste0(cemitool_outdir, "/", class, "/"),
                    verbose=verbose, plot=plot, apply_vst=apply_vst, filter=filter, 
                    network_type=network_type, class = class, cor_method = cor_method,
                    count_filter=count_filter)

# Save the cem object
saveRDS(cem, paste0(cemitool_outdir, "/", class, "/", analysis_name, ".cem.", class, ".rds"))

# if not class == "time_point" then run cemitool again with class = "time_point"
if (!class == "time_point") {
  cem_time_point <- run_cemitool(dds, gmt_in, int_df, directory= paste0(cemitool_outdir, "/time_point/"),
                    verbose=verbose, plot=plot, apply_vst=apply_vst, filter=filter, 
                    network_type=network_type, class = "time_point", cor_method = cor_method,
                    count_filter=count_filter)
  
  # Save the cem object
  saveRDS(cem_time_point, paste0(cemitool_outdir, "/time_point/", analysis_name, ".cem.time_point.rds"))
}

# Close the log file
sink()
sink(type="message")