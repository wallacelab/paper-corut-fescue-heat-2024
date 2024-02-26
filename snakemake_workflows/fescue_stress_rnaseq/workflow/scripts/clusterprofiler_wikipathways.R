# Purpose: Run clusterProfiler on the DESeq2 results to get WikiPathways enrichment
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(ggsci)
library(org.At.tair.db)
library(enrichplot)
library(DOSE)

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
n_category <- snakemake@params[["n_category"]]

# Define directories
plot_dir <- snakemake@output[["plots_dir"]]
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}

# For each result_name, load the corresponding res and res.shrink
result_names <- resultsNames(dds)[!resultsNames(dds) == "Intercept"]

# Function to extract the TAIR_hit values from the subset of transcripts that are DE
get_tair_hits <- function(df, res, dds) {
  
  print("Number of DE transcripts:")
  print(length(res))
  
  # Extract the subset of the dataframe where target_id is in the vector list
  subset_df <- df[df$target_id %in% res & df$TAIR_hit != "", ]
  print(subset_df)
  # Remove the ".X" pattern from the TAIR_hit column
  subset_df$TAIR_hit <- sub("\\.[0-9]$", "", subset_df$TAIR_hit)
  print(subset_df$TAIR_hit)
  # Extract the unique TAIR_hit values from this subset
  tair_hits_list <- unique(subset_df$TAIR_hit)
  print(tair_hits_list)
  return(tair_hits_list)
}

for (res_name in result_names) {
    res_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds")
    res <- readRDS(res_file_path)
    
    print(res)
    
    res_shrink_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds")
    res_shrink <- readRDS(res_shrink_file_path)
    
    print(res_shrink)

    # Create the directory if doesn't exist
    if (!dir.exists(snakemake@output[["tables_dir"]])) {
        dir.create(snakemake@output[["tables_dir"]], recursive = TRUE, showWarnings = FALSE)
    }

    res_shrink_signif <- res_shrink[res_shrink$padj < padj_threshold & !is.na(res_shrink$padj) & abs(res_shrink$log2FoldChange) >= lfc_threshold, ]

    res_shrink_signif_genes <- rownames(res_shrink_signif)

    # Get TAIR_hit values for significant DEGs
    signif_degs_tair_hits <- get_tair_hits(tcp_to_gene_df, res_shrink_signif_genes, dds)

    signif_degs_tair_hits_entrez = bitr(signif_degs_tair_hits, fromType="TAIR", toType="ENTREZID", OrgDb="org.At.tair.db")
    print(signif_degs_tair_hits_entrez$ENTREZID)

    # Check if gene ids in signif_degs_tair_hits_entrez$ENTREZID is in the org.At.tair.db
    # If not, remove them

    enrich_wikip <- enrichWP(signif_degs_tair_hits_entrez$ENTREZID, organism = "Arabidopsis thaliana")
    print(enrich_wikip)

    # if enrich_wikip NUL then skip the rest of the code
    if (is.null(enrich_wikip)) {
        next
    }

    # Plot the results
    mytheme <- theme(panel.border=element_blank(),
        panel.grid.major=element_line(linetype='dotted', colour='#808080'),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.x = element_line())
    
    enrich_wikip_plot <- ggplot(enrich_wikip, showCategory=10,
                        aes(GeneRatio, fct_reorder(Description, GeneRatio), fill=p.adjust)) +
                        geom_col() +
                        scale_x_continuous(expand=c(0,0)) +
                        scale_fill_gradientn(colours=c('#b3eebe', "#46bac2", '#371ea3'),
                                        guide=guide_colorbar(reverse=TRUE)) +
                        theme_dose(12) +
                        mytheme +
                        xlab("GeneRatio") + ylab(NULL) +
                        ggtitle("WikiPathways")
    
    # enrich_wikip_plot <- generate_GO_enrich_dotplot(go_results=enrich_wikip, 
    #                                                 n_category=n_category, 
    #                                                 font_size=14, 
    #                                                 label_format=45, 
    #                                                 title="")

    # Save the plot
    ggsave(paste0(plot_dir, "/", "wikip.bar_plot", analysis_name, ".", res_name, ".png"), enrich_wikip_plot, width=12, height=10, dpi=350)

    # Write the results to file
    write.csv(enrich_wikip,
            paste0(snakemake@output[["tables_dir"]], "/wikipathways_enrichment.", analysis_name, ".", res_name, ".csv"), 
            row.names = TRUE)

}

sink(type="message")
sink()
