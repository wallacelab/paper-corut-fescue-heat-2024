# Purpose: Generate exploratory plots for the DESeq2 analysis
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

# Check if necessary libraries are installed
required_packages <- c("DESeq2", "ggplot2", "pheatmap", "gridExtra", "ggpubr", "ggsci", "RColorBrewer", "reshape2", "PCAtools")
stopifnot(all(required_packages %in% installed.packages()))

# Load the libraries
lapply(required_packages, require, character.only = TRUE)

# Source the additional functions script using the path from Snakemake's params
funcs_script <- snakemake@params[['funcs_script']]
stopifnot(file.exists(funcs_script))
source(funcs_script)

# Load the dds and vsd objects
dds_file <- snakemake@input[['dds_subset']]
vsd_file <- snakemake@input[['vsd_subset']]
stopifnot(file.exists(dds_file), file.exists(vsd_file))

dds <- readRDS(dds_file)
vsd <- readRDS(vsd_file)

# Extract the analysis name from the input file for use in plot titles and saving files
analysis_name <- gsub("dds_(.*).rds", "\\1", basename(dds_file))

color_by_var <- snakemake@params[['color_by']]
print(color_by_var)
shape_by_var <- snakemake@params[['shape_by']]
print(shape_by_var)

print(colnames(colData(dds)))
# Check if color_by_var and shape_by_var are valid columns in the data
stopifnot(color_by_var %in% colnames(colData(dds)), shape_by_var %in% colnames(colData(dds)))

# ===================================================== #
# ================  Exploratory Plots  ================ #
# ===================================================== #

# Set the color palette
cbPalette <- c("#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00",  "#999999")

# Set the number of top genes to plot
ntop <- snakemake@params[['n_top']]

# if shape_by_var is "treatment", then set the shape legend title to "Endophyte status"
# if shape_by_var is "condition", then set the shape legend title to "Condition"
if (shape_by_var == "treatment") {
    shape_legend_title <- "Endophyte status"
} else if (shape_by_var == "condition") {
    shape_legend_title <- "Condition"
} else {
    shape_legend_title <- shape_by_var
}

# Plot PCA for PC1 and PC2
pca.pc12 <- plot_pca(vsd, ntop=ntop, pcs=c("PC1", "PC2"), 
                     color_by=color_by_var, shape_by=shape_by_var, 
                     legend_titles= list(shape=shape_legend_title, color="Time points"), 
                     color_order=c("HR0", "HR2", "HR12", "D1" ,"D3", "W1"), color_pal=cbPalette)
print("Saving PCA plot for PC1 and PC2")
ggsave(snakemake@output[['pca_pc1_pc2']], pca.pc12, width=14, height=10, units="in", dpi=350)

# Plot PCA for PC2 and PC3
pca.pc23 <- plot_pca(vsd, ntop=ntop, pcs=c("PC2", "PC3"), 
                     color_by=color_by_var, shape_by=shape_by_var, 
                     legend_titles= list(shape=shape_legend_title, color="Time points"),
                     color_order=c("HR0", "HR2", "HR12", "D1" ,"D3", "W1"), color_pal=cbPalette)

# Plot PCA for PC3 and PC4 
pca.pc34 <- plot_pca(vsd, ntop=ntop, pcs=c("PC3", "PC4"), 
                     color_by=color_by_var, shape_by=shape_by_var, 
                     legend_titles= list(shape=shape_legend_title, color="Time points"), 
                     color_order=c("HR0", "HR2", "HR12", "D1" ,"D3", "W1"), color_pal=cbPalette)

# Plot the combined PCA plot
pca.combined <- ggpubr::ggarrange(pca.pc12, pca.pc34,
                            ncol = 2, nrow = 1, align = "hv",
                            common.legend = TRUE, legend = "right")

print("Saving combined PCA plot")
ggsave(snakemake@output[['pca_combined']], pca.combined, width=28, height=14, dpi=350)

# Eigen Correlation Plot
print("Saving eigen correlation plot")
pcs <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)
png(snakemake@output[['eigencor']], width=14, height=10, units="in" ,res=350)
eigencorplot(pcs, metavars = c(shape_by_var, color_by_var, 'plate'))
dev.off()

# Box Plot - Cook Distances
print("Saving Cook's distance box plot")
png(snakemake@output[['cooks']], width=20, height=14, units="in", res=350)
plot_cooks_dist(dds, analysis_name)
dev.off()

# Plot Dispersion Estimates
print("Saving dispersion estimates plot")
png(snakemake@output[['disp_est']], width=10, height=6, units="in", res=350)
plotDispEsts(dds)
dev.off()

# Desnsity Plot
png(snakemake@output[['density']], width=10, height=8, units="in", res=350)
plot_density(dds)
dev.off()

# Heatmap - Sample Distances
print("Calculating sample distances")
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- sub(".*?_", "", rownames(t(assay(vsd))))
colnames(sample_dist_matrix) <- NULL
print("Sample distance matrix:")
print(head(sample_dist_matrix))
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

sample_dists_heat <- pheatmap(sample_dist_matrix, 
         clustering_distance_rows=as.dist(sample_dists), 
         clustering_distance_cols=as.dist(sample_dists), 
         col=colors, fontsize_row=10)

print("Saving sample distance heatmap")
png(snakemake@output[['dist_heatmap']], width=12, height=8, units="in", res=350)
sample_dists_heat
dev.off()

# Heatmap - Count Matrix
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
annotation_data <- as.data.frame(colData(vsd)[, c(shape_by_var, color_by_var)])
count_matrix_heat <- pheatmap(assay(vsd)[select,], 
         cluster_rows=FALSE, show_rownames=FALSE, 
         cluster_cols=FALSE, annotation_col=annotation_data)
print("Saving count matrix heatmap")
png(snakemake@output[['count_heatmap']], width=12, height=8, units = "in", res=350)
count_matrix_heat
dev.off()

# Heatmap - Gene Clusters
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 20)
mat <- assay(vsd)[topVarGenes,]
mat <- mat - rowMeans(mat)
annotation_data <- as.data.frame(colData(vsd)[, c(shape_by_var, color_by_var)])
gene_cluster_heatmap <- pheatmap(mat, annotation_col=annotation_data)
print("Saving gene cluster heatmap")
png(snakemake@output[['gene_cluster_heatmap']], width=12, height=8, units = "in", res=350)
gene_cluster_heatmap
dev.off()

# Close logging
sink(type="message")
sink()