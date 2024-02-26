# Purpose: Plot PCA for DESeq2 analysis
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("tidyverse")
library("DESeq2")
library("genefilter")
library("ggpubr")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
vsd <- vst(dds, blind = FALSE)

# Set color palette
cbPalette <- c("#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00",  "#999999")

# Plot PCA - PC1 vs PC2
pca.plot <- plotPCA(vsd, intgroup = snakemake@wildcards[["variable"]]) +
    geom_point(aes(color = !!sym(snakemake@wildcards[["variable"]])),
				size = 4) +
    theme_bw() + 
  	theme(plot.title = element_text(size = 16, face = "bold"),
        	plot.caption = element_text(size = 12, face = "italic"),
		  	axis.text.y = element_text(size = 14),
		  	axis.text.x = element_text(size = 14),
			axis.title.x = element_text(size = 16),
			axis.title.y = element_text(size = 16),
			legend.title=element_text(size=16), 
    		legend.text=element_text(size=14),
			panel.grid.major = element_blank(),
    		panel.grid.minor = element_blank()) +
	labs(color = snakemake@wildcards[["variable"]]) +
	scale_colour_manual(values=cbPalette)

ggsave(snakemake@output[[1]], plot=pca.plot, width=16, height=10, dpi=350)

# Plot PCA - PC3 vs PC4

ntop <- 500 # number of genes to consider
rv <- rowVars(assay(vsd)) # get row-wise variance
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))] # select top n genes
mat <- t( assay(vsd)[select, ] ) # transpose matrix
pca <- prcomp(mat) # perform PCA

# Calculate percent variance
percentVar <- pca$sdev^2 / sum(pca$sdev^2)

# Add sample information to PCA
plot_pca <- data.frame(pca$x, colData(dds))

# Get the column name to convert
shape_var <- snakemake@config[["params"]][["deseq2"]][["pca"]][["shape"]]

# Check if shape_var is not null
if (!is.null(shape_var)) {
  	# Check if shape_var is in the column names of plot_pca
  	if (shape_var %in% colnames(plot_pca)) {
  	  # Convert the column to a factor
  	  	plot_pca[[shape_var]] <- as.factor(plot_pca[[shape_var]])
  	} else {
  	  	stop(paste("The column", shape_var, "does not exist in plot_pca dataframe"))
  	}
}

if(is.factor(plot_pca[[shape_var]])) {
  print(paste(shape_var, "is a factor"))
} else {
  print(paste(shape_var, "is not a factor"))
}

# Change the levels of the factor if shape_ver == "time_point"
if (shape_var == "time_point") {
  	plot_pca[[shape_var]] <- factor(plot_pca[[shape_var]], 
									levels = c("HR0", "HR2", "HR12", "D1" ,"D3", "W1"))
}

if (snakemake@wildcards[["variable"]] == "treatment") {
  	# Set the levels of the factor
  	cbPalette <- c("#F0E442", "#0072B2")
}

# Calculate xlim	
xlim_pc1 <- c(min(plot_pca$PC1, na.rm = TRUE)-0.5, max(plot_pca$PC1, na.rm = TRUE)+0.5)

xlim_pc3 <- c(min(plot_pca$PC3, na.rm = TRUE)-0.5, max(plot_pca$PC3, na.rm = TRUE)+0.5)

# Plot PCA with pc1 and pc2
pca.plot.pc1_2 <- ggplot(plot_pca, aes(x = PC1, y = PC2)) +
  	geom_point(aes(color = !!sym(snakemake@wildcards[["variable"]]),
				shape = !!sym(snakemake@config[["params"]][["deseq2"]][["pca"]][["shape"]])), 
				size = 4.5) +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
	# add xlim and ylim
	scale_x_continuous(limits = xlim_pc1) +		
  	# scale_x_continuous(expand = c(0.3,  0.3)) +
	theme_bw() + 
  	theme(plot.title = element_text(size = rel(1.5)),
		  	axis.text.y = element_text(size = 20),
		  	axis.text.x = element_text(size = 20),
			axis.title.x = element_text(size = 22),
			axis.title.y = element_text(size = 22),
  	    legend.title = element_text(size = 20),
  	    legend.text = element_text(size = 18),
		panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank()) +
	labs(color = snakemake@wildcards[["variable"]],
		shape= snakemake@config[["params"]][["deseq2"]][["pca"]][["shape"]]) +
  	scale_colour_manual(values=cbPalette)

# Plot PCA with pc3 and pc4
pca.plot.pc3_4 <- ggplot(plot_pca, aes(x = PC3, y = PC4)) +
  	geom_point(aes(color = !!sym(snakemake@wildcards[["variable"]]),
				shape = !!sym(snakemake@config[["params"]][["deseq2"]][["pca"]][["shape"]])), 
				size = 4.5) +
      xlab(paste0("PC3: ",round(percentVar[3] * 100),"% variance")) +
      ylab(paste0("PC4: ",round(percentVar[4] * 100),"% variance")) +
	scale_x_continuous(limits = xlim_pc3) +
  	# scale_x_continuous(expand = c(0.3,  0.3)) +
    theme_bw() + 
  	theme(plot.title = element_text(size = rel(1.5)),
		  	axis.text.y = element_text(size = 20),
		  	axis.text.x = element_text(size = 20),
			axis.title.x = element_text(size = 22),
			axis.title.y = element_text(size = 22),
  	    legend.title = element_text(size = 20),
  	    legend.text = element_text(size = 18),
		panel.grid.major = element_blank(),
  	  	panel.grid.minor = element_blank()) +
	labs(color = snakemake@wildcards[["variable"]],
		shape= snakemake@config[["params"]][["deseq2"]][["pca"]][["shape"]]) +
  	scale_colour_manual(values=cbPalette)

# Combine plots
combined.plot <- ggarrange(pca.plot.pc1_2, pca.plot.pc3_4,
                            ncol = 2, nrow = 1, 
                            common.legend = TRUE, legend = "right")

# Save the plot
ggsave(snakemake@output[[2]], plot=combined.plot, width=16, height=8, dpi=600)