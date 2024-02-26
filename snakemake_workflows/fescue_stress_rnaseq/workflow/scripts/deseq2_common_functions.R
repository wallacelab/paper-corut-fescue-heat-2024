# Date: 2023-11-02
# Author: Adnan Kivanc Corut

# This script contains common functions used in the DESeq2 analysis

# ======================================================================================== #
#        Functions for Differential Expression Analysis of FescueHeatStress Project        #
# ======================================================================================== #

# ====== Differential Expression Analysis Functions ====== #

# Plot Cooks Distance function
#"""
#This function generates a boxplot of Cooks distances
#Parameters:
#    object: DESeq2 object
#    title: title of the plot
#"""
plot_cooks_dist <- function(object, title, outlier_highlight=FALSE) {
  par(mar = c(12, 5, 2, 2), las = 2, cex = 1.5, cex.axis = 0.5)
  
  # Log transformation of the data
  log_data <- log10(assays(object)[["cooks"]])
  
  # Add a main title and labels
  main_title <- title
  xlab <- ""
  ylab <- "Log10(cooks distances)"
  
  # Create boxplot with color
  boxplot(log_data, range = 0, las=2, notch = TRUE,
          main = main_title, xlab = xlab, ylab = ylab,
          col = "lightblue", outline = TRUE)
  
  if (outlier_highlight == TRUE){
    # Adding outliers in different color
    outlier_values <- boxplot.stats(log_data)$out  # getting outlier values
    points(which(log_data %in% outlier_values), outlier_values, col = "red", pch = 19, cex = 2)
  }
}

# Plot PCA function
#"""
#This function generates a PCA plot
#Parameters:
#    object: variance stabilized or rlogged DESeq2 object
#    ntop: number of genes to consider
#   pcs: principal components to plot
#   color_by: color by
#   shape_by: shape by
#   color_pal: color palette
#"""
plot_pca <- function(object, ntop=500, pcs=c("PC1", "PC2"), color_by="condition", shape_by="condition", 
                      legend_titles = NULL, color_order = NULL, shape_order = NULL, color_pal){

  ntop <- 500
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t(assay(object)[select, ])
  pca <- prcomp(mat)
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  
  # Add sample information to PCA
  plot_pca_data <- data.frame(pca$x, colData(object))
  
  # Set x and y limits
  xlim <- c(min(plot_pca_data[, pcs[[1]]]), max(plot_pca_data[, pcs[[1]]]))
  ylim <- c(min(plot_pca_data[, pcs[[2]]]), max(plot_pca_data[, pcs[[2]]]))

  if (is.null(legend_titles) || is.null(legend_titles$color)) {
    legend_titles$color <- color_by
  }
  if (is.null(legend_titles) || is.null(legend_titles$shape)) {
    legend_titles$shape <- shape_by
  }

  # Assuming plot_pca_data is the data frame you're plotting
  if (!is.null(color_order)) {
    plot_pca_data[[color_by]] <- factor(plot_pca_data[[color_by]], levels = color_order)
  }
  if (!is.null(shape_order)) {
    plot_pca_data[[shape_by]] <- factor(plot_pca_data[[shape_by]], levels = shape_order)
  }
  
  # Main PCA Plot
  pMain <- ggplot(plot_pca_data, aes(!!sym(pcs[[1]]), !!sym(pcs[[2]]), color = !!sym(color_by), shape = !!sym(shape_by))) +
    geom_point(size = 6) +
    scale_colour_manual(values=color_pal) +
    theme_bw() +
    labs(colour = legend_titles$color, shape = legend_titles$shape) +
    scale_x_continuous(limits = xlim) + scale_y_continuous(limits = ylim) +
    theme(legend.position = 'right', legend.box = 'horizontal',
          legend.direction = 'vertical',
          legend.key.height = unit(0.2, 'cm'),
          legend.key.width = unit(0.1, 'cm'),
          legend.title = element_text(size = rel(1.6)),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.text = element_text(size = rel(1.4)),
          axis.text.y = element_text(size = 22),
          axis.text.x = element_text(size = 22),
          axis.title.x = element_text(size = 24), 
		      axis.title.y = element_text(size = 24),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  if (pcs[[1]] == "PC1" & pcs[[2]] == "PC2") {
    pMain <- pMain +
      xlab(paste0(pcs[[1]], ": ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0(pcs[[2]], ": ",round(percentVar[2] * 100),"% variance"))
  } else if (pcs[[1]] == "PC2" & pcs[[2]] == "PC3") {
    pMain <- pMain +
      xlab(paste0(pcs[[1]], ": ",round(percentVar[2] * 100),"% variance")) +
      ylab(paste0(pcs[[2]], ": ",round(percentVar[3] * 100),"% variance"))
  } else if (pcs[[1]] == "PC3" & pcs[[2]] == "PC4") {
    pMain <- pMain +
      xlab(paste0(pcs[[1]], ": ",round(percentVar[3] * 100),"% variance")) +
      ylab(paste0(pcs[[2]], ": ",round(percentVar[4] * 100),"% variance"))
  } else {
    pMain <- pMain +
      xlab("") +
      ylab("") 
  }

  xlim.update <- layer_scales(pMain)$x$get_limits()
  ylim.update <- layer_scales(pMain)$y$get_limits()

  # Top Density Plot
  pTop <- ggplot(plot_pca_data, aes(!!sym(pcs[[1]]), fill = !!sym(color_by), linetype = !!sym(shape_by))) +
    geom_density(linewidth = 0.2, alpha = 0.5) + ylab('Density') +
    theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = rel(1.5)),
        axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(), legend.position = 'none') +
    guides(fill = "none", color = "none", linetype = "none") + 
    scale_fill_manual(values=color_pal) +
    scale_x_continuous(limits = xlim.update)

  # Right Density Plot
  pRight <- ggplot(plot_pca_data, aes(!!sym(pcs[[2]]), fill = !!sym(color_by), linetype = !!sym(shape_by))) +
    geom_density(linewidth = 0.2, alpha = 0.5) +  coord_flip() + ylab('Density') +
    theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank(), legend.position = 'none') +
    guides(fill = "none", color = "none", linetype = "none") +
    scale_fill_manual(values=color_pal) +
    scale_x_continuous(limits = ylim.update)

  legend <- get_legend(pMain)
  
  # Create and return the arranged ggplot object
  combined_plot <- ggarrange(pTop, legend, pMain + theme(legend.position = 'none'), pRight,
             ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  
  return(combined_plot)   
}

# Helper function to extract the legend
get_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# Plot PCAs function (plots PC1vs.PC2, PC2vs.PC3, and PC3vs.PC4)
#"""
#This function generates a PCA plot
#Parameters:
#    sample_var: sample variable to color by
#    object: variance stabilized or rlogged DESeq2 object
#"""
plot_PCAs <- function(sample_var, object) {
  
    ntop <- 500 # number of genes to consider
    rv <- rowVars(assay(object)) # get row-wise variance
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))] # select top n genes
    mat <- t( assay(object)[select, ] ) # transpose matrix
    pca <- prcomp(mat) # perform PCA
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    # Add sample information to PCA
    plot_pca_data <- data.frame(pca$x, colData(object))
    ## Define sample colors depending on the sample variable
    
    if (sample_var == "genotype") {
        colors <- pal_npg()(4)
    } else if (sample_var == "condition") {
        colors <- pal_igv()(2)
    } else if (sample_var == "treatment") {
        colors <- pal_simpsons()(2)
    } else if (sample_var == "plate") {
        colors <- pal_jama()(4)
    } else if (sample_var == "time_point") {
        colors <- pal_lancet()(6)
    }

    plot1 <- ggplot(as.data.frame(plot_pca_data), aes(PC1, PC2, shape = !!rlang::sym(sample_var))) +
        geom_point(aes_string(color=sample_var), size = 4) +
        theme_bw() +
        theme(plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        scale_color_manual(values = colors)
    
    plot1 <- plot1 +
        xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
        ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))

    plot2 <- ggplot(as.data.frame(plot_pca_data), aes(PC2, PC3, shape = !!rlang::sym(sample_var))) +
        geom_point(aes_string(color=sample_var), size = 4) +
        theme_bw() +
        theme(plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        scale_color_manual(values = colors)
    
    plot2 <- plot2 +
        xlab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
        ylab(paste0("PC3: ",round(percentVar[3] * 100),"% variance"))

    plot3 <- ggplot(as.data.frame(plot_pca_data), aes(PC3, PC4, shape = !!rlang::sym(sample_var))) +
        geom_point(aes_string(color=sample_var), size = 4) +
        theme_bw() +
        theme(plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        scale_color_manual(values = colors)

    plot3 <- plot3 +
        xlab(paste0("PC3: ",round(percentVar[3] * 100),"% variance")) +
        ylab(paste0("PC4: ",round(percentVar[4] * 100),"% variance"))

    wplot <- ggarrange(plot1, plot2, plot3, ncol = 3, common.legend = TRUE, legend = "right")

    return(wplot)
}

# Plot density function
#"""
#This function generates a density plot of the log2 transformed counts
#Parameters:
#    object: DESeq2 object
#"""
plot_density <- function(object) {
  count.melt <- reshape2::melt(log2(counts(object) + 1))
  
  ggplot(data=count.melt, mapping=aes(x=value, color=Var2)) + 
    ylab("Density") +
    xlab("log2(counts + 1)") +
    theme_bw() +
    geom_density() + 
    theme(legend.position="none",
          axis.title = element_text(size = rel(2)),
          axis.text = element_text(size = rel(1.5)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

# Plot top N log2 fold change function
#"""
#This function generates a barplot of the top N log2 fold change
#Parameters:
#    res: DESeq2 results
#    top_n: number of genes to plot
#"""
plot_top_lfc_bar <- function(res, top_n=20) {
  
  # Make sure ggsci is installed
  if (!requireNamespace("ggsci", quietly = TRUE)) {
    install.packages("ggsci")
  }

  # Convert DESeqResults to a standard dataframe
  res <- as.data.frame(res)
  
  # Extract top 20 upregulated and 20 downregulated genes
  top_genes <- res %>%
    arrange(log2FoldChange) %>%
    head(top_n)
  
  bottom_genes <- res %>%
    arrange(-log2FoldChange) %>%
    head(top_n)
  
  # Combine and sort the genes
  combined_genes <- rbind(top_genes, bottom_genes)
  
  # color palette
  my_colors <- pal_d3()(2)
  # change the order of the colors
  my_colors <- c(my_colors[2], my_colors[1])
  
  # Plot
  p <- ggplot(combined_genes, aes(x=reorder(row.names(combined_genes), log2FoldChange), y=log2FoldChange)) +
    geom_bar(stat='identity', aes(fill=(log2FoldChange > 0)), width=0.7) +
    scale_fill_manual(values=my_colors, 
                  labels=c("Down-regulated", "Up-regulated"),
                  name="") +
    coord_flip() +
    geom_hline(yintercept=0, linetype="solid", color="black") +
    labs(title="", x="Transcripts", y=bquote(~Log[2]~ "fold change")) +
    theme_classic() +
    # don't show y axis line and ticks
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.title.x = element_text(size = 18), 
		      axis.title.y = element_text(size = 18),
          legend.title=element_blank(),
          legend.text=element_text(size=16),
          legend.position="top")
  
  return(p)

}

# Plot number of DEGs (upregulated and downregulated) function'
#"""
#This function generates a barplot of the number of DEGs
#Parameters:
#    res_list: list of DESeq2 results
#    padj_cutoff: adjusted p-value cutoff
#    lfc_cutoff: log2 fold change cutoff
#"""
plot_deg_bar <- function(res_list, padj_cutoff = 0.05, lfc_cutoff = 0.58) {

    # Make sure ggsci is installed
  if (!requireNamespace("ggsci", quietly = TRUE)) {
    install.packages("ggsci")
  }
  
  my_colors <- pal_d3()(2)

  # Calculate the number of upregulated and downregulated genes for each table
  up_counts <- vapply(res_list, function(df) {
    sum(!is.na(df$padj) & !is.na(df$log2FoldChange) & df$padj <= padj_cutoff & df$log2FoldChange >= lfc_cutoff)
  }, numeric(1))
  
  down_counts <- vapply(res_list, function(df) {
    sum(!is.na(df$padj) & !is.na(df$log2FoldChange) & df$padj <= padj_cutoff & df$log2FoldChange <= -lfc_cutoff)
  }, numeric(1))
  
  # Setting the order of the 'Table' column based on the order in the list
  data <- data.frame(
    Table = factor(rep(names(up_counts), 2), levels = names(res_list)),
    Count = c(up_counts, -down_counts),
    Type = factor(c(rep("Up-regulated", length(up_counts)), rep("Down-regulated", length(down_counts))),
                  levels = c("Up-regulated", "Down-regulated"))
  )
  
  # Plot using ggplot2
  p <- ggplot(data, aes(x = Table, y = Count, fill = Type)) +
    geom_bar(stat = "identity", width=0.6) +
    geom_hline(yintercept = 0, color = "darkgrey", linetype = "longdash") +
    geom_text(aes(label = abs(Count), fontface = "bold"), position = position_nudge(y = ifelse(data$Count > 0, -max(up_counts)*0.1, max(down_counts)*0.1)), size = 6) +
    scale_fill_manual(values=my_colors, 
                  labels=c("Up-regulated", "Down-regulated")) +
    labs(caption = "*padj <= 0.05, |log2FoldChange| >= 0.58", y = "Number of DEGs", x = "Comparisons", title = NULL) +  # Removed title
    theme_bw() +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.subtitle = element_text(size = 16),
      plot.caption = element_blank(),
      # axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_blank(),
      legend.title=element_blank(),
      legend.text = element_text(size = 20),
      legend.position = "top",  # Moved legend to top
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Volcano plot function
#"""
#This function generates a volcano plot from DESeq2 results
#Parameters:
#    deseq2.res: DESeq2 results
#    title: title of the plot
#    subtitle: subtitle of the plot
#    y_upper_lim: upper limit of y-axis
#"""
generate_volcano_plot <- function(deseq2.res, p_cutoff=10e-3, lfc_cutoff=0.58, title="", subtitle="", caption=""){

  # get y-axis limits based on min p-value
  y_upper_lim <- min(deseq2.res$pvalue)
  
  # Find max and min log2 fold change values
  max_log2FC <- max(deseq2.res$log2FoldChange)
  min_log2FC <- min(deseq2.res$log2FoldChange)

  # Define the color palette
  mypal = pal_nejm()(5)
  
  # Generate the volcano plot
  volcano_plot <- EnhancedVolcano(deseq2.res,
                  lab = rownames(deseq2.res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = title,
                  subtitle = subtitle,
                  selectLab = c(""),
                  caption = bquote(~Log[2]~ "fold change cutoff, " ~ .(lfc_cutoff) ~ "; p-value cutoff, " ~ .(p_cutoff)),
                  pCutoff = p_cutoff,
                  FCcutoff = lfc_cutoff,
                  pointSize = 4,
                  xlim = c(min_log2FC - 0.5, max_log2FC + 0.5),
                  ylim = c(0, -log10(y_upper_lim)),
                  col=mypal,
                  colAlpha = 1,
                  legendPosition = "bottom",
                  #drawConnectors=TRUE,
                  #widthConnectors = 0.6,
                  #boxedLabels = FALSE,
                  #colConnectors = 'black',
                  labFace = 'bold',
                  labSize = 6,
                  #shape = c(4, 6, 18, 16),
                  #labCol = 'black',
                  legendLabSize = 16,
                  legendIconSize = 8,
                  border = 'full',
				          gridlines.major = FALSE,
  				        gridlines.minor = FALSE)

	return(volcano_plot)
}

plot_exp_heatmap <- function(dds, vsd, res, top_n, time_point, alternate_color=FALSE) {
  
  # Subset the dds object by tissue type, if given
  if (!is.null(time_point)) {
    dds <- dds[,dds$time_point == time_point]
    vsd <- vsd[,vsd$time_point == time_point]
  }

  # Get the significant genes
  sig_genes <- subset(res, padj < 0.05) 
  # Sort by absolute log2FoldChange
  sig_genes <- sig_genes[order(sig_genes$padj, decreasing = FALSE),]  
  top_genes <- head(sig_genes, top_n) # get the top n genes

  # Extract the gene names of the top genes
  gene_names <- rownames(top_genes)

  # Subset the VSD matrix to only include the top genes
  vsd_top_genes <- vsd[gene_names,]

  # Scale the rows of the gene counts matrix to have a mean of 0 and a standard deviation of 1
  vsd_counts_scaled <- t(scale(t(assay(vsd_top_genes))))

  # Center the data
  vsd_counts_scaled <- (vsd_counts_scaled - rowMeans(vsd_counts_scaled)) / rowSds(vsd_counts_scaled)

  # Get the metadata information for Time_point, Genotype, and Location
  ha = HeatmapAnnotation(df = colData(dds)[,c("treatment", "condition")],
                          col = list(treatment = c("infected" = "orange", "uninfected" = "purple"),
                                     condition= c("NC" = "deepskyblue3", "HS" = "brown2")),
                         annotation_name_side = "right",
                         annotation_legend_param = list(treatment = list(title = 'Treatment',
                                                                          legend_direction = 'vertical',
                                                                          title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                                                          labels_gp = gpar(fontsize = 14)),
                                                        condition = list(title = 'Condition',
                                                                          legend_direction = 'vertical'),
                                                                          title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                                                          labels_gp = gpar(fontsize = 14)),
                                                                          annotation_name_gp= gpar(fontsize = 16))

  ## FC's of DEGs as row annotation
  FCs <- data.frame(signif(2**(top_genes$log2FoldChange), digits = 3))
  rownames(FCs)<-rownames(top_genes)
  FCs<-data.frame("FC"=apply(FCs, 1, function(x){if(x>1|x==1){paste(">=1")} else {paste("<1")}}))

  left_ans <- rowAnnotation(
     FC = FCs$FC,
     col = list("FC" = c(">=1" = "darkgoldenrod4", "<1" = "darkgoldenrod1")),
     annotation_name_gp= gpar(fontsize = 16),
     annotation_legend_param = list(
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 16, 
                    fontface = "bold"), 
      labels_gp = gpar(fontsize = 14)))

  # If altertate_color is false, use the following color scheme if TRUE, use the other color scheme 
  if (alternate_color == FALSE) {
    myPalette <- colorRampPalette(rev(c("#D73027", "#FC8D59", 
      "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100)
  } else {
    myPalette <- colorRampPalette(c('forestgreen', 'black', 'purple'))(100)
  }
  
  myBreaks <- seq(-3, 3, length.out = 100)
  
  heatmap_col = colorRamp2(myBreaks, myPalette)

  # Remove the plate name from the sample names 
  sample_names <- colnames(vsd_counts_scaled)
  sample_names <- gsub("-plate[0-9]", "", sample_names)
  colnames(vsd_counts_scaled) <- sample_names

  htmap = Heatmap(vsd_counts_scaled, 
          col = colorRamp2(c(-4, -0.0001, 00001, 4), c("darkblue", "lightblue", "lightsalmon", "darkred")), 
          # cluster_rows = TRUE, cluster_columns = FALSE,
          name = "Z-Score",
          top_annotation = ha,
          left_annotation = left_ans,
          row_km = 2,
          column_km = 2,
          row_title = NULL,
          column_title = NULL,
          column_names_gp = gpar(fontsize = 14),
          row_names_gp = grid::gpar(fontsize = 14),
          heatmap_legend_param = list(
                            legend_direction = 'vertical',
                            title_gp = gpar(col = "black", fontsize = 16, fontface="bold"),
                            labels_gp = gpar(fontsize= 14)))
  draw(htmap, merge_legend = TRUE, column_dend_side = "top",
        padding = unit(c(10, 10, 2, 16), "mm")) #bottom, left, top, right paddings
}