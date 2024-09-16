# Exploratory data analysis visualization
# RNASEq Utils

save_it <- function(image_object, filepath, filename, resolution=300, w=800, h=650, format=NULL){
  if (tolower(format) == "png"){
    png(paste0(filepath,"/",filename,".png"), res=resolution, width=w, height=h)
    print(image_object)
    dev.off()
  }
  if (tolower(format) == "pdf"){
    pdf(paste0(filepath,"/",filename,".pdf"), width=w/96, height=h/96)
    print(image_object)
    dev.off()
  }
  }

total_counts_ggplot <- function(counts_data, groupby=NULL, design, type = "", fontsize=18) {
  counts <- counts_data
  memo <- ""
  
  if (ncol(counts) < 31) {
    x_axis_labels <- fontsize
  } else {
    x_axis_labels <- fontsize-4
  }
  
  plot_data <- data.frame(
    sample = as.factor(colnames(counts)),
    counts = colSums(counts) / 1e6,
    if (is.character(groupby)){
      grouping = lapply(design[groupby], as.character)
    }else{
      grouping = NULL
    }
  )
  print(plot_data)
  plot <- ggplot2::ggplot(
    data = plot_data,
    ggplot2::aes(x = sample, y = counts, fill = as.character(flatten(design["condition"])))
  )
  
  plot <- plot +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = fontsize), 
      legend.key.size = ggplot2::unit(2, 'cm'),
      legend.title = ggplot2::element_text(size = 0),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = fontsize
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      axis.text.y = ggplot2::element_text(
        size = fontsize
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = fontsize,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = paste("Total", type, "Read Counts (Millions)", memo),
      y = paste(type, "Counts (Millions)")
    )
  
  return(plot)
}
eda_boxplot <- function(processed_data, design, fontsize=18) {
  counts <- as.data.frame(processed_data)
  memo <- ""
  
  if (ncol(counts) < 31) {
    x_axis_labels <- fontsize
  } else {
    x_axis_labels <- fontsize-4
  }
  
  longer_data <- tidyr::pivot_longer(
    data = counts,
    colnames(counts),
    names_to = "sample",
    values_to = "expression"
  )
  
  longer_data$grouping <- rep(design$condition, nrow(counts))
  
  plot <- ggplot2::ggplot(
    data = longer_data,
    ggplot2::aes(x = sample, y = expression, fill = grouping)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = fontsize), 
      legend.key.size = ggplot2::unit(2, 'cm'),
      legend.title = ggplot2::element_text(size = 0),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = fontsize
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      axis.text.y = ggplot2::element_text(
        size = fontsize
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = fontsize,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = paste("Distribution of Transformed Data", memo),
      y = "Transformed Expression"
    )
  
  return(plot)
}


plot_count_distribution <- function(matrix, design, normalized=TRUE) {
  # make a colour vector
  statusCol <- as.numeric(factor(design$condition)) + 1
  # Check distributions of samples using boxplots
  plot <- boxplot(matrix, 
          xlab="", 
          ylab="",
          las=2,
          col=statusCol)
  # Let's add a blue horizontal line that corresponds to the median logCPM
  if(normalized){
    abline(h=median(as.matrix(matrix)), col="blue")
  }
  return(plot)
}

density_SD <- function(sds, n_genes_max, fontsize=18){
  #' Density plot of data standard deviation
  #'
  #' Draw a density plot of the standard deviation in the
  #' data. The function also adds a  vertical red lines for a specified range of
  #' genes.
  #'
  #' @param data Data matrix that has been through pre-processing
  #' @param n_genes_max Integer for the upper limit of gene range
  #' for visualization.
  #' @export
  #' @return Formatted density plot of the standard deviation
  #' distribution.
  
  
  sds <- apply(sds[, 1:dim(sds)[2]], 1, sd)
  # just for making the plot nicer, if not needed comment it out
  max_sd <- mean(sds) + 4 * sd(sds) 
  sds[sds > max_sd] <- max_sd
  sds <- as.data.frame(sds)
  xintercept_cut <- sort(sds$sds, decreasing = TRUE)[n_genes_max]
  plot <- ggplot2::ggplot(sds, ggplot2::aes(x = sds)) +
    ggplot2::geom_density(color = "darkblue", fill = "lightblue") +
    ggplot2::labs(
      title = "Standard Deviations of All Genes",
      y = "Density",
      x = "Standard Deviation"
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = fontsize), 
      legend.key.size = ggplot2::unit(2, 'cm'),
      legend.title = ggplot2::element_text(size = 0),
      plot.title = ggplot2::element_text(
        color = "black",
        size = fontsize,
        face = "bold",
        hjust = .5
      ),
      axis.text.x = ggplot2::element_text(
        size = fontsize
      ),
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = fontsize
      ),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = fontsize
      )) +
      ggplot2::geom_vline(
      ggplot2::aes(xintercept = xintercept_cut),
      color = "red",
      linetype = "dashed",
      linewidth = 1)
  return(plot)
}

# Heatmap of most variable genes
variance_heatmap <- function(data, n_genes_max){
  #' Heatmap of the highly variable genes 
  #' Exploratory data analysis (EDA)
  #'
  #' @param data transformed Data matrix that has been through pre-processing
  #' @param n_genes_max Integer for the upper limit of gene range
  #' for visualization.
  #' @export
  #' @return Formatted density plot of the standard deviation
  #' distribution.
  
  library(gplots)
  library(RColorBrewer)
  # We estimate the variance for each row in the logcounts matrix
  countVar <- apply(data, 1, var)
  # Get the row numbers for the top 500 most variable genes
  highVar <- order(countVar, decreasing=TRUE)[1:n_genes_max]
  # Subset logcounts matrix
  hmDat <- data[highVar,]
  
  # Get some nicer colours
  mypalette <- brewer.pal(11, "RdYlBu")
  morecols <- colorRampPalette(mypalette)
  # Set up colour vector for celltype variable
  col.cell <- c("purple","orange")[design$condition]
  
  # Plot the heatmap
  plot <- heatmap.2(hmDat, 
                    col=rev(morecols(50)),
                    trace="column", 
                    main="Top 500 most variable genes across samples",
                    ColSideColors=col.cell,scale="row")
  
  return(plot)
  
}

get_pc_variance <- function(data) {
  # subset data if more than 100 columns
  if (ncol(data) > 100) {
    part <- 1:100
    data <- data[, part]
  }
  # pca
  pca.object <- prcomp(t(data))
  
  # var proportions vector
  prop_var <- summary(pca.object)$importance[2, ] * 100
  
  return(prop_var |> round(1))
}


plot_PCA <- function(transformed_data, 
                     design,
                     PCAx = 1,
                     PCAy = 2,
                     selected_color = "Names",
                     selected_shape = "Names",
                     fontsize=18,
                     title='PCA',
                     point_size=10, label=design$caseid){
  
  
  if (ncol(transformed_data) < 31) {
    x_axis_labels <- fontsize
  } else {
    x_axis_labels <- fontsize-4
  }
  
  pca.object <- prcomp(t(transformed_data))
  
  # 5 pc's or number of columns if <5
  npc <- min(5, ncol(transformed_data))
  pcaData <- as.data.frame(pca.object$x[, 1:npc])
  print(head(pcaData))
  pcaData$Names <- design$condition
  
  point_size <- point_size

  plot <- ggplot2::ggplot(
    data = pcaData,
    aes(
      x = PC1,
      y = PC2,
      color = selected_color,
      shape = selected_shape,
      group = selected_color)
    #label=selected_shape)
  ) +
    # Preferred shapes
    ggplot2::scale_shape_manual(
      values = c(
        15, 16, 17, 18, 8, 9, 3, 4, 7, 10, 11, 12, 13, 14,
        0, 1, 2, 5, 6, 19, 20, 30:100
      )
    ) +
    ggplot2::geom_point(size = point_size, stroke=3) + 
    geom_text(aes(label = label, color="black")) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right", # TODO no legend for large data
      legend.text = ggplot2::element_text(size = fontsize), 
      legend.key.size = ggplot2::unit(2, 'cm'),
      legend.title = ggplot2::element_text(size = 0),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = fontsize
      ),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = fontsize
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      axis.text.y = ggplot2::element_text(
        size = fontsize
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = fontsize,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = title,
      y = "Dimension 2",
      x = "Dimension 1"
    ) #+
  
  # selected principal components
  PCAxy <- c(as.integer(PCAx), as.integer(PCAy))
  percentVar <- get_pc_variance(transformed_data)[PCAxy] # round(100 * summary(pca.object)$importance[2, PCAxy], 0)
  plot <- plot + ggplot2::xlab(paste0("PC", PCAx, ": ", percentVar[1], "% Variance"))
  plot <- plot + ggplot2::ylab(paste0("PC", PCAy, ": ", percentVar[2], "% Variance"))
  return(plot)
}

################# DESeq2 DIFFERENTIAL GENE EXPRESSION ANALYSIS #################

plot_diff_count <- function(diff_genes_list_table, plot_coding=T, fontsize=16){
  
  # Input ind diff_list_genes dataframe
  # Output is a cout bar plot
  
  exp_col = c("red","blue")
  up <-  diff_genes_list_table[diff_genes_list_table$log2FoldChange > 0,]
  down <- diff_genes_list_table[diff_genes_list_table$log2FoldChange < 0,]
  
  up_coding <- diff_genes_list[(diff_genes_list$log2FoldChange > 0) & 
                                 (diff_genes_list$TYPE == "protein_coding"), ] %>% dim()
  up_non_coding <- diff_genes_list[(diff_genes_list$log2FoldChange > 0) & 
                                     (!diff_genes_list$TYPE == "protein_coding"), ] %>% dim()
  down_coding <- diff_genes_list[(diff_genes_list$log2FoldChange < 0) & 
                                   (diff_genes_list$TYPE == "protein_coding"), ] %>% dim()
  down_non_coding <- diff_genes_list[(diff_genes_list$log2FoldChange < 0) & 
                                       (!diff_genes_list$TYPE == "protein_coding"), ] %>% dim()
  
  
  Regulation <- c("up","down")
  Counts = c(as.numeric(up_coding[1]), 
             as.numeric(up_non_coding[1]), 
             as.numeric(down_coding[1]), 
             as.numeric(down_non_coding[1]))
  
  df <- as.data.frame(Regulation = c("up","down"),
                      coding=c(up_coding[1], down_coding[1]),
                      non_coding=c(up_non_coding[1], down_non_coding[1]))
  print(df)

  p<-ggplot(df, aes(fill=Trans, y=Counts, x=Regulation, label=Counts)) + 
    geom_bar(position="stack", stat="identity") +
    geom_text(colour="white", size = 5, position = position_stack(vjust = 0.5)) + 
      xlab("")
  return(p)
}

sig_genes_plot <- function(results, baseMean_cutoff=50, lFC_cutoff=0.58, padj_cutoff=0.05) {

  results <- as.data.frame(results)
  results <- results[results["baseMean"] > baseMean_cutoff,]
  results$calls <- 0
  results$calls[which(results$log2FoldChange > lFC_cutoff
                        & results$padj < padj_cutoff)] <- 1
  results$calls[which(results$log2FoldChange < -lFC_cutoff 
                        & results$padj < padj_cutoff)] <- -1
  Up <- apply(results, 2, function(x) sum(x == 1))
  Down <- apply(results, 2, function(x) sum(x == -1))
  stats <- rbind(Up, Down)
  gg <- reshape2::melt(stats)
  gg <- gg[11:length(rownames(gg)),]
  colnames(gg) <- c("Regulation", "Comparisons", "Genes")
  gg <- gg[gg["Comparisons"] == "calls",]
  print(gg)
  plot_bar <- ggplot2::ggplot(
    gg,
    ggplot2::aes(x = Comparisons, y = Genes, fill = Regulation)
  ) +
    ggplot2::geom_bar(position = "dodge", stat = "identity", width = 0.8) +
    #ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "top",
      axis.title.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
    ) +
    ggplot2::ylab("Number of differntially expressed genes") +
    ggplot2::geom_text(
      ggplot2::aes(label = Genes),
      position = ggplot2::position_dodge(width = 0.9),
      vjust = -0.5,
      hjust = 0
    ) +
    ggplot2::ylim(0, max(gg$Genes) * 1.1)
  
  return(plot_bar)
}


# Manual function for Normalization as done in DESeq
# “median of ratios normalization”
mor_normalization = function(data){
  library(dplyr)
  library(tibble)
  
  # take the log
  log_data = log(data) 
  
  # find the psuedo-references per sample by taking the geometric mean
  log_data = log_data %>% 
    rownames_to_column('gene') %>% 
    mutate (gene_averages = rowMeans(log_data)) %>% 
    filter(gene_averages != "-Inf")
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_data)
  
  # where to stop before the pseudo column 
  before_pseduo = pseudo_column - 1
  
  # find the ratio of the log data to the pseudo-reference
  ratios = sweep(log_data[,2:before_pseduo], 1, log_data[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  return(manually_normalized)
}


# Differential genes visualisation

plotDispEsts <- function( cds , main = "")
{
  
  plot(
    rowMeans( counts( cds, normalized=TRUE ) ),
    fitInfo(cds)$perGeneDispEsts,
    main = main, pch = 16, log="xy", cex.main = 4, cex = 2, cex.lab = 4, cex.axis = 4, mgp = c(6,0,0), mkh = 100)
  xg <- 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" , lwd = 3)
}



plotDE <- function( res, main = "" ){
  plot(
    res$baseMean,
    res$log2FoldChange,
    log="x", pch=16, , main = main,cex.main = 2, cex = 2, cex.lab = 4, cex.axis = 4, mgp = c(6,0,0), mkh = 10
    ,col = ifelse( res$padj < .05, "red", "black" ) )
}


#Enrichment Visualization

plot_GO <- function(enrichedGO_results, show=10, filename="", mutate=F, title="") {
  
  if (dim(enrichedGO_results)[1] != 0) {
    
    if (mutate) {
      fit <- mutate(enrichedGO_results, qscore = -log(p.adjust, base=10)) %>% 
        barplot(x="qscore", title=title)
    } else{
      fit <- plot(barplot(enrichedGO_results, showCategory = show, title=title))
    }
    if(str_length(filename) > 0) {
      png(filename, res=300, width=3000, height=3000)
      print(fit)
      dev.off()
    }
    
  }

}

# get normalized gene_counts
# example genelist_pattern >>> c("Myh11", "Cnn1", "Mylk", "Myocd", "Flna", 
# "Fn1", "Col[0-9]a*", "Col[0-9][0-9]a*")
get_normalized_counts_of_genes <- function(genelist_patterns, df_trans_data, gene_details){
  normalized_counts <- df_trans_data %>% 
    rownames_to_column("ENSEMBL") %>% 
    merge(y=gene_details, by="ENSEMBL", all.x=TRUE)
  normalized_counts <- normalized_counts %>% column_to_rownames("SYMBOL")
  df <- data.frame()
  for(gene in genelist_patterns){
    df <- rbind(df, normalized_counts[grep(gene,rownames(normalized_counts)),])
    }
  return(df)
}
  
get_enriched_geneset <- function(enriched_data, enrich_ID, diff_gene_list_for_kegg_gsea){
  # enriched data:  is obtained from gseKEGG function
  # ID: check for the ID of the pathway to plot
  # diff_gene_list_for_kegg_gsea: table obtained after merging enselbl mart anf diff_exp_genes
  
  df <- as.data.frame(enriched_data)
  ID_row <- df[df$ID == enrich_ID,]
  entrez_ids <- ID_row$core_enrichment %>% 
    strsplit("\\/")
  df_sliced <- diff_gene_list_for_kegg_gsea %>% 
    filter(diff_gene_list_for_kegg_gsea$ENTREZ_ID %in% unlist(entrez_ids))
  return(df_sliced$SYMBOL)
}

get_enriched_geneset_of_gsego <- function(enriched_data, description){
  # enriched data:  is obtained from gseGofunction
  # Description: check for the Descrition of the pathway to plot
  
  # return the list of genes in the descritption
  
  df <- as.data.frame(enriched_data)
  ID_row <- df[df$Description == description,]
  gene_ids <- ID_row$geneID %>% strsplit("\\/")
  print(description)
  return(gene_ids)
}