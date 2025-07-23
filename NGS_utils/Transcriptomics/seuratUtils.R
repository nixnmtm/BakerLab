# Single Cell/Spatial Utility functions
# Author: Nixon Raj

source("~/research/coding/BakerLab/NGS_utils/Transcriptomics/RNASeqUtils.R")
####Seurat Visium####

plot_average_exp_HeatMap <- function(obj, cluster_colname){
  
  n = length(levels(obj[[cluster_colname]][[cluster_colname]]))
  datalist = vector("list", length = n)
  for (i in 1:n){
    exp_data <- as.data.frame(get_highly_exp_genes(obj, cluster=i-1, genenames = F))
    exp_data$gene <- rownames(exp_data)
    exp_data <- exp_data %>% gather(key='sample', value='value', -gene)
    datalist[[i]] <- exp_data
  }
  
  alldata = do.call(rbind, datalist)
  ggplot(alldata, aes(sample, gene)) + geom_tile(aes(fill=value))
}  

library(viridis)
plotHeatMap <- function(df){
  og <- df
  df$gene <- rownames(df)
  df <- df %>% gather(key='cluster', value='value', -gene)
  ggplot(df, aes(cluster, gene, fill=value)) + 
    geom_tile() + 
    scale_fill_viridis()
}

plot_ComplexHeatMap <- function(obj, features, metadata_cluster_colname=NULL, cluster_color_map=NULL) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(scales)
  
  
  # Check which markers are missing in scale.data
  missing_genes <- setdiff(features, rownames(obj[["RNA"]]@layers$scale.data))
  
  if(length(missing_genes) == 0){
    expr_mat <- GetAssayData(obj, layer = "scale.data")[features, ]
  } else {
    message("Scaling raw data for missing genes: ", paste(missing_genes, collapse=", "))
    mat <- obj[["RNA"]]$data[features, ] %>% as.matrix()
    expr_mat <- t(scale(t(mat)))
  }
  
  # Replace extreme scaled values to avoid color compression
  expr_mat[expr_mat > 3] <- 3
  expr_mat[expr_mat < -3] <- -3
  
  # Extract cluster annotation
  cluster_anno <- obj@meta.data[[metadata_cluster_colname]]
  cluster_levels <- levels(factor(cluster_anno))
  
  # If no cluster_color_map provided, generate distinct colors automatically
  if (is.null(cluster_color_map)) {
    cluster_colors <- hue_pal()(length(cluster_levels))
    names(cluster_colors) <- cluster_levels
  } else {
    # Ensure provided mapping covers all clusters
    missing_clusters <- setdiff(cluster_levels, names(cluster_color_map))
    if (length(missing_clusters) > 0) {
      stop("The following clusters are missing in cluster_color_map: ", paste(missing_clusters, collapse = ", "))
    }
    cluster_colors <- cluster_color_map[cluster_levels]
  }
  
  # Create a diverging color function centered at 0
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  # Create a column annotation for clusters
  column_ha <- HeatmapAnnotation(
    cluster = factor(cluster_anno, levels = cluster_levels),
    col = list(cluster = cluster_colors),
    show_annotation_name = FALSE
  )
  
  # Plot heatmap
  Heatmap(
    expr_mat,
    name = "Expression",
    col = col_fun,
    top_annotation = column_ha,
    column_split = factor(cluster_anno, levels = cluster_levels),
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    cluster_column_slices = FALSE,
    column_title_gp = gpar(fontsize = 8),
    column_gap = unit(0.5, "mm"),
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_title_rot = 90,
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_quality = 4
  )
}



check_save_dea_data <- function(markers, path, filename, format="csv"){
  
  if (file.exists(filename)) {
    #Delete file if it exists
    file.remove(filename)
    write.table(markers, file = file.path(path,paste0(filename, ".", format)),
                sep = ",",
                append = F,
                col.names=T, 
                row.names = F)
  }else{
    write.table(markers, file = file.path(path,paste0(filename, ".", format)),
                sep = ",",
                append = F,
                col.names=T, 
                row.names = F)
  }
}


getNwriteDEG_df <- function(markers, path=NULL, file_name=NULL, FDR=F, pcut=1e-2, FCcut=1, 
                            rankbyFC=F,
                            rankbyPval=F,
                            rankbyAdjPval=F,
                            rankbyPctDiff=T){
  #:markers: results of FindMarkers
  #:path: relative path to save
  #:filename: filename to save
  
  ge <- markers %>% mutate(pct.diff=pct.1-pct.2)
  
  suffix="unfiltered"
  if (!is.null(path) & !is.null(file_name)){
    check_save_dea_data(markers = ge, path = path, 
                        filename = paste0("DEG_",file_name, "-", suffix))
  }
  
  suffix="filtered"
  if (FDR){
    ge <- ge %>%
      filter(p_val_adj_fdr < pcut & abs(avg_log2FC) > FCcut)
    print(paste0("Filtering by FDR adjusted pval cutoff: ", pcut))
    print(paste0("and average log2FC cutoff: ", FCcut))
  }else{
    ge <- ge %>%
      filter(p_val_adj < pcut & abs(avg_log2FC) > FCcut)
    print(paste0("Filtering by adjusted pval cutoff: ", pcut))
    print(paste0("and average log2FC cutoff: ", FCcut))
  }

  if (rankbyPctDiff){
    ge <- arrange(ge, desc(pct.diff))
  }
  
  if (rankbyFC){
    ge <- arrange(ge, desc(avg_log2FC))
  }
  
  if (rankbyPval){
    ge <- arrange(ge, p_val)
  }
  
  if (rankbyAdjPval){
    ge <- arrange(ge, p_val_adj)
  }
  
  if (!is.null(path) & !is.null(file_name)){
    check_save_dea_data(markers = ge, path = path, 
                        filename = paste0("DEG_",file_name,
                                          "_pcut-",pcut,"_FCcut-", FCcut, 
                                          "_", suffix))
  }
  return(ge)
}

compare_cond_clusters <- function(object, combined_ident_col, ident1, ident2) {
  # Check if combined_ident_col exists in metadata
  if (!(combined_ident_col %in% colnames(object@meta.data))) {
    stop(paste0("The combined identity column ", combined_ident_col, " is not found in object metadata. 
                Please create it before calling this function.
                Example: bl.int$cond_cluster <- paste(bl.int$histology, bl.int$collapsed_clusters, sep = '_'),
                         In this case the combined_ident_col is 'cond_cluster'"))
  }
  
  # Subset object to groups of interest
  cells_use <- rownames(object@meta.data)[object@meta.data[[combined_ident_col]] %in% c(ident1, ident2)]
  sub_obj <- subset(object, cells = cells_use)
  print(table(sub_obj[[combined_ident_col]]))
  
  # Set identity to combined identity column
  Idents(sub_obj) <- combined_ident_col
  
  # Run DE between groups
  markers <- FindMarkers(sub_obj, ident.1 = ident1, ident.2 = ident2)
  
  return(list(markers = markers, object = sub_obj))
}

getDEA_plotVP <- function(object, setIdent=NULL, ident1=NULL, ident2=NULL,  
                          pcutoff=1e-2, FCcutoff=1,
                          filename=NULL, filepath=NULL,
                          plotFDR=F, groupby=NULL,
                          connectors=F, subset_ident=NULL,
                          deg_col=c("firebrick1", "royalblue", "black"),
                          markGenes=NULL,
                          use_combined_ident=FALSE,
                          combined_ident_col=NULL) {
  
  
  
  
  
  if(use_combined_ident) {
    if(is.null(combined_ident_col) | is.null(ident1) | is.null(ident2)) {
      stop("Please provide combined_ident_col, group1 and group2 when use_combined_ident=TRUE")
    }
    
    # Call the compare_cond_clusters helper
    de_results <- compare_cond_clusters(object, combined_ident_col, ident1, ident2)
    print(de_results$markers)
    markers <- de_results$markers
    object_sub <- de_results$object
  } else {
    if (is.null(ident1) | is.null(ident2)){
      stop("Please provide both ident1 and ident2")
    }
    
    if (is.null(subset_ident) && is.null(groupby)){
      Idents(object) <- setIdent
      markers <- FindMarkers(object, ident.1 = ident1, ident.2= ident2)
      object_sub <- object
    } else if (!is.null(groupby)) {
      Idents(object) <- object[[groupby, drop=TRUE]]
      markers <- FindMarkers(object, ident.1 = ident1, ident.2 = ident2)
      object_sub <- object
    } else if (!is.null(subset_ident)){
      message("Subsetting based on cluster-", subset_ident)
      object_sub <- subset(object, idents=subset_ident)
      Idents(object_sub) <- setIdent
      markers <- FindMarkers(object_sub, ident.1 = ident1, ident.2 = ident2)
    }
  }
  
  # p-val adjustment and rest of plotting code...
  markers$p_val_adj_fdr <- p.adjust(markers$p_val, method='fdr')
  markers <- markers %>% as.data.frame() %>% rownames_to_column("gene")

    keyvals <- rep(deg_col[3], nrow(markers))  # Default color

    if (plotFDR){
      y = "p_val_adj_fdr"
      dir.create(file.path(filepath, "FDR"), recursive = TRUE)
      filepath <- file.path(filepath, "FDR")
      filename = paste0(filename, "_FDR")
      keyvals[markers$avg_log2FC > FCcutoff &
                markers$p_val_adj_fdr < pcutoff] <- deg_col[1]
      keyvals[markers$avg_log2FC < -FCcutoff &
                markers$p_val_adj_fdr < pcutoff] <- deg_col[2]

      getNwriteDEG_df(markers, path=filepath, FDR=T, file_name=filename,
                      pcut=pcutoff, FCcut=FCcutoff,
                      rankbyFC=F,rankbyPval=T,
                      rankbyAdjPval=F,rankbyPctDiff=F)

    }else{
      y = "p_val_adj"
      dir.create(file.path(filepath, "noFDR"), recursive = TRUE)
      filepath <- file.path(filepath, "noFDR")
      filename = paste0(filename, "_noFDR")
      keyvals[markers$avg_log2FC > FCcutoff &
                markers$p_val_adj < pcutoff] <- deg_col[1]
      keyvals[markers$avg_log2FC < -FCcutoff &
                markers$p_val_adj < pcutoff] <- deg_col[2]
      getNwriteDEG_df(markers, path=filepath, FDR=F, file_name=filename,
                      pcut=pcutoff, FCcut=FCcutoff,
                      rankbyFC=F,rankbyPval=T,
                      rankbyAdjPval=F,rankbyPctDiff=F)
    }

    suptitle = paste0("The dotted lines indicate","\n",
                      "Pval cutoff: ", pcutoff, "\n",
                      "FoldChange cutoff: ", FCcutoff)


    keyvals[is.na(keyvals)] <- deg_col[3]
    names(keyvals)[keyvals == deg_col[2]] <- paste0("Down-regulated","\n", ident2)
    names(keyvals)[keyvals == deg_col[1]] <- paste0("Up-regulated", "\n", ident1)
    names(keyvals)[keyvals == deg_col[3]] <- "NS"

    if (is.null(markGenes)){
      p <- EnhancedVolcano(markers,
                           lab=markers$gene,
                           title = filename,
                           subtitle = suptitle,
                           x='avg_log2FC',
                           y=y, FCcutoff = FCcutoff,
                           ylab = bquote('-'~Log[10]~ 'adjusted p_value'),
                           pCutoffCol = y,
                           pCutoff = pcutoff,
                           xlab = bquote('Average' ~Log[2]~ 'fold change'),
                           labSize = 4.0,
                           pointSize = 4,
                           colAlpha = 0.8,
                           legendLabSize = 12,
                           legendIconSize = 2.0,
                           widthConnectors = 0.75,
                           drawConnectors = connectors, arrowheads = F,
                           max.overlaps=10, colCustom = keyvals)
    }else{
      p <- EnhancedVolcano(markers, lab=markers$gene,
                           selectLab=markGenes,
                           title = filename,
                           subtitle = suptitle,
                           x='avg_log2FC',
                           y=y, FCcutoff = FCcutoff,
                           ylab = bquote('-'~Log[10]~ 'adjusted p_value'),
                           pCutoffCol = y,
                           pCutoff = pcutoff,
                           xlab = bquote('Average' ~Log[2]~ 'fold change'),
                           labSize = 4.0,
                           pointSize = 4,
                           colAlpha = 0.8,
                           legendLabSize = 12,
                           legendIconSize = 2.0,
                           widthConnectors = 0.75,
                           drawConnectors = connectors, arrowheads = F,
                           max.overlaps=10, colCustom = keyvals)
    }
    save_it(p, filepath, paste0("VP-", filename, "_pcut-", pcutoff,"_FCcut-", FCcutoff),
            format = "png", resolution=300, w=1200, h=1200)
    return(markers)
  }

# Function to exclude reads from cells with given UMIs > 1000
remove_contaminated_spots <- function(obj, genes2check, umi_count_cutoff=1000){
  
  obj$UMI_counts <- Matrix::colSums(GetAssayData(obj, 
                                                     assay = "Spatial", 
                                                     layer = "counts")[genes2check, , drop = FALSE])
  # Remove cells/spots with > 1000 UMIs of given genes
  obj <- subset(obj, subset = UMI_counts <= umi_count_cutoff)
  return(obj)
}

# Function to remove unwanted genes from a Seurat object
remove_unwanted_genes <- function(obj) {
  # Define lists of genes to remove
  hb_genes <- c("HBA1", "HBA2", "HBB", "HBG1", "HBG2", "HBD", "HBM")
  mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
  ribo_genes <- grep("^RPL|^RPS", rownames(obj), value = TRUE)
  
  # Combine all unwanted genes into a single list
  unwanted_genes <- unique(c(hb_genes, mt_genes, ribo_genes))
  
  # Get the counts matrix from the Spatial assay
  counts <- GetAssayData(obj, assay = "Spatial")
  
  # Print dimensions of the original counts matrix
  message("Original counts matrix dimensions: ", dim(counts))
  
  # Identify the rows (genes) that match the unwanted gene names
  unwanted_gene_indices <- which(rownames(counts) %in% unwanted_genes)
  
  # Print the number of unwanted genes found and their names
  message("Number of unwanted genes found: ", length(unwanted_gene_indices))
  message("Unwanted genes found: ", rownames(counts)[unwanted_gene_indices])
  
  # Remove unwanted genes from the counts matrix
  counts_filtered <- counts[-unwanted_gene_indices, ]
  
  # Print dimensions of the filtered counts matrix
  message("Filtered counts matrix dimensions: ", dim(counts_filtered))
  
  # Update the counts matrix in the original object
  obj[['Spatial']]@counts <- counts_filtered
  obj[['Spatial']]@data <- counts_filtered
  
  # Remove unwanted genes from the features list
  obj <- subset(obj, features = rownames(counts_filtered))
  
  return(obj)
}

process_Seurat <- function(obj, npcs=20, dims=20, res=0.8){
  #"Run all steps involved in a Seurat analysis"
  
  obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE) %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = T, npcs = npcs) %>%
    FindNeighbors(dims = 1:dims) %>%
    FindClusters(verbose = T, resolution = res) %>%
    RunUMAP(dims = 1:dims)# %>%
  #RunTSNE(dims = 1:dims)
  return(obj)
}

get_stDE_data <- function (markers, 
                           cluster_id = F, 
                           p_val_adj_cut = 1e-2, 
                           avg_log2FC_cut = 1,
                           top=10){
  if (cluster_id){
    markers <- markers[markers$cluster == cluster_id,]
  }
  sorted_slice <- markers %>%
    filter(p_val_adj < p_val_adj_cut & abs(avg_log2FC) > avg_log2FC_cut) %>% 
    arrange(p_val_adj, -avg_log2FC)
  rbind(slice_head(sorted_slice, n = top), 
        slice_tail(sorted_slice, n = top)) -> sliced_top
  return(sliced_top)
}

plot_HVF <- function(obj, top=50) {
  obj <- FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = 2000)
  # Identify the 10 most highly variable genes
  HVF_top <- head(x = VariableFeatures(object = wt), top)
  plot1 <- VariableFeaturePlot(object = obj)
  LabelPoints(plot = plot1, points = HVF_top, repel = TRUE, xnudge = 0, ynudge = 0)
  
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

get_highly_exp_genes <- function(obj, top_n=20, cluster=0, genenames=T, assay="RNA"){
  
  # To get the highly expressed genes in each cluster and get there
  avg.exp <- AverageExpression(obj, group.by = "seurat_clusters")
  avg.exp.df <- avg.exp[[assay]] %>% as.data.frame()
  genes <- apply(avg.exp.df[paste0("g",cluster)], 2, sort, decreasing=T) %>% head(n=top_n)
  if(genenames){
    return(rownames(genes))
  }else{
    return(genes)
  }
}

plotAverageExp_between_clusters <- function(obj, x_clus=NULL, y_clus=NULL, 
                                            assay="RNA", top=10){
  
  # compare avergae expression between cluster of an seurat object
  
  library(cowplot)
  theme_set(theme_cowplot())
  avg.exp <- log1p(AverageExpression(obj, verbose = FALSE)[[assay]]) %>% as.data.frame()
  
  xcol <- paste0("g",x_clus)
  ycol <- paste0("g",y_clus)
  x <- avg.exp[[xcol]]
  y <- avg.exp[[ycol]]
  data <- data.frame(x,y)
  rownames(data) <- rownames(avg.exp)
  data$diff <- data$x- data$y
  message("The overall highly expressed genes in this two clusters")
  #mx <- data[with(data, order(-pmax(x, y))),]
  mx <- data %>% filter(x > 2.5 & y > 2.5) %>% arrange(abs(diff))
  high_in_both <- mx %>% head(top)
  print(high_in_both)
  t <- arrange(data, desc(diff)) %>% head(top) %>% rownames()
  b <- arrange(data, diff) %>% head(top) %>% rownames()
  genes.to.label <- c(t,b)
  corr <- cor(data$x, data$y)
  p1 <- ggplot(as.matrix(data), aes(x,y)) + 
    geom_point() + geom_smooth(method = "lm") + 
    geom_abline() +
    scale_x_continuous(n.breaks = 10)  +
    scale_y_continuous(n.breaks = 10) +
    ggtitle(paste0("Correlation-Clus", 
                   x_clus, "=x",
                   " vs Clus",
                   y_clus, "=y",
                   ": ",
                   round(corr, 3))) +
    geom_text(x = 2.5, y = 0, label = lm_eqn(data), parse = TRUE)
  
  p1 <- LabelPoints(plot = p1, points = t, repel = TRUE, color="red")
  p1 <- LabelPoints(plot = p1, points = b, repel = TRUE, color="blue")
  p1 <- LabelPoints(plot = p1, points = rownames(high_in_both), repel = TRUE, color="green")
  p1
}

library(ggplot2)
library(RColorBrewer)
library(dplyr)  # For data manipulation

# Function to generate a consistent color palette for cell types
generate_color_palette <- function(cell_types) {
  # Generate a color palette based on the number of unique cell types
  num_cell_types <- length(unique(cell_types))
  
  # Use Set1 or another suitable palette for categorical data
  colors <- brewer.pal(min(num_cell_types, 12), "Set1")
  
  # Create a named vector of colors for the cell types
  color_mapping <- setNames(colors[1:num_cell_types], unique(cell_types))
  
  return(color_mapping)
}

anchorMapping <- function(reference, query, query.dims=15, anchor.labels,save.loc=FALSE,
                          metadata_colname="subclass",
                          normalization = "LogNormalize", reduction.method="cca"){
  DefaultAssay(query) <- "Spatial"
  anchors = FindTransferAnchors(reference, query = query,  normalization.method = normalization,
                                reduction = reduction.method)
  predictions.assay <- TransferData(anchorset = anchors, refdata = reference[[metadata_colname, drop = TRUE]], 
                                    prediction.assay=TRUE, weight.reduction=reduction.method, dims=1:query.dims)
  non.mapping <- c()
  for(i in 1:dim(predictions.assay)[1]){ if(sum(predictions.assay@data[i,])==0) 
    non.mapping <- c(non.mapping, rownames(predictions.assay[i]))}
  predictions.assay@misc$non.mapping <- non.mapping
  predictions.assay@misc$mapping <- setdiff(anchor.labels, non.mapping)
  query[["predictions"]] <- predictions.assay
  DefaultAssay(query) <- 'predictions'
  return(query)
}

total_celltype_proportion <- function(sp.obj, assay_name = "predictions", 
                                      layer_name = "data", 
                                      plot = FALSE, color_mapping = NULL,
                                      savepath=NULL, prefix=NULL) {
  
  #' Calculate and Plot Total Cell Type Proportions
  #'
  #' Computes the proportions of cell types from a spatial transcriptomics object that has been anchor-mapped 
  #' with a single-cell reference and includes a "predictions" assay. Optionally generates a plot to visualize 
  #' these proportions.
  #'
  #' @param sp.obj A spatial transcriptomics object containing anchor-mapped data with single-cell reference 
  #'               and a "predictions" assay.
  #' @param assay_name Character. Name of the assay from which to retrieve data (default is "predictions").
  #' @param layer_name Character. Name of the data layer in the assay (default is "data").
  #' @param plot Logical. If TRUE, generates a plot showing cell type proportions (default is FALSE).
  #' @param color_mapping Vector of colors for cell types. If NULL, colors are auto-generated.
  #' @param savepath Character. Path to save the plot image (optional).
  #' @param prefix Character. Prefix for the saved file name if savepath is provided (optional).
  #'
  #' @return A data frame with calculated proportions of each cell type.
  #'
  #' @details The function retrieves cell type prediction data from a spatial object that has been anchor-mapped 
  #'          to a single-cell reference dataset and includes a "predictions" assay. It calculates the proportions 
  #'          of each cell type and, if `plot = TRUE`, creates a stacked horizontal bar plot to visualize them.
  #'
  #' @examples
  #' total_celltype_proportion(sp.obj = my_spatial_data, plot = TRUE, savepath = "output")
  #'
  
  # Retrieve the prediction data
  pred <- GetAssayData(object = sp.obj, assay = assay_name, layer = layer_name)
  
  # Check if pred is a valid object
  if (is.null(pred) || nrow(pred) == 0) {
    stop("No prediction data found.")
  }
  
  # Exclude the last row (assumed to be max values)
  proportion_matrix <- pred[1:(nrow(pred) - 1), , drop = FALSE]
  
  # Calculate total proportions
  cell_type_totals <- rowSums(proportion_matrix)
  total_proportions <- cell_type_totals / sum(cell_type_totals)
  
  if (plot) {
    # Prepare the data for plotting
    
    proportion_df <- data.frame(
      CellType = names(total_proportions),
      Proportion = total_proportions,
      Percentage = paste0(round(total_proportions, 3)*100, "%")  # Convert to percentage for display
    )# %>%
    #arrange(desc(Proportion))  # Arrange in descending order by proportion
    
    proportion_df$Label <- paste0(proportion_df$Percentage, "(",proportion_df$CellType,")")
    
    # Reorder the CellType factor based on total proportions
    proportion_df$CellType <- factor(proportion_df$CellType, levels = proportion_df$CellType)
    
    # Check if proportion_df is valid
    if (nrow(proportion_df) == 0) {
      stop("No data available for plotting.")
    }
    
    # Generate color mapping if not provided
    if (is.null(color_mapping)) {
      color_mapping <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(proportion_df$CellType)))
    }
    
    # Create the stacked horizontal bar graph
    p <- ggplot(proportion_df, aes(x = "", y = Proportion, fill = CellType)) +
      geom_bar(stat = "identity", position = "stack") +                    # Create stacked bars
      geom_text(aes(label = Label),                                         # Add proportion labels
                position = position_stack(vjust = 0.5),                  # Position labels in the middle
                color = "black", size = 4, angle=90) +                             # Customize text color and size
      ggtitle("Proportion of Cell Types Across All Spatial Spots") + 
      labs(x = "", y = paste0("Proportion (%)")) +                                 # Axis labels
      theme_minimal() +                                                   # Simple theme
      scale_fill_manual(values = color_mapping) +                        # Use generated colors
      coord_flip() +                                                      # Flip coordinates for horizontal bars
      theme(legend.title = element_blank())                                # Remove legend title for clarity
    
    if(!is.null(savepath)){
      save_it(p, filepath=savepath, filename = paste0(prefix, "_total_celltype_proportions"), 
              format = "png", w=3000, h=1500)
    }
    
    
    # Print the plot
    print(p)
  }
  
  return(proportion_df)  # Return the calculated proportions
}

is_not_empty <- function(df) {
  !is.null(df) && nrow(df) > 0
}
#### GENE SET ENRICHMENT GO & KEGG ####
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
check_for_empty_enrichment <- function(enrichment_results) {
  if (length(enrichment_results) == 0 || nrow(enrichment_results) == 0) {
    message("No enrichment results found. Skipping this step.")
    return(FALSE)
  }
  return(TRUE)
}

# Function to map gene symbols to Entrez IDs
get_gentrez2gene_map_list <- function(genename_list) {
  # Get ENTREZID for the provided gene symbols
  entrez2gene_df <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = genename_list,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "SYMBOL"
  )
  
  # Remove duplicates and NA values
  entrez2gene_df <- entrez2gene_df[!duplicated(entrez2gene_df), ]
  entrez2gene_df <- na.omit(entrez2gene_df)
  
  # Map SYMBOL to ENTREZID
  map_entrez2gene <- entrez2gene_df$ENTREZID
  names(map_entrez2gene) <- entrez2gene_df$SYMBOL
  
  return(map_entrez2gene)
}

# Function to safely apply gene mapping to the dataframe
apply_entrez_mapping <- function(df, gene_column) {
  # Get the entrez mapping for the given gene list (from the gene column)
  entrez_map <- get_gentrez2gene_map_list(df[[gene_column]])
  
  # Match the SYMBOL from the dataframe with the names in the mapping
  df$entrez_id <- sapply(df[[gene_column]], function(gene) entrez_map[gene])
  
  return(df)
}

prepare_gene_lists <- function(dea_markers, pcut = 1e-2, FCcut = 1) {
  # Filter and prepare gene list
  deg <- dea_markers %>%
    as.data.frame() %>%
    filter(p_val_adj < pcut, abs(avg_log2FC) > FCcut) %>%
    arrange(p_val_adj)
  
  upregulated_genes <- deg[deg$avg_log2FC > 0,]$gene
  downregulated_genes <- deg[deg$avg_log2FC < 0,]$gene
  
  # Convert to ENTREZ IDs
  geneid.ls <- list(
    upregulated = upregulated_genes,
    downregulated = downregulated_genes
  )
  
  # Convert genes to ENTREZIDs
  geneid.ls <- lapply(geneid.ls, function(genes) {
    gene.df <- AnnotationDbi::select(org.Mm.eg.db, keys = genes, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
    unique(na.omit(gene.df$ENTREZID))
  })
  
  return(geneid.ls)
}


run_KEGG_enrichment <- function(geneid.ls, numCategory, prefix, filename, savepath) {
  compKEGG <- compareCluster(
    geneCluster = geneid.ls,
    fun = "enrichKEGG",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    organism = "mmu"
  )
  
  if (!check_for_empty_enrichment(compKEGG)) {
    return(NULL)
  }
  
  # Dot plot for KEGG
  g4 <- dotplot(compKEGG, showCategory = numCategory, title = paste0(prefix, "-", filename, "-KEGG Pathway Enrichment Analysis"))
  save_it(g4, savepath, paste0(filename, "-enrichKEGG"), format = "png", resolution = 300, w = 800, h = 1000)
  return(compKEGG)
}


run_GO_enrichment <- function(geneid.ls, ont, pvalueCutoff = 0.05, numCategory = 10,
                              filename = NULL, outpath = NULL) {
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  
  enrichment_results <- enrichGO(
    gene = geneid.ls, 
    OrgDb = "org.Mm.eg.db",  
    ont = ont, 
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "BH"
  )
  
  if (!check_for_empty_enrichment(enrichment_results@result)) {
    return(NULL)
  }
  
  # Calculate FoldEnrichment as per PMID: 34557778
  GeneRatio_freq <- as.numeric(sub("/\\d+", "", enrichment_results@result$GeneRatio)) / 
    as.numeric(sub("^\\d+/", "", enrichment_results@result$GeneRatio))
  BgRatio_freq <- as.numeric(sub("/\\d+", "", enrichment_results@result$BgRatio)) / 
    as.numeric(sub("^\\d+/", "", enrichment_results@result$BgRatio))
  
  enrichment_results@result <- enrichment_results@result %>%
    mutate(FoldEnrichment = GeneRatio_freq / BgRatio_freq)
  
  # Calculate Rich factor
  enrichment_results@result = mutate(enrichment_results@result, 
                       richFactor = Count/as.numeric(sub("/\\d+", "", BgRatio)))
  
  g <- dotplot(enrichment_results, showCategory = numCategory,
               title = paste0(filename, "-GO Enrichment Analysis"))
  
  save_it(g, outpath, paste0(filename, "-enrichGO"),
          format = "png", resolution = 300, w = 800, h = 1000)
  
  return(enrichment_results@result %>% head(numCategory))
}


run_CNET_plot <- function(go_result, fc, numCategory, ontology, prefix, filename, savepath) {
  if (!check_for_empty_enrichment(go_result)) return(NULL)
  
  # Get upregulated and downregulated genes
  upregulated_genes <- names(fc)[fc > 0]
  downregulated_genes <- names(fc)[fc < 0]
  
  # Create the CNET plot for upregulated genes
  cnet_up <- cnetplot(
    go_result,
    showCategory = numCategory,
    color.params = list(foldChange = fc[upregulated_genes], category = "black")
  ) + ggtitle(paste0(prefix, "-", filename, "Upregulated_GO_", ontology))
  
  # Save the CNET plot for upregulated genes
  save_it(cnet_up, savepath, paste0(filename, "_cnetplot-GO_", ontology, "_Upregulated"), format = "png", resolution = 300, w = 1000, h = 1500)
  
  # Create the CNET plot for downregulated genes
  cnet_down <- cnetplot(
    go_result,
    showCategory = numCategory,
    color.params = list(foldChange = fc[downregulated_genes], category = "black")
  ) + ggtitle(paste0(prefix, "-", filename, "Downregulated_GO_", ontology))
  
  # Save the CNET plot for downregulated genes
  save_it(cnet_down, savepath, paste0(filename, "_cnetplot-GO_", ontology, "_Downregulated"), format = "png", resolution = 300, w = 1000, h = 1500)
}


run_GO_KEGGenrichment <- function(dea_markers, pcut = 1e-2, FCcut = 1, numCategory = 10, savepath = NULL, filename = NULL, prefix = NULL) {
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  
  # Prepare gene lists
  geneid.ls <- prepare_gene_lists(dea_markers, pcut, FCcut)
  
  results <- list()
  
  # Run KEGG enrichment
  results$kegg <- run_KEGG_enrichment(geneid.ls, numCategory, prefix, filename, savepath)
  
  
  
  # Run GO enrichment and get the results
  results$ALL_up <- run_GO_enrichment(geneid.ls = geneid.ls$upregulated,
                                      ont="ALL", pvalueCutoff = pcut,
                                      numCategory = numCategory, 
                                      paste0("ALL_UP_",filename), 
                                      savepath)
  
  results$ALL_down <- run_GO_enrichment(geneid.ls = geneid.ls$downregulated,
                                        ont="ALL", pvalueCutoff = pcut,
                                        numCategory = numCategory, 
                                        paste0("ALL_DOWN_",filename), 
                                        savepath)

  results$BP_up <- run_GO_enrichment(geneid.ls$upregulated, 
                                     ont="BP", pvalueCutoff = pcut,
                                     numCategory = numCategory, 
                                     paste0("BP_UP_",filename), 
                                     savepath)

  results$BP_down <- run_GO_enrichment(geneid.ls$downregulated, 
                                       ont="BP", pvalueCutoff = pcut,
                                       numCategory = numCategory, 
                                       paste0("BP_DOWN_",filename), 
                                       savepath)

  results$CC_up <- run_GO_enrichment(geneid.ls$upregulated, 
                                     ont="CC", pvalueCutoff = pcut,
                                     numCategory = numCategory, 
                                     paste0("CC_UP_",filename), 
                                     savepath)

  results$CC_down <- run_GO_enrichment(geneid.ls$downregulated, 
                                       ont="CC", pvalueCutoff = pcut,
                                       numCategory = numCategory, 
                                       paste0("CC_DOWN_",filename),
                                       savepath)

  results$MF_up <- run_GO_enrichment(geneid.ls$upregulated, 
                                     ont="MF", pvalueCutoff = pcut,
                                     numCategory = numCategory, 
                                     paste0("MF_UP_",filename), 
                                     savepath)

  results$MF_down <- run_GO_enrichment(geneid.ls$downregulated, 
                                       ont="MF", pvalueCutoff = pcut,
                                       numCategory = numCategory, 
                                       paste0("MF_DOWN_",filename), 
                                       savepath)
  print(paste0("ALL DONE"))
  return(results)
}

# Load necessary libraries
library(ggplot2)

plot_top_GO_enrichment_from_df <- function(go_results_df,
                                           numTopTerms = 10,
                                           savepath = NULL,
                                           plot_by_x = "FoldEnrichment") {
  # Required columns
  required_cols <- c("qvalue", "Description", "GeneRatio", "BgRatio", "Count")
  if (!all(required_cols %in% colnames(go_results_df))) {
    stop(paste("Missing required columns:", paste(setdiff(required_cols, colnames(go_results_df)), collapse = ", ")))
  }
  
  # Calculate metrics
  go_results_df <- go_results_df %>%
    mutate(
      GeneRatio = as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("^\\d+/", "", GeneRatio)),
      BgRatio = as.numeric(sub("/\\d+", "", BgRatio)) / as.numeric(sub("^\\d+/", "", BgRatio)),
      richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)),
      FoldEnrichment = GeneRatio / BgRatio
    )
  
  # Select top terms based on qvalue
  top_go_terms <- go_results_df %>%
    arrange(qvalue) %>%
    slice_head(n = min(numTopTerms, nrow(.))) %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  
  # Y axis mapping and label
  y_col <- switch(plot_by_x,
                  "GeneRatio" = "GeneRatio",
                  "FoldEnrichment" = "FoldEnrichment",
                  "richFactor" = "richFactor",
                  stop("Invalid 'plot_by_x' value: must be 'GeneRatio', 'FoldEnrichment', or 'richFactor'."))
  
  y_label <- switch(plot_by_x,
                    "GeneRatio" = "Gene Ratio",
                    "FoldEnrichment" = "Fold Enrichment",
                    "richFactor" = "Rich Factor")
  
  # Define a threshold for label placement
  threshold <- 0.1 * max(top_go_terms[[y_col]], na.rm = TRUE)
  
  # Create a new column for hjust and text color based on value size
  top_go_terms <- top_go_terms %>%
    mutate(
      hjust_value = ifelse(.data[[y_col]] < threshold, -0.1, 1.05),
      text_color = ifelse(.data[[y_col]] < threshold, "black", "white")
    )
  
  # Build plot
  p <- ggplot(top_go_terms, aes(x = .data[[y_col]], 
                                y = reorder(Description, .data[[y_col]]), 
                                fill = qvalue)) +
    geom_bar(stat = "identity") +
    
    # geom_text(aes(label = Description, 
    #               hjust = hjust_value, 
    #               color = text_color),
    #           size = 5, show.legend = FALSE) +
    # 
    scale_color_identity() +  # Use colors as-is
    scale_fill_gradient(
      low = "red", high = "blue", name = "FDR",
      guide = guide_colourbar(reverse = TRUE)
    ) +
    
    theme_minimal() +
    labs(x = y_label, y = NULL) +
    theme(
      panel.grid.major = element_blank(),   # Remove major grid lines
      panel.grid.minor = element_blank(),   # Remove minor grid lines
      panel.background = element_blank(),   # Remove panel background
      axis.line = element_line(color = "black", linewidth = 2),  # Keep x and y axis lines
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_line(color = "black", linewidth = 1.2),# Show axis ticks
      axis.ticks.length.y = unit(-0.25, "cm"),  # Y-axis ticks inward
      axis.ticks.length.x = unit(0.25, "cm"),   # X-axis ticks default outward
    )
  return(p)
}

run_integration_analysis <- function(obj, integration, result_path, res = 2, fcut = 1, pcut = 0.05, cols=NULL) {
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(EnhancedVolcano)
    
    # ---- Step 1: Integration and Clustering ----
    obj <- FindNeighbors(obj, reduction = integration, dims = 1:30)
    obj <- FindClusters(obj, resolution = res, cluster.name = paste0(integration, "_clusters"))
    obj <- RunUMAP(obj, dims = 1:30, reduction = integration, reduction.name = paste0("umap.", integration))
    obj <- RunTSNE(obj, dims = 1:30, reduction = integration, reduction.name = paste0("tsne.", integration))
    obj <- JoinLayers(obj)
    
    # ---- Step 2: Define Colors ----
    if (!is.null(cols)){
      color.use <- cols 
      names(color.use) <- levels(obj)
    }else{
      color.use=NULL
    }
    
    # ---- Step 3: Dimensionality Reduction Plots ----
    p1 <- DimPlot(obj, reduction = paste0("umap.", integration), group.by = paste0(integration, "_clusters"),
                  combine = FALSE, label.size = 2, cols = color.use)
    p2 <- SpatialDimPlot(obj, pt.size.factor = 10, group.by = paste0(integration, "_clusters"),
                         combine = FALSE, label.size = 2, cols = color.use)
    
    wrapped_plot <- wrap_plots(c(p1, p2), ncol = 2, byrow = FALSE)
    
    # ---- Step 4: Create Directory for Results ----
    integration_save_path <- file.path(result_path, paste0(integration))
    dir.create(integration_save_path, recursive = TRUE)
    
    # Save wrapped plot
    save_it(wrapped_plot, integration_save_path, paste0(integration, "_res", res, "_UMAP_Spatial"), 
            format = "png", resolution = 300, w = 3000, h = 5000)

    # ---- Step 5: Identify Cluster Markers ----
    markers <- FindAllMarkers(obj) %>%
      mutate(pct.diff = pct.1 - pct.2) %>%
      group_by(cluster)
    
    # ---- Step 6: Extract and Print Top 5 Markers per Cluster ----
    top5 <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) %>%
      slice_head(n = 5) %>%
      ungroup()
    
    print(top5)  # Print top 5 markers in the console
    
    # Save top 5 markers as CSV
    write.csv(top5, file = file.path(integration_save_path, paste0(integration, "_res", res, "_Top5Markers.csv")))
    
    # ---- Step 7: Plot Complex Heatmap ----
    hm <- plot_ComplexHeatMap(obj, markers = top5, metadata_cluster_colname = paste0(integration, "_clusters"))
    save_it(hm, integration_save_path, paste0(integration, "_res", res,"_HeatMap"), 
            format = "png", resolution = 300, w = 3000, h = 5000)
    
    # ---- Step 8: Volcano Plots and DEG Export ----
    clusters <- unique(markers$cluster)
    
    for (i in clusters) {
      dea_clus <- markers[markers$cluster == i, ]
      suptitle <- paste0("Cluster ", i, " vs Others")
      
      p <- EnhancedVolcano(
        dea_clus,
        lab = dea_clus$gene,
        title = "Integrated Bladders",
        subtitle = suptitle,
        x = 'avg_log2FC',
        y = 'p_val_adj',
        FCcutoff = fcut,
        ylab = "p_val_adj",
        pCutoffCol = 'p_val_adj',
        pCutoff = pcut,
        xlab = bquote('Average' ~ Log[2] ~ 'fold change'),
        labSize = 5.0,
        pointSize = 3,
        colAlpha = 0.8,
        legendLabSize = 12,
        legendIconSize = 2.0,
        widthConnectors = 0.75,
        gridlines.major = FALSE,
        gridlines.minor = FALSE,
        drawConnectors = TRUE,
        max.overlaps = 20
      )
      
      save_it(p, integration_save_path, paste0(integration, "_res", res, "_Cluster", i, "vsOthers"), 
              format = "png", resolution = 300, w = 3000, h = 5000)
      
      # ---- Step 9: Filter DEGs and Save CSV ----
      filtered_dea <- markers %>%
        filter(avg_log2FC > 1, p_val_adj < 1e-2, cluster == i)
      
      write.csv(
        x = filtered_dea,
        file = file.path(integration_save_path, 
                         paste0("Filtered_", integration, "_res", res, "_Cluster", i, "vsOthers", "_DEG.csv"))
      )
    }
  return(obj)
}



# barplot_cell_proportion
# Cell_proportion
# https://github.com/Alexis-Varin/RightOmicsTools/blob/main/R/Barplot_Cell_Proportion.R

Barplot_Cell_Proportion = function(seurat_object,
                                   group.by = NULL,
                                   split.by = NULL,
                                   idents = NULL,
                                   group.idents = NULL,
                                   split.idents = NULL,
                                   colors = NULL,
                                   order.prop = NULL,
                                   order.group = NULL,
                                   order.split = NULL,
                                   order.colors = TRUE,
                                   alpha = 1,
                                   show.cellsum.label = TRUE,
                                   cellsum.label.size = 3,
                                   axis.text.size = 9,
                                   x.axis.angle = 60,
                                   x.axis.hjust = 1,
                                   y.axis.title.size = 11,
                                   legend.text.size = 9,
                                   legend.side = "bottom",
                                   show.legend = TRUE,
                                   split.plot.title.size = 24,
                                   prop.percent = TRUE,
                                   nrow = 1,
                                   unique.group.plot = TRUE,
                                   unique.split.plot = FALSE,
                                   output.data = FALSE) {
  
  ident1 = ident2 = sumpercent = nbcells = NULL
  
  if (is.null(group.by)) {
    prop.percent = FALSE
  }
  
  if (is.character(group.idents) & !is.character(group.by)) {
    group.idents = NULL
  }
  if (is.character(split.idents) & !is.character(split.by)) {
    split.idents = NULL
  }
  
  if (is.character(idents)) {
    if (is.character(group.idents)) {
      if (is.character(split.idents)) {
        activeident = Idents(seurat_object)
        seurat_object@meta.data$id_group_split = paste(as.character(seurat_object@active.ident),
                                                       as.character(seurat_object[[group.by]][ , 1, drop = TRUE]),
                                                       as.character(seurat_object[[split.by]][ , 1, drop = TRUE]), sep = "_")
        idents = paste(rep(idents, each = length(group.idents)), group.idents, sep = "_")
        idents = paste(rep(idents, each = length(split.idents)), split.idents, sep = "_")
        Idents(seurat_object) = "id_group_split"
        seurat_object = subset(seurat_object, idents = idents)
        Idents(seurat_object) = activeident
      }
      else {
        activeident = Idents(seurat_object)
        seurat_object@meta.data$id_group = paste(as.character(seurat_object@active.ident),
                                                 as.character(seurat_object[[group.by]][ , 1, drop = TRUE]), sep = "_")
        idents = paste(rep(idents, each = length(group.idents)), group.idents, sep = "_")
        Idents(seurat_object) = "id_group"
        seurat_object = subset(seurat_object, idents = idents)
        Idents(seurat_object) = activeident
      }
    }
    else if (is.character(split.idents)) {
      activeident = Idents(seurat_object)
      seurat_object@meta.data$id_split = paste(as.character(seurat_object@active.ident),
                                               as.character(seurat_object[[split.by]][ , 1, drop = TRUE]), sep = "_")
      idents = paste(rep(idents, each = length(split.idents)), split.idents, sep = "_")
      Idents(seurat_object) = "id_split"
      seurat_object = subset(seurat_object, idents = idents)
      Idents(seurat_object) = activeident
    }
    else {
      seurat_object = subset(seurat_object, idents = idents)
    }
  }
  else if (is.character(group.idents)) {
    if (is.character(split.idents)) {
      activeident = Idents(seurat_object)
      seurat_object@meta.data$group_split = paste(as.character(seurat_object[[group.by]][ , 1, drop = TRUE]),
                                                  as.character(seurat_object[[split.by]][ , 1, drop = TRUE]), sep = "_")
      group.idents = paste(rep(group.idents, each = length(split.idents)), split.idents, sep = "_")
      Idents(seurat_object) = "group_split"
      seurat_object = subset(seurat_object, idents = group.idents)
      Idents(seurat_object) = activeident
    }
    else {
      activeident = Idents(seurat_object)
      Idents(seurat_object) = as.character(seurat_object[[group.by]][ , 1, drop = TRUE])
      seurat_object = subset(seurat_object, idents = group.idents)
      Idents(seurat_object) = activeident
    }
  }
  else if (is.character(split.idents)) {
    activeident = Idents(seurat_object)
    Idents(seurat_object) = as.character(seurat_object[[split.by]][ , 1, drop = TRUE])
    seurat_object = subset(seurat_object, idents = split.idents)
    Idents(seurat_object) = activeident
  }
  
  if (is.null(colors)) {
    colors = hue_pal()(n = length(levels(Idents(seurat_object))))
  }
  
  if (isTRUE(order.colors)) {
    if (!is.null(names(colors))) {
      colors = colors[levels(Idents(seurat_object))]
    }
    if (is.character(order.prop)) {
      if (length(order.prop) > 1) {
        names(colors) = levels(Idents(seurat_object))
        colors = colors[order.prop]
      }
      else {
        if (order.prop == "reverse") {
          colors = rev(colors)
        }
      }
    }
  }
  
  idents.df = data.frame("ident1" = Idents(seurat_object))
  
  if (is.character(group.by)) {
    Idents(seurat_object) = group.by
    idents.df$ident2 = Idents(seurat_object)
  }
  
  table.list = proportion.plot = list()
  
  if (is.character(split.by)) {
    Idents(seurat_object) = split.by
    idents.df$ident3 = Idents(seurat_object)
    if (is.null(order.split)) {
      order.split = levels(Idents(seurat_object))
    }
    else {
      if (length(order.split) == length(levels(Idents(seurat_object)))) {
        seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.split)
      }
      else {
        if (length(order.split) == 1 & any(order.split == "reverse")) {
          order.split = rev(levels(Idents(seurat_object)))
          seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.split)
        }
        else {
          stop("order.split needs to be either 'reverse' or a character or numeric vector of same length as the number of identities")
        }
      }
    }
    levels.split.by = levels(Idents(seurat_object))
    for (i in levels.split.by) {
      Idents(seurat_object) = split.by
      split1.df = idents.df[which(idents.df$ident3 == i),]
      if (is.character(group.by)) {
        Idents(seurat_object) = group.by
        if (is.null(order.group)) {
          order.group = levels(Idents(seurat_object))
        }
        if (!is.null(order.group)) {
          if (length(order.group) > 1) {
            seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.group)
          }
          else {
            if (order.group == "reverse") {
              order.group = rev(levels(Idents(seurat_object)))
              seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.group)
            }
            else {
              stop("order.group needs to be either 'reverse' or a character vector")
            }
          }
        }
        k = l = 1
        table.df = sum.df = data.frame()
        for (j in levels(Idents(seurat_object))) {
          split2.df = split1.df[which(split1.df$ident2 == j),]
          if (!is.null(order.prop)) {
            if (is.character(order.prop)) {
              if (length(order.prop) > 1) {
                table.list[[i]][[j]] = as.data.frame(table(split2.df$ident1))
                colnames(table.list[[i]][[j]]) = c("ident1","nbcells")
                table.list[[i]][[j]]$percent = table.list[[i]][[j]]$nbcells/sum(table.list[[i]][[j]]$nbcells)
                table.list[[i]][[j]]$ident1 = factor(table.list[[i]][[j]]$ident1, levels = order.prop)
              }
              else {
                if (order.prop == "reverse") {
                  table.list[[i]][[j]] = as.data.frame(rev(table(split2.df$ident1)))
                  colnames(table.list[[i]][[j]]) = c("ident1","nbcells")
                  table.list[[i]][[j]]$percent = table.list[[i]][[j]]$nbcells/sum(table.list[[i]][[j]]$nbcells)
                }
                else {
                  stop("order.prop needs to be either 'reverse' or a character vector")
                }
              }
            }
            else {
              stop("order.prop needs to be either 'reverse' or a character vector")
            }
          }
          else {
            table.list[[i]][[j]] = as.data.frame(table(split2.df$ident1))
            colnames(table.list[[i]][[j]]) = c("ident1","nbcells")
            table.list[[i]][[j]]$percent = table.list[[i]][[j]]$nbcells/sum(table.list[[i]][[j]]$nbcells)
          }
          table.list[[i]][[j]]$ident2 = j
          l = l+1
          if (isTRUE(unique.group.plot & sum(table.list[[i]][[j]]$nbcells) > 0)) {
            table.df = rbind(table.df, table.list[[i]][[j]])
            k = k+1
            sum.df = rbind(sum.df,data.frame("ident1" = i,
                                             "sum" = sum(table.list[[i]][[j]]$nbcells),
                                             "ident2" = j,
                                             "sumpercent" = sum(table.list[[i]][[j]]$percent)))
            
          }
          else {
            if (isFALSE(unique.group.plot)) {
              if (isTRUE(prop.percent)) {
                if (isTRUE(j == levels(Idents(seurat_object))[1])) {
                  if (isTRUE(show.cellsum.label)) {
                    proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=percent, x=ident2, fill = ident1)) +
                      geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                        x = j,
                                                                                                        y="Relative number of cells") +
                      theme_bw() +
                      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                            axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                            axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                            axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                            plot.margin = margin(5.5,5.5,5.5,20.5))+
                      scale_fill_manual(values=colors)+
                      geom_text(data = data.frame("ident1" = i,
                                                  "sum" = sum(table.list[[i]][[j]]$nbcells),
                                                  "ident2" = j,
                                                  "sumpercent" = sum(table.list[[i]][[j]]$percent)), aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                                fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                                inherit.aes = F)+
                      scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)
                  }
                  else {
                    proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=percent, x=ident2, fill = ident1)) +
                      geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                        x = j,
                                                                                                        y="Relative number of cells") +
                      theme_bw() +
                      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                            axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                            axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                            axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                            plot.margin = margin(5.5,5.5,5.5,20.5))+
                      scale_y_continuous(expand= c(0,0), labels = percent)+
                      scale_fill_manual(values=colors)
                  }
                }
                else {
                  if (isTRUE(show.cellsum.label)) {
                    proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=percent, x=ident2, fill = ident1)) +
                      geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                        x = j,
                                                                                                        y="Relative number of cells") +
                      theme_bw() +
                      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                            axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                            axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                            axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                            plot.margin = margin(5.5,5.5,5.5,20.5))+
                      scale_fill_manual(values=colors)+
                      geom_text(data = data.frame("ident1" = i,
                                                  "sum" = sum(table.list[[i]][[j]]$nbcells),
                                                  "ident2" = j,
                                                  "sumpercent" = sum(table.list[[i]][[j]]$percent)), aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                                fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                                inherit.aes = F)+
                      scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), breaks = NULL)
                  }
                  else {
                    proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=percent, x=ident2, fill = ident1)) +
                      geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                        x = j,
                                                                                                        y="Relative number of cells") +
                      theme_bw() +
                      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                            axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                            axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                            axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                            plot.margin = margin(5.5,5.5,5.5,20.5))+
                      scale_y_continuous(expand= c(0,0), breaks = NULL)+
                      scale_fill_manual(values=colors)
                  }
                }
              }
              else {
                if (isTRUE(show.cellsum.label)) {
                  proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=nbcells, x=ident2, fill = ident1)) +
                    geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                      x = j,
                                                                                                      y="Number of cells") +
                    theme_bw() +
                    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                          axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                          axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                          axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side)+
                    scale_fill_manual(values=colors)+
                    geom_text(data = data.frame("ident1" = i,
                                                "sum" = sum(table.list[[i]][[j]]$nbcells),
                                                "ident2" = j,
                                                "sumpercent" = sum(table.list[[i]][[j]]$percent)), aes(label = sum, x = ident2, y = sum), stat = "summary",
                              fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                              inherit.aes = F)+
                    scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))
                }
                else {
                  proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=nbcells, x=ident2, fill = ident1)) +
                    geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                      x = j,
                                                                                                      y="Number of cells") +
                    theme_bw() +
                    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                          axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                          axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                          axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side)+
                    scale_y_continuous(expand= c(0,0))+
                    scale_fill_manual(values=colors)
                }
              }
            }
          }
        }
        if (isTRUE(unique.group.plot)) {
          if (isTRUE(prop.percent)) {
            if (isTRUE(show.cellsum.label)) {
              proportion.plot[[i]] = ggplot(table.df, aes(fill=ident1, y=percent, x=factor(ident2, levels = order.group))) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                  y="Relative number of cells") +
                geom_text(data = sum.df, aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                          fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                          inherit.aes = F)+
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                      legend.position = legend.side, aspect.ratio = l/k, plot.margin = margin(5.5,65.5,5.5,5.5))+
                scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)+
                scale_fill_manual(values=colors)+
                ggtitle(i)
            }
            else {
              proportion.plot[[i]] = ggplot(table.df, aes(fill=ident1, y=percent, x=factor(ident2, levels = order.group))) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                  y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                      legend.position = legend.side, aspect.ratio = l/k, plot.margin = margin(5.5,65.5,5.5,5.5))+
                scale_y_continuous(expand= c(0,0), labels = percent)+
                scale_fill_manual(values=colors)+
                ggtitle(i)
            }
          }
          else {
            if (isTRUE(show.cellsum.label)) {
              proportion.plot[[i]] = ggplot(table.df, aes(fill=ident1, y=nbcells, x=factor(ident2, levels = order.group))) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                  y="Number of cells") +
                geom_text(data = sum.df, aes(label = sum, x = ident2, y = sum), stat = "summary",
                          fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                          inherit.aes = F)+
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                      legend.position = legend.side, aspect.ratio = l/k, plot.margin = margin(5.5,65.5,5.5,5.5))+
                scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))+
                scale_fill_manual(values=colors)+
                ggtitle(i)
            }
            else {
              proportion.plot[[i]] = ggplot(table.df, aes(fill=ident1, y=nbcells, x=factor(ident2, levels = order.group))) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                  y="Number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                      legend.position = legend.side, aspect.ratio = l/k, plot.margin = margin(5.5,65.5,5.5,5.5))+
                scale_y_continuous(expand= c(0,0))+
                scale_fill_manual(values=colors)+
                ggtitle(i)
            }
          }
        }
      }
      else {
        if (!is.null(order.prop)) {
          if (is.character(order.prop)) {
            if (length(order.prop) > 1) {
              table.list[[i]] = as.data.frame(table(split1.df$ident1))
              colnames(table.list[[i]]) = c("ident1","nbcells")
              table.list[[i]]$percent = table.list[[i]]$nbcells/sum(table.list[[i]]$nbcells)
              table.list[[i]]$ident1 = factor(table.list[[i]]$ident1, levels = order.prop)
            }
            else {
              if (order.prop == "reverse") {
                table.list[[i]] = as.data.frame(rev(table(split1.df$ident1)))
                colnames(table.list[[i]]) = c("ident1","nbcells")
                table.list[[i]]$percent = table.list[[i]]$nbcells/sum(table.list[[i]]$nbcells)
              }
              else {
                stop("order.prop needs to be either 'reverse' or a character vector")
              }
            }
          }
          else {
            stop("order.prop needs to be either 'reverse' or a character vector")
          }
        }
        else {
          table.list[[i]] = as.data.frame(table(split1.df$ident1))
          colnames(table.list[[i]]) = c("ident1","nbcells")
          table.list[[i]]$percent = table.list[[i]]$nbcells/sum(table.list[[i]]$nbcells)
        }
        table.list[[i]]$ident3 = i
        
        if (isTRUE(prop.percent)) {
          if (isTRUE(show.cellsum.label)) {
            proportion.plot[[i]] = ggplot(table.list[[i]], aes(y=percent, x=ident1, fill = ident1)) +
              geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Relative number of cells") +
              geom_text(aes(label = nbcells, x = ident1, y = percent), stat = "identity",
                        vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size)+
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                    legend.position = legend.side, plot.margin = margin(5.5,65.5,5.5,5.5))+
              scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)+
              scale_fill_manual(values=colors)+
              NoLegend()+
              ggtitle(i)
          }
          else {
            proportion.plot[[i]] = ggplot(table.list[[i]], aes(y=percent, x=ident1, fill = ident1)) +
              geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Relative number of cells") +
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                    legend.position = legend.side, plot.margin = margin(5.5,65.5,5.5,5.5))+
              scale_y_continuous(expand= c(0,0), labels = percent)+
              scale_fill_manual(values=colors)+
              NoLegend()+
              ggtitle(i)
          }
        }
        else {
          if (isTRUE(show.cellsum.label)) {
            proportion.plot[[i]] = ggplot(table.list[[i]], aes(y=nbcells, x=ident1, fill = ident1)) +
              geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Number of cells") +
              geom_text(aes(label = nbcells, x = ident1, y = nbcells), stat = "identity",
                        vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size)+
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                    legend.position = legend.side, plot.margin = margin(5.5,65.5,5.5,5.5))+
              scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))+
              scale_fill_manual(values=colors)+
              NoLegend()+
              ggtitle(i)
          }
          else {
            proportion.plot[[i]] = ggplot(table.list[[i]], aes(y=nbcells, x=ident1, fill = ident1)) +
              geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Number of cells") +
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                    legend.position = legend.side, plot.margin = margin(5.5,65.5,5.5,5.5))+
              scale_y_continuous(expand= c(0,0))+
              scale_fill_manual(values=colors)+
              NoLegend()+
              ggtitle(i)
          }
        }
      }
    }
    
    if (isFALSE(unique.group.plot) & is.character(group.by)) {
      for (i in 1:length(table.list)) {
        k = 1
        table.tmp = table.list[[i]]
        for (j in 1:length(table.tmp)) {
          if (sum(table.tmp[[j]]$nbcells) == 0) {
            proportion.plot[[i]][[k]] = NULL
            table.list[[i]][[k]] = NULL
          }
          else
            k = k+1
        }
      }
      
      for (i in 1:length(proportion.plot)) {
        for (j in 2:length(proportion.plot[[i]])) {
          proportion.plot[[i]][[j]] = proportion.plot[[i]][[j]]+
            theme(axis.title.y = element_blank())
        }
      }
      if (isTRUE(unique.split.plot)) {
        Idents(seurat_object) = split.by
        temp.plot = list()
        temp.title = ""
        for (i in levels(Idents(seurat_object))) {
          temp.plot = c(temp.plot,proportion.plot[[i]])
          temp.title = paste0(temp.title,i)
          if (i != levels(Idents(seurat_object))[length(levels(Idents(seurat_object)))]) {
            temp.title = paste0(temp.title," vs ")
            temp.plot = c(temp.plot,list(plot_spacer()))
          }
        }
        temp.plot = wrap_plots(temp.plot, nrow = nrow, guides = "collect")+
          plot_annotation(title = temp.title, theme = theme(plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                                                            legend.position = legend.side))
        proportion.plot = temp.plot
        if (isFALSE(show.legend)) {
          proportion.plot = proportion.plot+
            plot_annotation(theme = theme(legend.position = "none"))
        }
      }
      else {
        Idents(seurat_object) = split.by
        for (i in levels(Idents(seurat_object))) {
          proportion.plot[[i]] = wrap_plots(proportion.plot[[i]], nrow = nrow, guides = "collect")+
            plot_annotation(title = i, theme = theme(plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                                                     legend.position = legend.side))
          if (isFALSE(show.legend)) {
            proportion.plot[[i]] = proportion.plot[[i]]+
              plot_annotation(theme = theme(legend.position = "none"))
          }
        }
      }
      if(isTRUE(output.data)) {
        return(table.list)
      }
      else {
        return(proportion.plot)
      }
    }
    
    else {
      if (isTRUE(unique.split.plot)) {
        Idents(seurat_object) = split.by
        proportion.plot[[length(levels(Idents(seurat_object)))]] = proportion.plot[[length(levels(Idents(seurat_object)))]]+
          theme(plot.margin = margin(5.5,5.5,5.5,5.5))
        proportion.plot = wrap_plots(proportion.plot, nrow = nrow, guides = "collect")+
          plot_annotation(theme = theme(legend.position = legend.side))
        if (isFALSE(show.legend)) {
          proportion.plot = proportion.plot+
            plot_annotation(theme = theme(legend.position = "none"))
        }
      }
      else {
        for (i in 1:length(proportion.plot)) {
          proportion.plot[[i]] = proportion.plot[[i]]+
            theme(plot.margin = margin(5.5,5.5,5.5,5.5))
          if (isFALSE(show.legend)) {
            proportion.plot[[i]] = proportion.plot[[i]]+
              NoLegend()
          }
        }
      }
      if(isTRUE(output.data)) {
        for (i in 1:length(table.list)) {
          k = 1
          table.tmp = table.list[[i]]
          for (j in 1:length(table.tmp)) {
            if (sum(table.tmp[[j]]$nbcells) == 0) {
              table.list[[i]][[k]] = NULL
            }
            else
              k = k+1
          }
        }
        return(table.list)
      }
      else {
        return(proportion.plot)
      }
    }
  }
  
  if (is.character(group.by)) {
    Idents(seurat_object) = group.by
    if (is.null(order.group)) {
      order.group = levels(Idents(seurat_object))
    }
    if (!is.null(order.group)) {
      if (length(order.group) > 1) {
        seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.group)
      }
      else {
        if (order.group == "reverse") {
          order.group = rev(levels(Idents(seurat_object)))
          seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.group)
        }
        else {
          stop("order.group needs to be either 'reverse' or a character vector")
        }
      }
    }
    table.df = sum.df = data.frame()
    for (j in levels(Idents(seurat_object))) {
      split1.df = idents.df[which(idents.df$ident2 == j),]
      if (!is.null(order.prop)) {
        if (is.character(order.prop)) {
          if (length(order.prop) > 1) {
            table.list[[j]] = as.data.frame(table(split1.df$ident1))
            colnames(table.list[[j]]) = c("ident1","nbcells")
            table.list[[j]]$percent = table.list[[j]]$nbcells/sum(table.list[[j]]$nbcells)
            table.list[[j]]$ident1 = factor(table.list[[j]]$ident1, levels = order.prop)
          }
          else {
            if (order.prop == "reverse") {
              table.list[[j]] = as.data.frame(rev(table(split1.df$ident1)))
              colnames(table.list[[j]]) = c("ident1","nbcells")
              table.list[[j]]$percent = table.list[[j]]$nbcells/sum(table.list[[j]]$nbcells)
            }
            else {
              stop("order.prop needs to be either 'reverse' or a character vector")
            }
          }
        }
        else {
          stop("order.prop needs to be either 'reverse' or a character vector")
        }
      }
      else {
        table.list[[j]] = as.data.frame(table(split1.df$ident1))
        colnames(table.list[[j]]) = c("ident1","nbcells")
        table.list[[j]]$percent = table.list[[j]]$nbcells/sum(table.list[[j]]$nbcells)
      }
      table.list[[j]]$ident2 = j
      if (isTRUE(unique.group.plot & sum(table.list[[j]]$nbcells) > 0)) {
        table.df = rbind(table.df, table.list[[j]])
        sum.df = rbind(sum.df,data.frame("sum" = sum(table.list[[j]]$nbcells),
                                         "ident2" = j,
                                         "sumpercent" = sum(table.list[[j]]$percent)))
        
      }
      else {
        if (isTRUE(prop.percent)) {
          if (isTRUE(j == levels(Idents(seurat_object))[1])) {
            if (isTRUE(show.cellsum.label)) {
              proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                  y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_fill_manual(values=colors)+
                geom_text(data = data.frame("sum" = sum(table.list[[j]]$nbcells),
                                            "ident2" = j,
                                            "sumpercent" = sum(table.list[[j]]$percent)), aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                          fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                          inherit.aes = F)+
                scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)
            }
            else {
              proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                  y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_y_continuous(expand= c(0,0), labels = percent)+
                scale_fill_manual(values=colors)
            }
          }
          else {
            if (isTRUE(show.cellsum.label)) {
              proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                  y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_fill_manual(values=colors)+
                geom_text(data = data.frame("sum" = sum(table.list[[j]]$nbcells),
                                            "ident2" = j,
                                            "sumpercent" = sum(table.list[[j]]$percent)), aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                          fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                          inherit.aes = F)+
                scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), breaks = NULL)
            }
            else {
              proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                  y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_y_continuous(expand= c(0,0), breaks = NULL)+
                scale_fill_manual(values=colors)
            }
          }
        }
        else {
          if (isTRUE(show.cellsum.label)) {
            proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=nbcells, x=ident2, fill = ident1)) +
              geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                y="Number of cells") +
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side)+
              scale_fill_manual(values=colors)+
              geom_text(data = data.frame("sum" = sum(table.list[[j]]$nbcells),
                                          "ident2" = j,
                                          "sumpercent" = sum(table.list[[j]]$percent)), aes(label = sum, x = ident2, y = sum), stat = "summary",
                        fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                        inherit.aes = F)+
              scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))
          }
          else {
            proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=nbcells, x=ident2, fill = ident1)) +
              geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                                y="Number of cells") +
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side)+
              scale_y_continuous(expand= c(0,0))+
              scale_fill_manual(values=colors)
          }
        }
      }
    }
    if (isFALSE(unique.group.plot)) {
      for (i in 2:length(proportion.plot)) {
        proportion.plot[[i]] = proportion.plot[[i]]+
          theme(axis.title.y = element_blank())
      }
      proportion.plot = wrap_plots(proportion.plot, nrow = nrow, guides = "collect")+
        plot_annotation(theme = theme(legend.position = legend.side))
      if (isFALSE(show.legend)) {
        proportion.plot = proportion.plot+
          plot_annotation(theme = theme(legend.position = "none"))
      }
      if(isTRUE(output.data)) {
        k = 1
        table.tmp = table.list
        for (j in 1:length(table.tmp)) {
          if (sum(table.tmp[[j]]$nbcells) == 0) {
            table.list[[k]] = NULL
          }
          else
            k = k+1
        }
        return(table.list)
      }
      else {
        return(proportion.plot)
      }
    }
    else {
      if (isTRUE(prop.percent)) {
        if (isTRUE(show.cellsum.label)) {
          proportion.plot = ggplot(table.df, aes(fill=ident1, y=percent, x=factor(ident2, levels = order.group))) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                              y="Relative number of cells") +
            geom_text(data = sum.df, aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                      fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                      inherit.aes = F)+
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side)+
            scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)+
            scale_fill_manual(values=colors)
        }
        else {
          proportion.plot = ggplot(table.df, aes(fill=ident1, y=percent, x=factor(ident2, levels = order.group))) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                              y="Relative number of cells") +
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side)+
            scale_y_continuous(expand= c(0,0), labels = percent)+
            scale_fill_manual(values=colors)
        }
      }
      else {
        if (isTRUE(show.cellsum.label)) {
          proportion.plot = ggplot(table.df, aes(fill=ident1, y=nbcells, x=factor(ident2, levels = order.group))) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                              y="Number of cells") +
            geom_text(data = sum.df, aes(label = sum, x = ident2, y = sum), stat = "summary",
                      fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                      inherit.aes = F)+
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side)+
            scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))+
            scale_fill_manual(values=colors)
        }
        else {
          proportion.plot = ggplot(table.df, aes(fill=ident1, y=nbcells, x=factor(ident2, levels = order.group))) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                              y="Number of cells") +
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side)+
            scale_y_continuous(expand= c(0,0))+
            scale_fill_manual(values=colors)
        }
      }
      if (isFALSE(show.legend)) {
        proportion.plot = proportion.plot+
          NoLegend()
      }
      if(isTRUE(output.data)) {
        return(table.df)
      }
      else {
        return(proportion.plot)
      }
    }
  }
  
  else {
    if (!is.null(order.prop)) {
      if (is.character(order.prop)) {
        if (length(order.prop) > 1) {
          table.df = as.data.frame(table(idents.df$ident1))
          colnames(table.df) = c("ident1","nbcells")
          table.df$percent = table.df$nbcells/sum(table.df$nbcells)
          table.df$ident1 = factor(table.df$ident1, levels = order.prop)
        }
        else {
          if (order.prop == "reverse") {
            table.df = as.data.frame(rev(table(idents.df$ident1)))
            colnames(table.df) = c("ident1","nbcells")
            table.df$percent = table.df$nbcells/sum(table.df$nbcells)
          }
          else {
            stop("order.prop needs to be either 'reverse' or a character vector")
          }
        }
      }
      else {
        stop("order.prop needs to be either 'reverse' or a character vector")
      }
    }
    else {
      table.df = as.data.frame(table(idents.df$ident1))
      colnames(table.df) = c("ident1","nbcells")
      table.df$percent = table.df$nbcells/sum(table.df$nbcells)
    }
    sum.df = data.frame("sum" = sum(table.df$nbcells),
                        "ident1" = levels(Idents(seurat_object)),
                        "sumpercent" = sum(table.df$percent))
    if (isTRUE(prop.percent)) {
      if (isTRUE(show.cellsum.label)) {
        proportion.plot = ggplot(table.df, aes(y=percent, x=ident1, fill = ident1)) +
          geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Relative number of cells") +
          geom_text(aes(label = nbcells, x = ident1, y = percent), stat = "identity",
                    vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size)+
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size),
                legend.position = legend.side)+
          scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)+
          scale_fill_manual(values=colors)+
          NoLegend()
      }
      else {
        proportion.plot = ggplot(table.df, aes(y=percent, x=ident1, fill = ident1)) +
          geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Relative number of cells") +
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size),
                legend.position = legend.side)+
          scale_y_continuous(expand= c(0,0), labels = percent)+
          scale_fill_manual(values=colors)+
          NoLegend()
      }
    }
    else {
      if (isTRUE(show.cellsum.label)) {
        proportion.plot = ggplot(table.df, aes(y=nbcells, x=ident1, fill = ident1)) +
          geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Number of cells") +
          geom_text(aes(label = nbcells, x = ident1, y = nbcells), stat = "identity",
                    vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size)+
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size),
                legend.position = legend.side)+
          scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))+
          scale_fill_manual(values=colors)+
          NoLegend()
      }
      else {
        proportion.plot = ggplot(table.df, aes(y=nbcells, x=ident1, fill = ident1)) +
          geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Number of cells") +
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size),
                legend.position = legend.side)+
          scale_y_continuous(expand= c(0,0))+
          scale_fill_manual(values=colors)+
          NoLegend()
      }
    }
    if (isFALSE(show.legend)) {
      proportion.plot = proportion.plot+
        NoLegend()
    }
    if(isTRUE(output.data)) {
      return(table.df)
    }
    else {
      return(proportion.plot)
    }
  }
}

addSmallLegend <- function(pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  # Reduces the font size of legend and the size of dots in the legend
  list(guides(
    shape = guide_legend(override.aes = list(size = pointSize)),
    color = guide_legend(override.aes = list(size = pointSize))
  ),
  theme(
    legend.title = element_text(size = textSize),
    legend.text = element_text(size = textSize),
    legend.key.size = unit(spaceLegend, "lines")
  ))
}

plot_avg_exp_barplot <-  function(avg.exp.mat){
  
  #Takes in a matrix with rownames and two columns
  #plots the expression into barplot
  
  df <- as.data.frame(as.matrix(avg.exp.mat))
  
  # Add gene names as a column
  df$gene <- rownames(df)
  
  # Load reshape2 for long-format conversion
  library(reshape2)
  df_long <- melt(df, id.vars = "gene", variable.name = "sample", value.name = "expression")
  
  # Load ggplot2 for plotting
  library(ggplot2)
  
  # Plot
  ggplot(df_long, aes(x = gene, y = expression, fill = sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    labs(title = "Gene expression comparison", x = "Gene", y = "Expression")
}



#' Modified from https://github.com/mojaveazure/seurat-disk/issues/172
#' Patch of SeuratDisk::SaveH5Seurat function
#'
#' The "Assay5" class attribute "RNA" needs to be converted to a standard "Assay"
#' class for compatibility with SeuratDisk. 
#'
#' @param object the Seurat object
#' @param filename the file path where to save the Seurat object
#' @param verbose SaveH5Seurat verbosity
#' @param overwrite whether to overwrite an existing file
#' 
#' @return NULL
SaveH5SeuratSpatialObject <- function(
    object,
    dims=30,
    filename,
    verbose = TRUE,
    overwrite = TRUE
) {
  
  object <- object %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs=dims) %>%
    RunUMAP(dims = 1:dims) 
  
  # add copy of "RNA" 
  object[["RNA"]] <- CreateAssayObject(counts = object[["Spatial"]]$counts)
  object@reductions$pca@cell.embeddings <- as.matrix(object@reductions$pca@cell.embeddings)
  object@reductions$umap@cell.embeddings <- as.matrix(object@reductions$umap@cell.embeddings)
  
  object@images[["slice1"]] <- object@images[["slice1"]]
  
  # switch default assay
  DefaultAssay(object) <- "RNA"
  # remove original
  #object[["RNA"]] <- NULL
  # export
  SaveH5Seurat(object, filename = filename, overwrite, verbose)
  
  return(NULL)
}


####SPATA utils#### 
convert2Seurat <- function (spata_obj, 
                            image_path, 
                            slice_name = "slice1", 
                            image_name = "tissue_lowres_image.png"){
  # convert a Spata object into seurat
  # using transformSpataToSeurat() doesnt add image to the seurat object
  # so created this function to manually add image and create a seurat obj
  
  seurat_obj <- SPATA2::transformSpataToSeurat(spata_obj, 
                                               NormalizeData=F, 
                                               FindVariableFeatures=F, 
                                               SCTransform = F,
                                               ScaleData=F, 
                                               RunPCA=F, 
                                               FindNeighbors=F, 
                                               FindClusters=F, 
                                               RunTSNE = F,
                                               RunUMAP = F)
  
  image_file <- file.path(image_path, image_name)
  print(image_file)
  # Ensure image_file is single string
  if (length(image_file) != 1) {
    warning("image_file has length > 1. Using the first entry.")
    image_file <- image_file[1]
  }
  
  if (file.exists(image_file)) { 
    seurat_obj <- base::tryCatch({
      lowres <- Read10X_Image(image_path, image.name = image_name)
      seurat_obj@images[[slice_name]] <- lowres
      seurat_obj@images[[slice_name]]@assay <- "Spatial"
      seurat_obj@images[[slice_name]]@key <- paste0(slice_name, "_")
      base::message("The SpatialImage is manually added to the Seurat Object as ", slice_name)
      seurat_obj
    }, error = function(error) {
      msg <- tryCatch(as.character(error), error=function(e) "Unknown error in error handler")
      base::warning("Error in adding image slice manually: ", paste(msg, collapse = " | "))
      seurat_obj
    })
  } else{
    warning("Image file does not exist: ", image_file)
  }
  return(seurat_obj)
}


get_all_histology_subset <- function(spata2_obj){
  "Create a subset based on barcodes present in organ segmentation defined, returns all segmentations except unnamed"
  fdata <- spata2_obj@fdata[[spata2_obj@samples]]
  segment_barcodes <- fdata[fdata$histology != "unnamed",]$barcodes
  subset_obj <- subsetByBarcodes(spata2_obj, barcodes = segment_barcodes, verbose = T)
  return(subset_obj)
}

get_organ_subset <- function(spata2_obj, histology_name){
  "Create a subset based on histology name given"
  fdata <- spata2_obj@fdata[[spata2_obj@samples]]
  if (isString(histology_name)){
    segment_barcodes <- fdata[fdata$histology == histology_name,]$barcodes
  } else {
    segment_barcodes <- fdata[fdata$histology %in% histology_name,]$barcodes
  }
  subset_obj <- subsetByBarcodes(spata2_obj, barcodes = segment_barcodes, verbose = T)
  return(subset_obj)
}

plot_surface_by_histology <- function(spata2_obj, histology_name=NULL) {
  # plotSurface based on the histology name/names provided
  if(is.null(histology_name)){
    plotSurface(spata2_obj, color_by="histology", pt_clrp = "lo", display_title=T)
  }else{
    spata2_subset <- get_organ_subset(spata2_obj, histology_name)
    plotSurface(spata2_subset, color_by="histology", pt_clrp = "lo", display_title=T)
  }
}

merge_histology <- function(spata2_obj, histology_names_to_merge, new_histology_name){
  # Merge histology types into a single type 
  # (eg: histology_names_to_merge = c(gut1, gut2, gut3) to new_histology_name = "Gut")
  
  fdata <- spata2_obj@fdata[[spata2_obj@samples]]
  paste("Before Merging:")
  print(unique(fdata$histology))
  fdata <- fdata %>% mutate(across(`histology`, as.character)) # change to char 
  fdata$histology[fdata$histology %in% histology_names_to_merge] <- new_histology_name
  
  fdata <- fdata %>% mutate(across(`histology`, as.factor)) # change it back to factor
  spata2_obj@fdata[[spata2_obj@samples]] <- fdata
  nw_fdata <- spata2_obj@fdata[[spata2_obj@samples]]
  paste("After Merging:")
  print(unique(nw_fdata$histology))
  return(spata2_obj)
}

get_organ_barcodes <- function(obj, histology_name){
  # Get barcode of the given histology type
  obj <- get_organ_subset(obj, histology_name = histology_name)
  obj_barcodes <- obj@fdata[[obj@samples]]$barcodes
  return(obj_barcodes)
}

calc_helper <- function(object,genes){
  counts = GetAssayData(object, assay = "RNA", layer = "counts")
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
    print(paste0("Number of cells expressing the gene", sum(counts[genes,]>0)))
    print(paste0("Total Number of cells", ncells))
  }else{return(NA)}
}

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

remove_HB_genes <- function(obj){
  # Remove all Hb genes from the object before any processing
  counts <- GetAssayData(obj, assay = "Spatial")
  print(dim(counts))
  print(length(which(grepl("Hb", rownames(counts)))))
  counts <- counts[-(which(grepl("Hb", rownames(counts)))),]
  print(dim(counts))
  obj <- subset(obj, features = rownames(counts))
  return(obj)
}



####  RCTD  ####

# copied the function fropm spacexr and added w and h
plot_cond_occur_mod <- function (cell_type_names, resultsdir, weights, puck, w=10, h=10, weight_cutoff=NULL) 
{
  occur <- numeric(length(cell_type_names))
  names(occur) = cell_type_names
  for (i in 1:length(cell_type_names)) {
    cell_type = cell_type_names[i]
    if (!is_null(weight_cutoff)){
      my_cond = weights[, cell_type] > weight_cutoff
    }
    else{
      my_cond = weights[, cell_type] > pmax(0.25, 2 - log(puck@nUMI,2) / 5)
    }
    occur[cell_type] = sum(my_cond)
  }
  df <- reshape2::melt(as.list(occur))
  colnames(df) = c("Count", "Cell_Type")
  pdf(file.path(resultsdir, "cell_type_occur.pdf"), width = w, height = h)
  plot <- ggplot2::ggplot(df, ggplot2::aes(x = Cell_Type, y = Count, 
                                           fill = Cell_Type)) + ggplot2::geom_bar(stat = "identity") + 
    ggplot2::theme_minimal()
  invisible(print(plot))
  dev.off()
  return(plot)
}

create_SpatialRNA_puck <- function(obj, layer="counts", image="s1"){
  # Get gene expression counts
  sp_counts <- LayerData(obj, assay = "Spatial", layer = layer)
  sp_coords <- GetTissueCoordinates(obj, image = image) 
  
  if ("cell" %in% names(sp_coords)){
    sp_coords <- sp_coords %>% select(-cell)
  }
  colnames(sp_coords) <- c("x", "y")
  sp_coords[is.na(colnames(sp_coords))] <- NULL
  
  # Create a query in spacexr
  query <- spacexr::SpatialRNA(sp_coords, sp_counts, colSums(sp_counts))
  query <- restrict_puck(query, colnames(query@counts))
  return(query)  
}


multi_rctd_2_cell_prop <- function(multi_rctd_obj){
  #' Convert Multi-Mode RCTD Results to a Cell Proportion Matrix
  #'
  #' This function processes the results from an RCTD object that was run in "MULTI" mode and 
  #' extracts the cell type proportions for each spot. It returns a matrix where rows represent
  #' spots and columns represent cell types, with values corresponding to the proportions of each
  #' cell type per spot.
  #'
  #' @param multi_rctd_obj An RCTD object that has been run in "MULTI" mode. 
  #'   This object contains the cell type information and results needed for the matrix generation.
  #'   For example, this can be created using `run.RCTD()` with `doublet_mode = "multi"`.
  #'
  #' @return A matrix of dimensions \code{total_spots} by \code{length(cell_types)}. The rows correspond to 
  #'   the spots in the spatial data, and the columns correspond to the cell types. Each entry in the matrix
  #'   is the proportion of a cell type for a given spot. Cell types that were not detected in a spot will have a value of 0.
  #'   The row names of the matrix correspond to the spot IDs, and the column names represent the cell types.
  #'
  #' @examples
  #' \dontrun{
  #' # Assuming `multi_myRCTD` is a multi-mode RCTD object:
  #' cell_prop_matrix <- multi_rctd_2_cell_prop(multi_myRCTD)
  #' head(cell_prop_matrix) # View the first few rows of the matrix
  #' }
  #'
  #' @export
  
  cell_types <- multi_rctd_obj@cell_type_info$info[[2]]
  total_spots <- length(multi_rctd_obj@results)
  
  # Create an empty matrix
  multi_mode_cell_prop <- matrix(0, nrow = total_spots, ncol = length(cell_types), 
                                 dimnames = list(NULL, cell_types))
  for (i in seq_along(multi_rctd_obj@results)) {
    sub_weights <- multi_rctd_obj@results[[i]]$sub_weights
    conf_list <- multi_rctd_obj@results[[i]]$conf_list
    
    # Match the names of sub_weights with the column names (cell_types)
    # and place the values in the corresponding columns of the matrix
    matching_indices <- match(names(conf_list), cell_types)
    
    # Fill the matrix at row 'i' with the sub_weights values in the matching columns
    multi_mode_cell_prop[i, matching_indices] <- sub_weights
  }
  
  rownames(multi_mode_cell_prop) <- colnames(multi_rctd_obj@spatialRNA@counts)
  return(multi_mode_cell_prop) 
  
}

normalize <- function(x){ 
  # normalize a vector
  return(x/sum(x)) 
}

RCTD_celltype_average_prop <- function(multi_rctd_obj, normalize=F, celltype=NULL){
  #' Calculate Average Proportion of Cell Types in a RCTD object 
  #'
  #' Computes the average proportion of cell types from a RCTD object (MULTI mode), optionally normalizing the 
  #' result and returning a specific cell type's proportion if specified.
  #'
  #' @param multi_rctd_obj A multi-sample RC-TD object containing cell type proportion data.
  #' @param normalize Logical. If TRUE, normalizes the average proportions (default is FALSE).
  #' @param celltype Character. Specific cell type to retrieve the average proportion for (optional).
  #'
  #' @return Numeric vector of average proportions for all cell types, or a single value if `celltype` is specified.
  #'
  #' @details This function retrieves cell type proportions across multiple samples, calculates the average 
  #'          proportions for each cell type, and can return the average of a specified cell type. Zeros are 
  #'          converted to NA for accurate averaging.
  #'
  #' @examples
  #' celltype_average_prop(multi_rctd_obj = my_multi_rctd_data, normalize = TRUE, celltype = "B_cells")
  #'
  #'
  
  # First get the celltype proportions
  multi_mode_cell_prop <- as.data.frame(multi_rctd_2_cell_prop(multi_rctd_obj))
  multi_mode_cell_prop[multi_mode_cell_prop == 0] <- NA 
  # Calculate average proportions, only pixels with expression of the celltype is considered
  avg_prop <- colMeans(multi_mode_cell_prop, na.rm = TRUE) %>% replace_na(0)
  if (normalize){
    avg_prop <- normalize(avg_prop)
  }
  if (is.null(celltype)){ 
    return(avg_prop)
  }
  else{
    return(avg_prop[[celltype]])
  }
}

#### STDeconvolve ####

seurat_to_spe <- function(seu, sample_id, img_id) {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(GetTissueCoordinates(seu)[,c("x","y")])
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}

plotNget_genenames <- function(deconGexp, deconProp, posi, 
                               ens2gene_map, cluster_id, top=20,
                               groups=NULL, group_cols=NULL, plot_it=F, print_it=T){
  # get top gene names from celltype clusters deconvoluted from STdeconvolve
  # also plot the Topic
  all_topgenes <- topGenes(deconGexp, n=top)
  topgenesC <- all_topgenes[[cluster_id]]
  geneNames <- c()
  if (all(startsWith(names(topgenesC), "ENS"))){
    for(ensid in names(topgenesC)){
      if (print_it){
        (ens2gene_map[[ensid]])
      }
      geneNames <- c(geneNames, ens2gene_map[[ensid]])
    }
  }else{
    if(print_it){
      print(names(topgenesC))
    }
    geneNames <- topgenesC
  }
  p <- vizTopic(theta = deconProp, pos = posi, topic = cluster_id, 
                plotTitle = paste0("Topic",cluster_id),
                size = 6, stroke = 0.5, alpha = 1,
                low = "white",
                high = "red")
  if (!is.null(groups)){
    p1 <- vizAllTopics(deconProp, posi[rownames(deconProp),], r=55,
                       groups = groups, 
                       group_cols = group_cols)
  }else{
    p1 <- vizAllTopics(deconProp, posi[rownames(deconProp),], r=55, lwd = 0)
  }
  print(geneNames)
  if(plot_it){
    p+p1
  }
  
  #return(geneNames)
}

get_geneSymbols_map <- function(se){
  geneSymbols <- se@rowRanges@elementMetadata$symbol
  names(geneSymbols) <- names(se@rowRanges)
  return(geneSymbols)
}

runLDA_STDeconvolve <- function(se, numb.of.celltypes=NULL, model.Savepath=NULL, filename=NULL){
  
  #geneSymbols <- get_geneSymbols_map(se)
  
  ## this is the genes x barcode sparse count matrix
  cd <- se@assays@data@listData$counts
  #rownames(cd) <- geneSymbols[rownames(cd)]
  pos <- SpatialExperiment::spatialCoords(se)
  colnames(pos) <- c('x', 'y')
  counts <- cleanCounts(cd , min.lib.size = 100, min.reads = 10)
  hist(log10(Matrix::colSums(counts)+1))
  overdis <- restrictCorpus(counts, removeAbove=0.95, removeBelow = 0.05, nTopOD = NA)
  ldas <- fitLDA(t(as.matrix(overdis)), Ks = numb.of.celltypes)
  
  # put all needed variables in a list to return for future use
  decon <- list()
  decon$counts <- counts
  decon$model <- ldas
  decon$pos <- pos
  decon$corpus <- overdis
  
  qs::qsave(ldas, file = file.path(model.Savepath, paste0(filename,"_model.qs")))
  qs::qsave(pos, file = file.path(model.Savepath, paste0(filename,"_pos.qs")))
  return(decon)
}

getCellProportions <- function(LDAmodel){
  ## extract deconvolved cell-type proportions (theta) 
  ## and transcriptional profiles (beta)
  optLDA <- optimalModel(models = LDAmodel, opt = "min")
  cellprops <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
  return(cellprops)
}

plotSTDeconvolve <- function(LDAmodel, pos=NULL, cols=NULL, file_path=NULL, 
                             file_name=NULL){
  results <- getCellProportions(LDAmodel)
  deconProp <- results$theta
  deconGexp <- results$beta
  
  if (is.null(cols)){
    cols=rainbow(ncol(deconProp))
  }
  
  if (!is.null(file_path)){
    p <- vizAllTopics(deconProp, pos[rownames(deconProp),], 
                      r=55, lwd=0, topicCols = cols) + 
      ggtitle(file_name) + 
      ggplot2::guides(colour = "none")
    save_it(p, filepath = file_path, filename=file_name, 
            format = "pdf", resolution=300)
  }else{
    vizAllTopics(deconProp, pos[rownames(deconProp),], 
                 r=55, lwd=0, topicCols = cols) + 
      ggtitle(file_name) + 
      ggplot2::guides(colour = "none")
  }
  
}

plotvizGeneCount <- function(se, genename=NULL, filepath=NULL) {
  cd <- se@assays@data@listData$counts
  #rownames(cd) <- geneSymbols[rownames(cd)]
  pos <- SpatialExperiment::spatialCoords(se)
  colnames(pos) <- c('x', 'y')
  counts <- cleanCounts(cd , min.lib.size = 100, min.reads = 10)
  df <- merge(as.data.frame(pos), as.data.frame(t(as.matrix(counts))), by = 0)
  if (!is.null(filepath)){
    p <- vizGeneCounts(df = df, gene = genename,
                       size = 3, stroke = 0.1, 
                       plotTitle = genename, winsorize = 0.05, 
                       showLegend = TRUE)
    save_it(p, filepath = filepath, 
            filename=paste0(genename,"_geneCount"), format = "pdf", resolution=300)
  }else{
    vizGeneCounts(df = df, gene = genename,
                  size = 3, stroke = 0.1, 
                  plotTitle = genename, winsorize = 0.05, 
                  showLegend = TRUE)
  }}


save_to_fetch_gene_details <- function(deconGexp, top=20, 
                                       savepath=NULL, 
                                       filename="deconGexp.csv"){
  df <- enframe(topGenes(deconGexp, n=top), 
                name = "list_id", 
                value = "values") %>%
    unnest_longer(values)
  colnames(df) <- c("topic", "deconGexp", "gene")
  df$topic <- as.numeric(df$topic)
  write.csv(df, 
            file = file.path(savepath, paste0(filename)),
            row.names = F, 
            quote = FALSE)
}


vlnplot_w_significance <- function(obj, gene_signature, file_name=NULL, group_by = NULL, 
                                   test_sign, y_max,w=14, h=8){
  plot_case <- function(signature){
    library(ggplot2)
    library(ggpubr)
    
    VlnPlot(obj, features = signature,
            pt.size =  0.5, #ncol=3,
            group.by = group_by, 
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
    ) + 
      stat_compare_means(comparisons = test_sign, label = "p.signif") + 
      NoLegend()
    #+stat_compare_means(comparisons = test_sign, label = "p.format")
  }
  print(map(gene_signature, plot_case) %>% cowplot::plot_grid(plotlist = .))
  if (!is_null(file_name)){
    file_name <- paste0(file_name, "_VlnPlot.png")
    ggsave(file_name, width = w, height = h)
  }
}

#### other common utils ####

# from Human to Mouse
convertHumanGeneList <- function(x){
  
  library("biomaRt")
  #human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  hs_mart <- useEnsembl("ensembl","hsapiens_gene_ensembl", mirror = "useast", host = "www.ensembl.org")
  mm_mart <- useEnsembl("ensembl","mmusculus_gene_ensembl", mirror = "useast", host = "www.ensembl.org")
  #mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                    values = x , mart = hs_mart, attributesL = c("mgi_symbol"), martL = mm_mart, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  return(humanx)
}

convert_genes <- function(gene_list, species = "human_to_mouse") {
  # Load the biomaRt package

  library(biomaRt)

  # Initialize biomaRt for human and mouse gene conversion
  ensembl_human <- useEnsembl("ensembl","hsapiens_gene_ensembl", mirror = "useast", host = "www.ensembl.org")
  ensembl_mouse <- useEnsembl("ensembl","mmusculus_gene_ensembl", mirror = "useast", host = "www.ensembl.org")
  
  # Check species and set up appropriate conversion
  if (species == "human_to_mouse") {
    # Map human gene symbols to mouse gene symbols
    result <- getLDS(
      attributes = c("hgnc_symbol"),           # Human gene symbols
      filters = "hgnc_symbol",                 # Use human gene symbols as input
      values = gene_list,                      # Input list of human genes
      mart = ensembl_human,                    # Human mart
      attributesL = c("mgi_symbol"),           # Get mouse gene symbols (MGI)
      martL = ensembl_mouse                    # Mouse mart
    )
    
    # Return only valid mouse genes (removing NA values)
    valid_mouse_genes <- result[!is.na(result$mgi_symbol), ]
    return(valid_mouse_genes$mgi_symbol)
    
  } else if (species == "mouse_to_human") {
    # Map mouse gene symbols to human gene symbols
    result <- getLDS(
      attributes = c("mgi_symbol"),            # Mouse gene symbols
      filters = "mgi_symbol",                  # Use mouse gene symbols as input
      values = gene_list,                      # Input list of mouse genes
      mart = ensembl_mouse,                    # Mouse mart
      attributesL = c("hgnc_symbol"),          # Get human gene symbols (HGNC)
      martL = ensembl_human                    # Human mart
    )
    
    # Return only valid human genes (removing NA values)
    valid_human_genes <- result[!is.na(result$hgnc_symbol), ]
    return(valid_human_genes$hgnc_symbol)
    
  } else {
    stop("Invalid species. Please specify 'human_to_mouse' or 'mouse_to_human'.")
  }
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


isString <- function(input) {
  is.character(input) & length(input) == 1
}

toCommaString <- function(genelist){
  # converting vector
  vec2 <- shQuote(genelist, type = "cmd")
  # combining elements using , 
  comma_vec2 <- paste(vec2, collapse = ", ")
  return(cat(comma_vec2))
}

df2longdf <- function(df){
  df$gene <- rownames(df)
  df <- df %>% gather(key='sample', value='value', -gene)
  return(df)
}


