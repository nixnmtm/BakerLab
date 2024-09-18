# Single Cell/Spatial Utility functions
# Author: Nixon Raj


####SPATA utils####  

convert2Seurat <- function (spata_obj, image_path){
  # convert a Spata object into seurat
  # using transformSpataToSeurat() doesnt add image to the seurat object
  # so created this function to manually add image and create a seurat obj
  
  seurat_obj <- SPATA2::transformSpataToSeurat(spata_obj, NormalizeData=F, 
                                               FindVariableFeatures=F, SCTransform = F,
                                               ScaleData=F, RunPCA=F, FindNeighbors=F, 
                                               FindClusters=F, RunTSNE = F,
                                               RunUMAP = F)
  seurat_obj <- base::tryCatch({
    lowres <- Read10X_Image(image_path, image.name = "tissue_lowres_image.png")
    seurat_obj@images$slice1 = lowres
    seurat_obj@images$slice1@assay = "Spatial"
    seurat_obj@images$slice1@key = "slice1_"
    base::warning("The SpatialImage is manually added to the Seurat Object")
    seurat_obj
  }, error = function(error) {
    base::warning("Error in adding image slice manually")
  })
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


####Seurat Visium####

plot_average_exp_HeatMap <- function(obj, cluster_colname){
  
  n = length(levels(bl[[cluster_colname]][[cluster_colname]]))
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

getNwriteDEG_df <- function(markers, path=NULL, file_name=NULL, pcut=1e-2, FCcut=1, 
                            rankbyFC=F,
                            rankbyPval=F,
                            rankbyAdjPval=F, 
                            write=F){
  #:markers: results of FindMarkers
  #:path: relative path to save
  #:filename: filename to save
  
  ge <- markers %>%
    as.data.frame() %>%
    filter(p_val < pcut & abs(avg_log2FC) > FCcut)
  
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
    write.table(ge, 
                file = file.path(merged_save_path, paste0(file_name, ".csv")),
                sep = ",", 
                row.names = F, 
                quote = FALSE)
  }
  return(ge)
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

mapgenes_mm_hs <- function(genes, mouse2human=F, human2mouse=F){
  # Converts genes from mouse 2 human or from human 2 mouse 
  library(Orthology.eg.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  if(mouse2human){
    gns <- mapIds(org.Mm.eg.db, genes, "ENTREZID", "SYMBOL")
    mapped <- select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
    naind <- is.na(mapped$Homo_sapiens)
    hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
    out <- data.frame(Mouse_symbol = genes, mapped, Human_symbol = NA)
    out$Human_symbol[!naind] <- hsymb
  }
  if(human2mouse){
    gns <- mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
    mapped <- select(Orthology.eg.db, gns, "Mus_musculus", "Homo_sapiens")
    naind <- is.na(mapped$Mus_musculus)
    msymb <- mapIds(org.Mm.eg.db, as.character(mapped$Mus_musculus[!naind]), "SYMBOL", "ENTREZID")
    out <- data.frame(Human_symbol = genes, mapped, Mouse_symbol = NA)
    out$Mouse_symbol[!naind] <- msymb
  }
  
  return(out)
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


check_save_dea_data <- function(markers, path, filename){
  
  if (file.exists(filename)) {
    #Delete file if it exists
    file.remove(filename)
    write.table(markers, file = file.path(path,filename),
                sep = ",",
                append = F,
                col.names=T, 
                row.names = F)
  }else{
    write.table(markers, file = file.path(path,filename),
                sep = ",",
                append = F,
                col.names=T, 
                row.names = F)
  }
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


vlnplot_w_significance <- function(obj, gene_signature, file_name, group_by = NONE, 
                                   test_sign, y_max,w=14, h=8){
  plot_case <- function(signature){
    library(ggplot2)
    library(ggpubr)
    
    VlnPlot(obj, features = signature,
            pt.size =  0.5, 
            group.by = group_by, 
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
    ) + 
      stat_compare_means(comparisons = test_sign, label = "p.signif") + 
      NoLegend()
    #+stat_compare_means(comparisons = test_sign, label = "p.format")
  }
  print(map(gene_signature, plot_case) %>% cowplot::plot_grid(plotlist = .))
  file_name <- paste0(file_name, "_VlnPlot.png")
  ggsave(file_name, width = w, height = h)
}


