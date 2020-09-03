#' First steps of scRNA QC and filtering
#'
#' This function allows you explore basic QC of scRNA-seq data.
#' @param no defaults.
#' @keywords QC, scRNA, explore, Seurat
#' @export
#' @examples
#' seuratQC()


seuratQC <- function(ids, nam, pmito, nf_max, nf_min){
  library(Seurat)
  library(ggplot2)
  pdf(file = sprintf("%s_plots.pdf",nam),width = 11,height = 8.5)

  d10x.data <- sapply(ids, function(i){
    d10x <- Read10X(file.path(i,"/outs/filtered_feature_bc_matrix/"))
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
    d10x
  })

  experiment.data <- do.call("cbind", d10x.data)

  seurat.obj <- CreateSeuratObject(
    experiment.data,
    project = nam,
    min.cells = 10,
    min.features = 200,
    names.field = 3,
    names.delim = "\\-")

  mito.genes.190708 <- grep("^mt-", rownames(seurat.obj), value = T)
  percent.mito <- Matrix::colSums(GetAssayData(seurat.obj,slot="counts")[mito.genes.190708,]) /
    Matrix::colSums(GetAssayData(seurat.obj,slot="counts"))
  seurat.obj$percent.mito <- percent.mito

  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  # Visualize QC metrics as a violin plot
  print(VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,group.by = "orig.ident"))
  vln <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,group.by = "orig.ident")
  ggsave(filename=sprintf("%s_vln_plot.pdf",nam), vln)

  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

  pM <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mito",group.by = "orig.ident")
  pM <- pM + geom_hline(yintercept=pmito,linetype="dashed")

  nFRNA <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")
  nFRNA <- nFRNA + geom_hline(yintercept=nf_min,linetype="dashed") + geom_hline(yintercept=nf_max,linetype="dashed")

  print(CombinePlots(plots = list(pM, nFRNA)))
  cp <- CombinePlots(plots = list(pM, nFRNA))
  ggsave(filename=sprintf("%s_scatter_plots.pdf",nam), cp)


  #filter data based on features too few and too many, and percent mitochondrial reads. Set at conservative 10%.
  print("...subsetting...")
  nf_minx <<- nf_min
  nf_maxx <<- nf_max
  pmitox <<- pmito
  seurat.obj <- subset(seurat.obj, nFeature_RNA > nf_minx & nFeature_RNA < nf_maxx & percent.mito < pmitox)

  print("here1")
  table(seurat.obj$orig.ident)

  print("here2")
  seurat.obj <- NormalizeData(object = seurat.obj,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)

  seurat.obj <- ScaleData(seurat.obj)

  seurat.obj<-FindVariableFeatures(seurat.obj)

  #plot top 10 variable genes
  top10 <- head(VariableFeatures(seurat.obj), 10)
  plot1 <- VariableFeaturePlot(seurat.obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot2)
  ggsave(filename=sprintf("%s_variable_genes.pdf",nam), plot2)

  print("run PCA")
  seurat.obj <- RunPCA(seurat.obj,do.print=F)
  print("run PCA HM")
  print(PCHeatmap(seurat.obj))#,pc.use = 1,do.balanced = FALSE))
  pch <- PCHeatmap(seurat.obj)#,pc.use = 1,do.balanced = FALSE)
  ggsave(filename=sprintf("%s_pca_heatmap.pdf",nam), pch)

  print("run PCA DIM")
  print(DimPlot(seurat.obj,reduction = "pca"))
  pca <- DimPlot(seurat.obj,reduction = "pca")
  ggsave(filename=sprintf("%s_pca_plot.pdf",nam), pca)

  #build tree
  print("run buildTree")
  seurat.obj <- BuildClusterTree(seurat.obj)
  print(PlotClusterTree(seurat.obj))
  ct <- PlotClusterTree(seurat.obj)
  ggsave(filename=sprintf("%s_cluster_tree_plot.pdf",nam), ct)

  #Clustering
  seurat.obj <- FindNeighbors(seurat.obj)
  seurat.obj <- FindClusters(seurat.obj, resolution = .25)
  #dim reduction using TSNE
  seurat.obj <- RunTSNE(seurat.obj,dims.use = 1:11,max_iter=2000)

  #tsne by cluster
  print(DimPlot(object = seurat.obj, reduction = "tsne"))
  tsne <- DimPlot(object = seurat.obj, reduction = "tsne")
  ggsave(filename=sprintf("%s_tsne_cluster.pdf",nam), tsne)

  #tsne by original identity
  print(DimPlot(object = seurat.obj, reduction = "tsne",group.by = "orig.ident"))
  tsne.orig <- DimPlot(object = seurat.obj, reduction = "tsne",group.by = "orig.ident")
  ggsave(filename=sprintf("%s_tsne_orig.pdf",nam), tsne.orig)

  #tsne split by cluster
  print(DimPlot(object = seurat.obj, reduction = "tsne", split.by = "ident"))
  tsne.split <- DimPlot(object = seurat.obj, reduction = "tsne", split.by = "ident")
  ggsave(filename=sprintf("%s_tsne_split.pdf",nam), tsne.split)

  #UMAP
  #to use umap you need to activate python27 environment for reticulate to find umap pythong libraries.
  library(reticulate)
  use_condaenv(condaenv="python27", conda="/gpfs/tools/anaconda/bin/conda")

  seurat.obj <- RunUMAP(seurat.obj,dims=1:25)
  print(DimPlot(object = seurat.obj, reduction = "umap"))
  umap.clust <- DimPlot(object = seurat.obj, reduction = "umap")
  ggsave(filename=sprintf("%s_umap_clust.pdf",nam), umap.clust)

  print(DimPlot(object = seurat.obj, reduction = "umap",group.by = "orig.ident"))
  umap.orig <- DimPlot(object = seurat.obj, reduction = "umap",group.by = "orig.ident")
  ggsave(filename=sprintf("%s_umap_orig.pdf",nam), umap.orig)

  dev.off()
  return(seurat.obj)
}
