#' A Data Exploration Function
#'
#' This function allows you explore basic QC of RNA-seq data without running differential expression testing or GO analysis.
#' @param no defaults.
#' @keywords QC, explore, DESeq2, edgeR
#' @export
#' @examples
#' myExplore()


myExplore <- function(counts,design,colData1,colData2,formulaDE2,ref,name,anno) {
  suppressPackageStartupMessages(library("DESeq2"))
  suppressPackageStartupMessages(library("stringr"))
  suppressPackageStartupMessages(library("gplots"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("pheatmap"))
  suppressPackageStartupMessages(library("cluster"))
  suppressPackageStartupMessages(library("VennDiagram"))
  pdf(sprintf("%s_plots.pdf",name),width = 7,height = 7)

  print("#########################################################")
  print("Starting DESeq2 General QC")
  print("#########################################################")

  dds <- DESeqDataSetFromMatrix(countData = counts*2,
                                colData = design,
                                design = as.formula(formulaDE2))
  dds[[colData1]] <- relevel(dds[[colData1]], ref = ref)

  summary(rowSums(counts(dds)))
  hist(log10(rowSums(counts(dds))))
  abline(v=log10(10),col="red")

  keep <- rowSums(counts(dds)) >= 10
  print("Genes passing expression filter:")
  print(table(keep))


  #cRef <- counts[,row.names(design[design[[colData1]] == ref,])]
  #cOth <- counts[,row.names(design[design[[colData1]] != ref,])]

  #tryCatch(makeVenns(overlap=list(row.names(cRef[rowSums(cRef) >= 10,]),row.names(cOth[rowSums(cOth) >= 10,])),
  #          nam=name,
  #          cats = c(levels(droplevels(as.factor(design[[colData1]])))[1], levels(droplevels(as.factor(design[[colData1]])))[2]),
  #          geneAnno = anno,
  #         byy = 0),error=function(e) return(NULL) )

  dds <- dds[keep,]
  dds <- suppressMessages(DESeq(dds))

  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

  pca <- prcomp(t(assay(vsd)))
  summaryPCA <- summary(pca)
  print("PCA summary:")
  print(summaryPCA)

  loadings = as.data.frame(pca$x[])

  print("#########################################################")
  print("Starting QC graphing of VSD transformed raw count data.")
  print("#########################################################")

  print("...graphing PCA...")
  pca12 <- ggplot(data = loadings,
                  aes(x=PC1,
                      y=PC2,
                      color=as.factor(design[[colData1]]),
                      shape=as.factor(design[[colData2]])))
  print(pca12 +
          geom_point(size=17) +
          scale_shape_manual(values=seq(21, 21+length(levels(as.factor(design[[colData2]])))-1, 1)) +
          theme_bw() +
          labs(x=sprintf("PC1 (%s of variance)",paste(round(summaryPCA$importance[2,1]*100, 2), "%", sep="")),
               y=sprintf("PC2 (%s of variance)",paste(round(summaryPCA$importance[2,2]*100, 2), "%", sep=""))) +
          geom_text(aes(label=row.names(design))) +
          ggtitle(sprintf("PCA of normalized raw count data - %s",name)) + labs(color = colData1,shape=colData2))

  print("...graphing Euclidean distance matrix...")
  sampleDists <- dist(t(assay(vsd)))
  df <- as.data.frame(colData(vsd)[,c(colData1,colData2)])
  sampleDistMatrix <- as.matrix(sampleDists)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           show_rownames=FALSE,
           main="Euclidean distance of normalized counts",
           annotation_col=df,annotation_row=df
  )

  print("...graphing heatmap of top 100 expressed genes...")
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:100]

  df <- as.data.frame(colData(vsd)[,c(colData1,colData2)])
  pheatmap(assay(vsd)[select,],
           cluster_rows=TRUE,
           show_rownames=TRUE,
           cluster_cols=TRUE,
           annotation_col=df,
           main=sprintf("%s - Heatmap of Top 100 Expressed Genes",name),
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "correlation",
           clustering_method = "ward.D2",
           treeheight_row=25,
           fontsize = 6,
           color = colorRampPalette(c("green", "black", "red"))(50))

  dev.off()

  return(vsd)
}
