#' A Differential Expression Function
#'
#' This function allows you explore basic QC of RNA-seq data and run differential expression testing and GO analysis.
#' @param no defaults.
#' @keywords QC, explore, DESeq2, edgeR
#' @export
#' @examples
#' myDiffExp()


myDiffExp <- function(counts,design,colData1,colData2,formulaEdgR,formulaDE2,ref,name,org,anno,db) {
  library("edgeR")
  library("DESeq2")
  #library("stringr")
  library("gplots")
  library("ggplot2")
  library("pheatmap")
  library("cluster")
  library("VennDiagram")
  library("WebGestaltR")
  pdf(sprintf("%s_plots.pdf",name),width = 7,height = 7)

  print("#########################################################")
  print("Starting DESeq2 General QC")
  print("#########################################################")

  design <- droplevels(design)

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

  cRef <- counts[,row.names(design[design[[colData1]] == ref,])]
  cOth <- counts[,row.names(design[design[[colData1]] != ref,])]

  #  makeVenns(overlap=list(row.names(cRef[rowSums(cRef) >= 10,]),row.names(cOth[rowSums(cOth) >= 10,])),
  #            nam=name,
  #            cats = c(levels(droplevels(design[[colData1]]))[1], levels(droplevels(design[[colData1]]))[2]),
  #            geneAnno = anno,
  #            byy = 0)

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

  if (!is.null(formulaEdgR)) {

    print("#########################################################")
    print("Starting differential expression analysis using EdgeR.")
    print("#########################################################")

    y <- DGEList(counts=counts[keep,]*2,group=design[,colData1])
    y$samples$group <- relevel(y$samples$group, ref=ref)
    y <- calcNormFactors(y)
    designEdgeR <- model.matrix(as.formula(formulaEdgR))
    y <- estimateDisp(y,designEdgeR)
    fit <- glmQLFit(y,designEdgeR)
    qlf <- glmQLFTest(fit,coef=2)#coef=1#, lfc=1)
    topTags(qlf)
    res <- topTags(qlf,n=Inf)
    print("Differentially expressed genes:")
    print(table(p.adjust(qlf$table$PValue, method = "BH", n = length(qlf$table$PValue)) < 0.05 & abs(qlf$table$logFC) > 1))

    print("...graphing MA plot of differential expression results...")
    qlf$table$clr <- ifelse((p.adjust(qlf$table$PValue, method = "BH", n = length(qlf$table$PValue)) < 0.05
                             & abs(qlf$table$logFC) > 1), "red", "black")
    plotSmear(object =qlf,col=qlf$table$clr,main="MA Plot - Significant Genes Highlighted")#,ylim=c(-2,2))

    results <- merge(x = anno, y = as.data.frame(res), by.x= 0, by.y=0)
    write.table(results,file = sprintf("%s.csv",name),sep=",",quote = F,row.names = F, col.names = T)
    openxlsx::write.xlsx(x = results,file = sprintf("%s.xlsx",name),firstRow=T,keepNA=T)


    resSig <- subset(res[[1]],FDR < 0.05 & (logFC > 1 | logFC < -1))
    #resSigUp <- subset(res[[1]],FDR < 0.05 & logFC > 1)
    #resSigDown <- subset(res[[1]],FDR < 0.05 & logFC < -1)
    logcpm <- cpm(y, prior.count=2, log=TRUE)

    gsea.df <- data.frame(gene_names = row.names(res),logFC = res$table$logFC)
    gsea.df <- gsea.df[order(gsea.df$logFC),]

    if(dim(resSig)[1] > 1){
      df <- as.data.frame(colData(vsd)[,c(colData1,colData2)])
      pheatmap(logcpm[row.names(resSig),],
               cluster_rows=TRUE,
               show_rownames=F,
               cluster_cols=TRUE,
               annotation_col=df,
               main=sprintf("%s - Heatmap of Differentially Expressed Genes",name),
               scale = "row",
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "correlation",
               clustering_method = "ward.D2",
               treeheight_row=25,
               fontsize = 8,
               color = colorRampPalette(c("green", "black", "red"))(50))
    }

    dev.off()

    sapply(db,webG_X <- function(x) {
      webG_2(counts,resSig,org,x,name,gsea.df)})
  } else {
    print("#########################################################")
    print("Starting differential expression analysis using DESeq2")
    print("#########################################################")

    resultsNames(dds)
    print(resultsNames(dds))
    res <- results(dds,alpha=0.05, lfcThreshold=log2(2))
    summary(res)#, alpha=0.05)

    print("...printing differential expression results...")
    results <- merge(x = anno, y = as.data.frame(res), by.x= 0, by.y=0)
    print(head(results))
    write.table(results,file = sprintf("%s.csv",name),sep=",",quote = F,row.names = F, col.names = T)
    openxlsx::write.xlsx(x = results,file = sprintf("%s.xlsx",name),firstRow=T,keepNA=T)

    print("...graphing MA plot of differential expression results...")
    DESeq2::plotMA(res,alpha = 0.05)#, ylim=c(-2,2))

    resSig <- subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
    resSigUp <- subset(res,padj < 0.05 & log2FoldChange > 1)
    resSigDown <- subset(res,padj < 0.05 & log2FoldChange < -1)

    gsea.df <- data.frame(gene_names = row.names(res),log2FoldChange = res$log2FoldChange)
    gsea.df <- gsea.df[order(gsea.df$log2FoldChange),]

    if(dim(resSig)[1] > 1){
      df <- as.data.frame(colData(vsd)[,c(colData1,colData2)])
      pheatmap(assay(vsd)[row.names(resSig),],
               cluster_rows=TRUE,
               show_rownames=F,
               cluster_cols=TRUE,
               annotation_col=df,
               main=sprintf("%s - Heatmap of Differentially Expressed Genes",name),
               scale = "row",
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "correlation",
               clustering_method = "ward.D2",
               treeheight_row=25,
               fontsize = 8,
               color = colorRampPalette(c("green", "black", "red"))(50))
    }

    dev.off()

    sapply(db,webG_X <- function(x) {
      webG_2(counts,resSig,org,x,name,gsea.df)})

  }
  graphics.off()
  return(results)
}
