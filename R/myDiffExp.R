#' A Differential Expression Function
#'
#' This function allows you explore basic QC of RNA-seq data and run differential expression testing and GO analysis.
#' @param no defaults.
#' @keywords QC, explore, DESeq2, edgeR
#' @export
#' @examples
#' myDiffExp()


myDiffExp <- function(ddsIn=NULL,counts=NULL,design=NULL,colData1,colData2,formulaEdgR,formulaDE2,
                      ref,name,org,anno,db,runWebG=TRUE,fdr=0.05,lfc=c(1,-1),printToScreen=F) {
  # library("edgeR")
  # library("DESeq2")
  # library("stringr")
  # library("gplots")
  # library("ggplot2")
  # library("pheatmap")
  # library("cluster")
  # library("VennDiagram")
  # library("WebGestaltR")
  if(printToScreen == F){
  pdf(sprintf("%s_plots.pdf",name),width = 7,height = 7)}

  print("#########################################################")
  print("Starting DESeq2 General QC")
  print("#########################################################")

  if(!(is.null(ddsIn))){
    design <- as.data.frame(colData(ddsIn))
    counts <- counts(ddsIn,normalized=F)
    rm(ddsIn)
  }
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

  #cRef <- counts[,row.names(design[design[[colData1]] == ref,])]
  #cOth <- counts[,row.names(design[design[[colData1]] != ref,])]

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
  saveRDS(object = dds,file = sprintf("%s_dds.RDS",name))
  if (!is.null(formulaEdgR)) {

    print("#########################################################")
    print("Starting differential expression analysis using EdgeR.")
    print("#########################################################")
    #set base condition design for making designEdgeR
    design[,colData1] <- factor(design[,colData1])
    design[,colData1] <- relevel(design[,colData1], ref=ref)

    #create edgeR object
    y <- DGEList(counts=counts[keep,]*2,group=design[,colData1])

    #set group levels fo edgeR object
    y$samples$group <- relevel(y$samples$group, ref=ref)

    #calculate norm factors using TMM method using edgeR
    y <- calcNormFactors(y,method="TMM")

    #output normalized count table to edegeRExp
    edegeRExp <- cpm(y)

    if(stringr::str_starts(row.names(as.data.frame(edegeRExp))[1],"NM") == TRUE | stringr::str_starts(row.names(as.data.frame(edegeRExp))[1],"NR") == TRUE ){
      edegeRExp <- merge(x = anno, y = as.data.frame(edegeRExp), by.x= 0, by.y=0)
    }else{edegeRExp <- merge(x = anno, y = as.data.frame(edegeRExp), by.x= "gene_name", by.y=0)}

    openxlsx::write.xlsx(x = edegeRExp,file = sprintf("%s_edgeR_norm_counts.xlsx",name),firstRow=T,keepNA=T)

    designEdgeR <- model.matrix(as.formula(formulaEdgR))
    y <- estimateDisp(y,designEdgeR)
    fit <- glmQLFit(y,designEdgeR)
    qlf <- glmQLFTest(fit)
    topTags(qlf)
    res <- topTags(qlf,n=Inf)
    print("head of results")
    print(res)
    print("Differentially expressed genes:")
    print(table(p.adjust(qlf$table$PValue, method = "BH", n = length(qlf$table$PValue)) < fdr & abs(qlf$table$logFC) > lfc[1]))

    print("...graphing MA plot of differential expression results...")
    qlf$table$clr <- ifelse((p.adjust(qlf$table$PValue, method = "BH", n = length(qlf$table$PValue)) < fdr
                             & abs(qlf$table$logFC) > lfc[1]), "red", "black")
    plotSmear(object =qlf,col=qlf$table$clr,main="MA Plot - Significant Genes Highlighted")#,ylim=c(-2,2))

    if(stringr::str_starts(row.names(as.data.frame(res))[1],"NM") == TRUE | stringr::str_starts(row.names(as.data.frame(res))[1],"NR") == TRUE){
      results <- merge(x = anno, y = as.data.frame(res), by.x= 0, by.y=0)
    }else{print("here");results <- merge(x = anno, y = as.data.frame(res), by.x= "gene_name", by.y=0)}

    write.table(results,file = sprintf("%s.csv",name),sep=",",quote = F,row.names = F, col.names = T)
    openxlsx::write.xlsx(x = results,file = sprintf("%s.xlsx",name),firstRow=T,keepNA=T)


    resSig <- subset(res[[1]],FDR < fdr & (logFC > lfc[1] | logFC < lfc[2]))
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

    if(printToScreen==F){dev.off()}

    resSigUp <- subset(res[[1]],FDR < fdr & logFC > lfc[1] )
    resSigDown <- subset(res[[1]],FDR < fdr & logFC < lfc[2])

    if(runWebG==TRUE){
    if(dim(resSigUp)[1] > 10 |  dim(resSigDown)[1] > 10){
    sapply(db,webG_X <- function(x) {
      webG_2(counts,resSig,org,x,name,gsea.df)})
    }else {print("Not enough significant results to run GO term analysis.")}}else{print("Skipping WebG.")}

    saveRDS(object = fit,file = sprintf("%s_edgeR.RDS",name))
  } else {
    print("#########################################################")
    print("Starting differential expression analysis using DESeq2")
    print("#########################################################")

    resultsNames(dds)
    print(resultsNames(dds))
    res <- results(dds,alpha=fdr, lfcThreshold=log2(2))
    summary(res)#, alpha=fdr)

    print("...printing differential expression results...")
    if(stringr::str_starts(row.names(as.data.frame(res))[1],"NM") == TRUE | stringr::str_starts(row.names(as.data.frame(res))[1],"NR") == TRUE){
      results <- merge(x = anno, y = as.data.frame(res), by.x= 0, by.y=0)
    }else{print("gene_name");results <- merge(x = anno, y = as.data.frame(res), by.x= "gene_name", by.y=0)}
    print(head(results))
    write.table(results,file = sprintf("%s.csv",name),sep=",",quote = F,row.names = F, col.names = T)
    openxlsx::write.xlsx(x = results,file = sprintf("%s.xlsx",name),firstRow=T,keepNA=T)

    print("...graphing MA plot of differential expression results...")
    DESeq2::plotMA(res,alpha = fdr)#, ylim=c(-2,2))

    resSig <- subset(res,padj < fdr & (log2FoldChange > lfc[1] | log2FoldChange < lfc[2]))
    resSigUp <- subset(res,padj < fdr & log2FoldChange > lfc[1])
    resSigDown <- subset(res,padj < fdr & log2FoldChange < lfc[2])

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

    if(printToScreen==F){dev.off()}

    if(runWebG==TRUE){
    if(dim(resSigUp)[1] > 10 |  dim(resSigDown)[1] > 10){
      sapply(db,webG_X <- function(x) {
        webG_2(counts,resSig,org,x,name,gsea.df)})
    } else {print("Not enough significant results to run GO term analysis.")}}else{print("Skipping WebG.")}

      saveRDS(object = dds,file = sprintf("%s_dds.RDS",name))
  }
  if(printToScreen==F){graphics.off()}
  return(results)
}
