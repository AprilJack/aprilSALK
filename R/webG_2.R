#' A GO Analysis Function 2.0
#'
#' This function takes in gene lists and runs ORA and optionally GSEA analysis on them, against user defined databases.
#' @param GSEA input is default NULL.
#' @keywords GO, GSEA, ORA, WebGestalt
#' @export
#' @examples
#' webG_2()

webG_2 <- function(counts,resSig,org,db,name,gsea.df=NULL){
  #Use WebG
  #sapply(db,webgX <- function(x) {
  #  webG(counts,resSig,org,x,name)})
  #C: the number of reference genes in the category
  #O: the number of genes in the user gene list and also in the category
  #E: The expected number in the category
  #R: ratio of enrichment
  #library(ggplot2)
  #library("WebGestaltR")#,lib.loc = "/gpfs/analyses/april/April/lib/R/library/")
  print("#########################################################")
  print(sprintf("Starting GO Analysis for - %s using %s", name, db))
  print("#########################################################")
  ######### EdgeR Section #########
  if(!is.null(resSig$logFC)){
    if(dim(resSig[resSig$logFC > 0,])[1] > 10){
      print("...UP ORA EdgeR...")
      Up.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"UP","edgeR"),error=function(e) return(NULL) )
      if (!is.null(dim(Up.ORA)[1])) {
        print("...UP ORA EdgeR Graph...")
        gORA(Up.ORA,db,name,"UP","edgeR")
      }
    }
    if(dim(resSig[resSig$logFC < 0,])[1] > 10){
      print("...DOWN ORA EdgeR...")
      Down.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"DN","edgeR"),error=function(e) return(NULL) )
      if (!is.null(dim(Down.ORA)[1])) {
        print("...DOWN ORA EdgeR Graph...")
        gORA(Down.ORA,db,name,"DN","edgeR")
      }
    }
    if(dim(resSig)[1] > 10){
      print("...Combined ORA EdgeR...")
      Combined.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"Combined","edgeR"),error=function(e) return(NULL) )
      if (!is.null(dim(Combined.ORA)[1])) {
        print("...Combined ORA EdgeR Graph...")
        tryCatch(gORA(Combined.ORA,db,name,"Combined","edgeR"),error=function(e) return(NULL) )
      }
    }
    if(!is.null(gsea.df)){
    if(dim(gsea.df)[1] > 10){
      print("...GSEA EdgeR...")
      GSEA <- tryCatch(rGSEA(gsea.df,org,db,name,"edgeR"),error=function(e) return(NULL) )
      if (!is.null(GSEA[GSEA$FDR < 0.05,]$normalizedEnrichmentScore)) {
        print("...GSEA EdgeR Graph...")
        gGSEA(GSEA,db,name,"edgeR")
      }
    }
    }
  }else if(!is.null(resSig$meth.diff)){
    ######## MethylSeq section ########
    if(length(unique(resSig[resSig$meth.diff > 25,]$Nearest.Refseq))[1] > 10){
      print("...Hyper-methylated ORA MethylSeq...")
      Hyper.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"Hyper","MS"),error=function(e) return(NULL) )
      if (!is.null(dim(Hyper.ORA)[1])) {
        print("...Hyper ORA MethylSeq Graph...")
        gORA(Hyper.ORA,db,name,"Hyper","MS")
      }
    }
    if(length(unique(resSig[resSig$meth.diff < -25,]$Nearest.Refseq))[1] > 10){
      print("...Hypo-methylated ORA MethylSeq...")
      Hypo.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"Hypo","MS"),error=function(e) return(NULL) )
      if (!is.null(dim(Hypo.ORA)[1])) {
        print("...Hyper ORA MethylSeq Graph...")
        gORA(Hypo.ORA,db,name,"Hypo","MS")
      }
    }
  } else if(!is.null(resSig$log2FoldChange)){
    ##### DESeq2 section #######
    names(resSig)[names(resSig) == "log2FoldChange"] <- "logFC"
    if(dim(resSig[resSig$logFC > 0,])[1] > 10){
      print("...UP ORA DESeq2...")
      Up.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"UP","DESeq2"),error=function(e) return(NULL) )
      if (!is.null(dim(Up.ORA)[1])) {
        print("...UP ORA DESeq2 Graph...")
        gORA(Up.ORA,db,name,"UP","DESeq2")
      }
    }
    if(dim(resSig[resSig$logFC < 0,])[1] > 10){
      print("...DOWN ORA DESeq2...")
      Down.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"DN","DESeq2"),error=function(e) return(NULL) )
      if (!is.null(dim(Down.ORA)[1])) {
        print("...DOWN ORA DESeq2 Graph...")
        gORA(Down.ORA,db,name,"DN","DESeq2")
      }
    }
    if(dim(resSig)[1] > 10){
      print("...Combined ORA DESeq2...")
      Combined.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"Combined","DESeq2"),error=function(e) return(NULL) )
      if (!is.null(dim(Combined.ORA)[1])) {
        gORA(Combined.ORA,db,name,"Combined","DESeq2")
      }
    }
    if(!is.null(gsea.df)){
    if(dim(gsea.df)[1] > 10){
      print("...GSEA DESeq2 w Threads...")
      GSEA <- tryCatch(rGSEA(gsea.df,org,db,name,"DESeq2"),error=function(e) return(NULL) )
      if (!is.null(GSEA[GSEA$FDR < 0.05,]$normalizedEnrichmentScore)) {
        print("...GSEA DESeq2 Graph...")
        gGSEA(GSEA,db,name,"DESeq2")
      }
    }
    }
  }else if(!is.null(resSig$avg_logFC)){
    ##### Seurat3 section #######
    names(resSig)[names(resSig) == "avg_logFC"] <- "logFC"
    if(dim(resSig[resSig$logFC > 0,])[1] > 10){
      print("...UP ORA Seurat3...")
      Up.ORA <- tryCatch(rORA(gsea.df,resSig,org,db,name,"UP","Seurat"),error=function(e) return(NULL) )
      if (!is.null(dim(Up.ORA)[1])) {
        print("...UP ORA Seurat3 Graph...")
        gORA(Up.ORA,db,name,"UP","Seurat")
      }
    }
    if(dim(resSig[resSig$logFC < 0,])[1] > 10){
      print("...DOWN ORA Seurat3...")
      Down.ORA <- tryCatch(rORA(gsea.df,resSig,org,db,name,"DN","Seurat"),error=function(e) return(NULL) )
      if (!is.null(dim(Down.ORA)[1])) {
        print("...DOWN ORA Seurat3 Graph...")
        gORA(Down.ORA,db,name,"DN","Seurat")
      }
    }
    if(dim(resSig)[1] > 10){
      print("...Combined ORA Seurat3...")
      Combined.ORA <- tryCatch(rORA(gsea.df,resSig,org,db,name,"Combined","Seurat"),error=function(e) return(NULL) )
      if (!is.null(dim(Combined.ORA)[1])) {
        print("...Combined ORA Seurat3 Graph...")
        gORA(Combined.ORA,db,name,"Combined","Seurat")
      }
    }
    if(!is.null(gsea.df)){
    if(dim(gsea.df)[1] > 10){
      print("...GSEA Seurat3...")
      ddd <- data.frame(gene = row.names(gsea.df[order(gsea.df$avg_logFC),]),avg_logFC=gsea.df[order(gsea.df$avg_logFC),]$avg_logFC)
      ddd[ddd== -Inf] <- -1000
      ddd[ddd== Inf] <- 1000
      GSEA <- tryCatch(rGSEA(ddd,org,db,name,"Seurat"),error=function(e) return(NULL) )
      if (!is.null(GSEA[GSEA$FDR < 0.05,]$normalizedEnrichmentScore)) {
        print("...GSEA Seurat3 Graph...")
        gGSEA(GSEA,db,name,"Seurat")
      }
    }
    }
  }else if(!is.null(resSig$adj.P.Val)){
    if(dim(resSig[resSig$adj.P.Val > 0,])[1] > 10){
      print("...UP ORA SplineTimeR...")
      Up.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"UP","SplineTimeR"),error=function(e) return(NULL) )
      if (!is.null(dim(Up.ORA)[1])) {
        print("...UP ORA SplineTimeR Graph...")
        gORA(Up.ORA,db,name,"UP","SplineTimeR")
      }
    }
    if(dim(resSig[resSig$adj.P.Val < 0,])[1] > 10){
      print("...DOWN ORA SplineTimeR...")
      Down.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"DN","SplineTimeR"),error=function(e) return(NULL) )
      if (!is.null(dim(Down.ORA)[1])) {
        print("...DOWN ORA SplineTimeR Graph...")
        gORA(Down.ORA,db,name,"DN","SplineTimeR")
      }
    }
    if(dim(resSig)[1] > 10){
      print("...Combined ORA SplineTimeR...")
      Combined.ORA <- tryCatch(rORA(counts,resSig,org,db,name,"Combined","SplineTimeR"),error=function(e) return(NULL) )
      if (!is.null(dim(Combined.ORA)[1])) {
        print("...Combined ORA SplineTimeR Graph...")
        tryCatch(gORA(Combined.ORA,db,name,"Combined","SplineTimeR"),error=function(e) return(NULL) )
      }
    }
    if(!is.null(gsea.df)){
      if(dim(gsea.df)[1] > 10){
        print("...GSEA SplineTimeR...")
        GSEA <- tryCatch(rGSEA(gsea.df,org,db,name,"SplineTimeR"),error=function(e) return(NULL) )
        if (!is.null(GSEA[GSEA$FDR < 0.05,]$normalizedEnrichmentScore)) {
          print("...GSEA SplineTimeR Graph...")
          gGSEA(GSEA,db,name,"SplineTimeR")
        }
      }
    }
  } else{
    print("Please enter a valid DESeq2, EdgeR, Seurat, or MethylSeq object.") }
}




