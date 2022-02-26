#' A Bisulfite Sequencing Function
#'
#' This function takes in a raw methylKit object and does basic QC, then tests for differential methylation in 25bp windows with 15 bp step-size.
#' @param no defaults
#' @keywords cats
#' @export
#' @examples
#' meth25bp()


meth25bp <- function(methylseqObj = "methy seq obj",nam = "name of exp",gene.obj = gene.obj,org="organism name for WebGestaltR",db="GO db's to query"){
  # library("methylKit")
  # library("genomation")
  # library("ggplot2")
  # library("gplots")
  # library("pheatmap")
  # library("cluster")
  pdf(sprintf("%s_plots.pdf",nam),width = 7,height = 7)
  print("---------Filtering methyl object")
  filtered.methylseqObj = filterByCoverage(methylseqObj,lo.count=10,lo.perc=NULL,
                                           hi.count=NULL,hi.perc=99.9)

  print("---------Normalizing methyl object")
  normed.methylseqObj <- normalizeCoverage(filtered.methylseqObj)

  print("---------Tiling with 25bp windows with 15bp steps")
  tiles <- tileMethylCounts(normed.methylseqObj,win.size=25,step.size=15)

  print("---------Uniting strands of methyl object")
  meth <- unite(tiles, destrand=FALSE)

  print("---------Calculating samples correlations")
  getCorrelation(meth,plot=T)
  print("---------Clustering samples by correlation")
  clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
  print("---------PCA analysis")
  PCASamples(meth, screeplot=TRUE)
  PCASamples(meth)
  print("---------Calculating methylation percentages")
  percmeth <- percMethylation(meth)

  print("---------Calculating differential methylation")
  myDiff <- calculateDiffMeth(meth,mc.cores=4)
  diffMethPerChr(myDiff,plot=T,qvalue.cutoff=0.01, meth.cutoff=25,exclude=c("chrM","chrX","chrY"))

  print("---------Annotating hyper and hypo methylated using Genomation")
  myDiff.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
  diffAnnHyper <- annotateWithGeneParts(as(myDiff.hyper,"GRanges"),gene.obj)
  genomation::plotTargetAnnotation(diffAnnHyper,precedence=TRUE,main=sprintf("%s\ndifferential methylation\nannotation hyper",nam))

  myDiff.hypo <- getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
  diffAnnHypo <- annotateWithGeneParts(as(myDiff.hypo,"GRanges"),gene.obj)
  genomation::plotTargetAnnotation(diffAnnHypo,precedence=TRUE,main=sprintf("%s\ndifferential methylation\nannotation hypo",nam))

  print("---------Writing out Genomation annotation results")
  diffAnnAll <- annotateWithGeneParts(as(myDiff,"GRanges"),gene.obj)
  annotatedmeth <- data.frame(myDiff,percmeth,diffAnnAll@members)
  write.table(x = annotatedmeth,file = sprintf("%s_annotatedDiffMeth.txt", nam),sep = "\t",quote = F,row.names = F)
  openxlsx::write.xlsx(x = annotatedmeth,file = sprintf("%s_annotatedDiffMeth.xlsx", nam),firstRow=T,keepNA=T)


  print("---------Annotating hyper and hypo methylated regions using HOMER")
  #annotate methyl regions using HOMER
  input <- sprintf("%s_annotatedDiffMeth.txt", nam)
  output <- sprintf("%s_annotatedDiffMeth_homerAnn.txt", nam)
  system(sprintf("source /gpfs/tools/.bashTools;
     annotatePeaks.pl %s mm10 > %s", input, output))
  annotatedmeth$PeakID <- paste("*-",row.names(annotatedmeth),sep="")
  annotatedHomermeth <- read.delim(output, header=T, stringsAsFactors=F, sep="\t")
  Final_annotatedHomermeth <- merge(annotatedmeth,annotatedHomermeth[,-c(2:4,6:7)],by.x="PeakID",by.y=1)
  write.table(x = Final_annotatedHomermeth,file = sprintf("%s_Final_annotatedDiffMeth.txt",nam),sep = "\t",quote = F,row.names = F)
  openxlsx::write.xlsx(x = Final_annotatedHomermeth,file = sprintf("%s_Final_annotatedDiffMeth.xlsx", nam),firstRow=T,keepNA=T)

  diffSigAnn <- subset(Final_annotatedHomermeth, qvalue < 0.01 & abs(meth.diff) > 25)
  getGenesNearCoordinateMouse

  sapply(db,webgX <- function(x) {
    webG_2(Final_annotatedHomermeth,diffSigAnn,org,x,nam)})

  graphics.off()

  Final_annotatedHomermeth50kb <- getGenesNearCoordinateMouse(Final_annotatedHomermeth,regionSize = 50000)
  write.table(x = Final_annotatedHomermeth50kb,file = sprintf("%s_Final_annotated_50kb.txt",nam),quote = FALSE,append = FALSE,sep = "\t",row.names = F)
  openxlsx::write.xlsx(x = Final_annotatedHomermeth50kb,file = sprintf("%s_Final_annotated_50kb.xlsx", nam),firstRow=T,keepNA=T)

  return(Final_annotatedHomermeth50kb)
}
