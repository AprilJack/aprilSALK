#' An Annotation Function
#'
#' This function finds all genes surrounding methylated bases +/- a variable region size in mouse.
#' @param regionSize How far +/- do you want to search for genes. Defaults to 50,000 bp.
#' @keywords mouse, annotation
#' @export
#' @examples
#' getGenesNearCoordinateMouse()


getGenesNearCoordinateMouse <- function(annMethResult, regionSize = 50000){
  if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
    stop("Package \"TxDb.Mmusculus.UCSC.mm10.knownGene\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop("Package \"org.Mm.eg.db\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(dplyr)
  library(GenomicRanges)
  library(org.Mm.eg.db)

  #create a data frame of gene ids and symbols
  eid <- as.data.frame(org.Mm.egSYMBOL)
  #creat a coordinate identifier
  annMethResult$coord <- do.call(paste, c(annMethResult[,2:4],sep=":"))
  #create data frame with 50 kb intervals we want to search
  df <- data.frame(chrom=annMethResult$chr,start=as.numeric(annMethResult$start)-regionSize,end=as.numeric(annMethResult$end)+regionSize)
  #make a GRanges object from the data frame
  dfG <- makeGRangesFromDataFrame(df)
  rm(df)
  #merge the GRanges with mm10 UCSC known genes
  dfGgenes <- mergeByOverlaps(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), dfG)
  rm(dfG)
  #create a data frame with the results
  dfGgenesDF <- data.frame(gene_id=dfGgenes$gene_id, methInterval = dfGgenes$dfG)
  rm(dfGgenes)
  #merge the GRanges gene results with symbol
  dfGgenesDFsym <- merge(dfGgenesDF,eid, by="gene_id",all.x=T)
  #add a column of unique coord
  dfGgenesDFsym$coord <- do.call(paste, c(list(dfGgenesDFsym$methInterval.seqnames,dfGgenesDFsym$methInterval.start+regionSize,dfGgenesDFsym$methInterval.end-regionSize),sep=":"))
  #collapse data frame by coord and list all genes found within
  result <- aggregate(symbol ~ coord, data = dfGgenesDFsym, paste, collapse = ",")

  annMethResult50kb <- merge(annMethResult,result,by="coord", all.x=T)
  return(annMethResult50kb)

}
