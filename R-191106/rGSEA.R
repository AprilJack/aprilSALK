#' A Gene Set Enrichment  Function
#'
#' This function allows you to preform GSEA analysis
#' @param no defaults.
#' @keywords GSEA
#' @export
#' @examples
#' rGSEA()

rGSEA <- function(gsea.df,org,db,name,analysis){
  igt <- "refseq_mrna"
  if(analysis == "Seurat"){ igt <- "genesymbol"}
  GSEA <- tryCatch(WebGestaltR(enrichMethod="GSEA",
                               organism=org,
                               enrichDatabase=db,
                               interestGeneType=igt,
                               projectName=sprintf("%s.%s.%s.GSEA",name,db,analysis),
                               interestGene=gsea.df, nThread = 16 ), error = function(err) {
                                 # error handler picks up where error was generated
                                 print(paste("MY_ERROR:  ",err))
                                 GSEA=NULL
                                 return(GSEA)
                               })
}
