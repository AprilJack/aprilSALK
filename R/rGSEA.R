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
  if(analysis == "Seurat"){
    igt <- "genesymbol"
    mSigDB_2019 <- c("/gpfs/analyses/april/April/ref/mSigDB_2019/c1.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c2.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c3.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c4.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c5.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c6.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c7.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/h.all.v7.0.symbols.gmt")

    mSigDB_2019_des <- c("/gpfs/analyses/april/April/ref/mSigDB_2019/c1.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c2.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c3.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c4.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c5.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c6.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c7.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/h.all.v7.0.symbols.des")

    GSEA <- tryCatch(WebGestaltR(enrichMethod="GSEA",
                                 organism="hsapiens",
                                 enrichDatabaseFile = mSigDB_2019,
                                 enrichDatabaseDescriptionFile = mSigDB_2019_des,
                                 enrichDatabaseType = "genesymbol",
                                 interestGeneType=igt,
                                 projectName=sprintf("%s.%s.%s.GSEA",name,db,analysis),
                                 interestGene=gsea.df,
                                 nThread = 16 ), error = function(err) {
                                   # error handler picks up where error was generated
                                   print(paste("MY_ERROR:  ",err))
                                   GSEA=NULL
                                   return(GSEA)
                                 })
  }else if(db == "mSigDB"){
    print("WTF")
    mSigDB_2019 <- c("/gpfs/analyses/april/April/ref/mSigDB_2019/c1.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c2.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c3.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c4.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c5.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c6.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/c7.all.v7.0.symbols.gmt",
                     "/gpfs/analyses/april/April/ref/mSigDB_2019/h.all.v7.0.symbols.gmt")

    mSigDB_2019_des <- c("/gpfs/analyses/april/April/ref/mSigDB_2019/c1.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c2.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c3.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c4.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c5.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c6.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/c7.all.v7.0.symbols.des",
                         "/gpfs/analyses/april/April/ref/mSigDB_2019/h.all.v7.0.symbols.des")

    GSEA <- tryCatch(WebGestaltR(enrichMethod="GSEA",
                                 organism="hsapiens",
                                 enrichDatabaseFile = mSigDB_2019,
                                 enrichDatabaseDescriptionFile = mSigDB_2019_des,
                                 enrichDatabaseType = "genesymbol",
                                 interestGeneType=igt,
                                 projectName=sprintf("%s.%s.%s.GSEA",name,db,analysis),
                                 interestGene=gsea.df,
                                 nThread = 16 ), error = function(err) {
                                   # error handler picks up where error was generated
                                   print(paste("MY_ERROR:  ",err))
                                   GSEA=NULL
                                   return(GSEA)
                                 })
    }
  else{
      print("here_in_GSEA2")
      GSEA <- tryCatch({WebGestaltR(enrichMethod="GSEA",
                               organism=org,
                               enrichDatabase=db,
                               interestGeneType=igt,
                               projectName=sprintf("%s.%s.%s.GSEA",name,db,analysis),
                               interestGene=gsea.df, nThread = 16 )}, error = function(err) {
                                 # error handler picks up where error was generated
                                 print(paste("MY_ERROR:  ",err))
                                 GSEA=NULL
                                 return(GSEA)})
    }
  return(GSEA)
}
