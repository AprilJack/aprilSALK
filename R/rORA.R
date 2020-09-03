#' A Over-representation Analysis Function
#'
#' This function allows you to preform ORA analysis
#' @param no defaults.
#' @keywords ORA
#' @export
#' @examples
#' rORA()

rORA <- function(counts,resSig,org,db,name,direction,analysis){

 igt <- "refseq_mrna"
 rgt <- "refseq_mrna"
 if(direction == "UP"){ ig <- row.names(resSig[resSig$logFC > 0,]); rg <- row.names(counts)}
 if(direction == "DN"){ ig <- row.names(resSig[resSig$logFC < 0,]); rg <- row.names(counts)}
 if(direction == "Combined"){ ig <- row.names(resSig); rg <- row.names(counts)}
 if(direction == "Hyper"){ ig <- unique(resSig[resSig$meth.diff > 25,]$Nearest.Refseq); rg <- referenceGene=unique(counts$Nearest.Refseq)}
 if(direction == "Hypo"){ ig <- unique(resSig[resSig$meth.diff < -25,]$Nearest.Refseq); rg <- referenceGene=unique(counts$Nearest.Refseq)}

 if(analysis == "Seurat"){
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
         ORA <- WebGestaltR::WebGestaltR(enrichMethod = "ORA",
                                              organism = "hsapiens",
                                              enrichDatabaseFile = mSigDB_2019,
                                              enrichDatabaseDescriptionFile = mSigDB_2019_des,
                                              enrichDatabaseType = "genesymbol",
                                              interestGene =  toupper(ig),
                                              interestGeneType = "genesymbol",
                                              referenceSet =  "genome_protein-coding",
                                              referenceGeneType = "genesymbol",
                                              nThreads = 16,
                                              projectName = sprintf("%s.%s.%s.%s.ORA",name,db,analysis,direction))
         return(ORA)
     }else if(db == "mSigDB" & org == "hsapiens"){
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
             tryCatch({ORA <- WebGestaltR::WebGestaltR(enrichMethod = "ORA",
                                             organism = "hsapiens",
                                             enrichDatabaseFile = mSigDB_2019,
                                             enrichDatabaseDescriptionFile = mSigDB_2019_des,
                                             enrichDatabaseType = "genesymbol",
                                             interestGene =  toupper(ig),
                                             interestGeneType = igt,
                                             referenceSet =  "genome_protein-coding",
                                             referenceGeneType = "genesymbol",
                                             nThreads = 16,
                                             projectName=sprintf("%s.%s.%s.%s.ORA",name,db,analysis,direction))
                       return(ORA)},error=function(e) return(NULL) )

         }
           else{
           tryCatch({
                 ORA <- WebGestaltR(enrichMethod="ORA",
                      organism=org,
                      enrichDatabase=db,
                      interestGeneType=igt,
                      referenceGeneType=rgt,
                      projectName=sprintf("%s.%s.%s.%s.ORA",name,db,analysis,direction),
                      referenceGene=rg,
                      nThreads = 16,
                      interestGene=ig)
                 return(ORA)},error=function(e) return(NULL) )

         }
}
