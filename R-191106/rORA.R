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

 if(analysis == "Seurat"){ igt <- "genesymbol"; rgt <- "genesymbol" }

 ORA <- WebGestaltR(enrichMethod="ORA",
                      organism=org,
                      enrichDatabase=db,
                      interestGeneType=igt,
                      referenceGeneType=rgt,
                      projectName=sprintf("%s.%s.%s.%s.ORA",name,db,analysis,direction),
                      referenceGene=rg,
                      interestGene=ig)
 return(ORA)
}
