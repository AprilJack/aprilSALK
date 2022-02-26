#' A Venn Diagram GO Analysis Function
#'
#' This function takes in a Venn diagram object from VennDiagram package, and runs WebGestaltR GO annotation.
#' @param no defaults
#' @keywords venn, GO, WebGestalt
#' @export
#' @examples
#' myGOKmeans_venn()

myGOKmeans_venn <- function(set_num,counts=NULL,VennDiagram_obj,name,org,db){
  library("WebGestaltR")

  if(!identical(VennDiagram_obj$..values..[[set_num]],character(0))){
    if(db=="MSigDB"){
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
      tryCatch(Up.ORA <- WebGestaltR::WebGestaltR(enrichMethod = "ORA",
                                                organism = "hsapiens",
                                                enrichDatabaseFile = mSigDB_2019,
                                                enrichDatabaseDescriptionFile = mSigDB_2019_des,
                                                enrichDatabaseType = "genesymbol",
                                                interestGene =  toupper(as.character(as.data.frame(VennDiagram_obj$..values..[set_num])[,1])),
                                                interestGeneType = "genesymbol",
                                                referenceSet =  "genome_protein-coding",
                                                referenceGeneType = "genesymbol",
                                                nThreads = 16,
                                                projectName=sprintf("%s.%s.set.%s.ORA",name,db,VennDiagram_obj$..set..[set_num])),error=function(e) return(NULL) )
    }
    else{
    Up.ORA <- tryCatch(WebGestaltR(enrichMethod="ORA",
                          organism=org,
                          enrichDatabase ="geneontology_Biological_Process",
                          interestGeneType="genesymbol",
                          #referenceGeneType="refseq_mrna",
                          projectName=sprintf("%s.set.%s.ORA",name,VennDiagram_obj$..set..[set_num]),
                          referenceSet =  "genome_protein-coding",
                          #referenceGene=as.vector(row.names(counts)),
                          interestGene=as.character(as.data.frame(VennDiagram_obj$..values..[set_num])[,1])),error=function(e) return(NULL) )
    }
    if (!is.null(dim(Up.ORA)[1])) {

      Up.ORA$labs <- paste(Up.ORA$geneSet,Up.ORA$description, sep=": ")
      Up.ORA.Sig <- Up.ORA[Up.ORA$FDR < 0.05,]
      Up.ORA.Sig <- Up.ORA.Sig[order(Up.ORA.Sig$FDR, decreasing= F),]
      Up.ORA.Sig$labs <- factor(Up.ORA.Sig$labs, levels = Up.ORA.Sig$labs[order(-log(Up.ORA.Sig$FDR))])
      p.Up <- ggplot(data=na.omit(Up.ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=enrichmentRatio)) +
        geom_bar(stat="identity",col="black") +
        scale_fill_gradient(low="lightskyblue", high="tomato") +
        geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3) +
        labs(title=sprintf("ORA Top 50 - %s Set %s",name,VennDiagram_obj$..set..[set_num]),
             subtitle="FDR < 5%")  +
        coord_flip() +
        theme_minimal() +
        theme(axis.title.y=element_blank())

      ggsave(filename = sprintf("%s_%s_ORA_Set_%s.pdf",name,db,VennDiagram_obj$..set..[set_num]), plot = p.Up, width = 11, height = 11)

    }
   return(Up.ORA)
  }
}
