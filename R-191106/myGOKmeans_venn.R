#' A Venn Diagram GO Analysis Function
#'
#' This function takes in a Venn diagram object from VennDiagram package, and runs WebGestaltR GO annotation.
#' @param no defaults
#' @keywords venn, GO, WebGestalt
#' @export
#' @examples
#' myGOKmeans_venn()

myGOKmeans_venn <- function(set_num,counts,VennDiagram_obj,name,org){
  library("WebGestaltR")
  
  if(!identical(VennDiagram_obj$..values..[[set_num]],character(0))){
    Up.ORA <- WebGestaltR(enrichMethod="ORA", 
                          organism=org,
                          enrichDatabase ="geneontology_Biological_Process",
                          interestGeneType="genesymbol",
                          referenceGeneType="refseq_mrna",
                          projectName=sprintf("%s.set.%s.ORA",name,VennDiagram_obj$..set..[set_num]),
                          referenceGene=as.vector(row.names(counts)),
                          interestGene=as.character(as.data.frame(VennDiagram_obj$..values..[set_num])[,1]))      
    
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
      
      ggsave(filename = sprintf("%s_ORA_Set_%s.pdf",name,VennDiagram_obj$..set..[set_num]), plot = p.Up, width = 11, height = 11) 
      
    }
  }
}