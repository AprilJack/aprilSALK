#' A Gene Set Enrichment Plotting Function
#'
#' This function allows you to  plot results from WebGestaltR GSEA results
#' @param no defaults.
#' @keywords GSEA
#' @export
#' @examples
#' gGSEA()

gGSEA <- function(GSEA,db,name,analysis){
  if(db == "mSigDB"){
    GSEA$labs <- GSEA$geneSet
  }else{
  GSEA$labs <- paste(GSEA$geneSet,GSEA$description, sep=": ")
  }
  GSEA.Sig <- GSEA[GSEA$FDR < 0.05,]
  GSEA.Sig <- GSEA.Sig[order(GSEA.Sig$normalizedEnrichmentScore, decreasing= F),]
  GSEA.Sig$labs <- factor(GSEA.Sig$labs, levels = GSEA.Sig$labs[order(GSEA.Sig$normalizedEnrichmentScore, decreasing= F)])
  p.GSEA <- ggplot(data=na.omit(GSEA.Sig[1:50,]), aes(x=labs, y=normalizedEnrichmentScore, fill=normalizedEnrichmentScore)) +
    geom_bar(stat="identity",col="black") +
    scale_fill_gradient(low="lightskyblue", high="tomato") +
    geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3)+
    labs(title=sprintf("GSEA Top 50 - %s %s %s",name,db,analysis), subtitle="FDR < 5%")  +
    coord_flip() +
    theme_minimal() +
    theme(axis.title.y=element_blank())

  ggsave(filename = sprintf("%s_%s_GSEA_%s.pdf",name,db,analysis), plot = p.GSEA, width = 11, height = 11)
}
