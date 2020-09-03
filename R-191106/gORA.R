#' A Over-representation Analysis Plotting Function
#'
#' This function allows you to  plot results from WebGestaltR ORA results
#' @param no defaults.
#' @keywords ORA
#' @export
#' @examples
#' gORA()

gORA <- function(ORA,db,name,direction,analysis){
 ORA$labs <- paste(ORA$geneSet,ORA$description, sep=": ")
 ORA.Sig <- ORA[ORA$FDR < 0.05,]
 ORA.Sig <- ORA.Sig[order(ORA.Sig$FDR, decreasing= F),]
 ORA.Sig$labs <- factor(ORA.Sig$labs, levels = ORA.Sig$labs[order(-log(ORA.Sig$FDR))])
 p <- ggplot(data=na.omit(ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)),fill=enrichmentRatio)) + #, fill=R
  geom_bar(stat="identity",col="black") +
  scale_fill_gradient(low="lightskyblue", high="tomato") +
  geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3) +
  labs(title=sprintf("ORA Top 50 - %s %s %s",name,direction,analysis), subtitle="FDR < 5%")  +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.y=element_blank())

 ggsave(filename = sprintf("%s_%s_ORA_%s_%s.pdf",name,db,analysis,direction), plot = p, width = 11, height = 11)
}
