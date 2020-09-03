#' A Function to Test for DE WITHIN Seurat Clusters
#'
#' This function allows you to find DE genes between two groups within a cluster, and then run GO analysis.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords seurat, scRNA, singlecell, DE
#' @export
#' @examples
#' withinClustDE()

withinClustDE <- function(seuratObj, clustIdent, groupBy="orig.ident", identGroup1, identGroup2,
                          minPCT=.25, geneAnno, nam, org,
                          db=c("geneontology_Biological_Process","geneontology_Molecular_Function","pathway_KEGG")){
  #differential expression testing between replicates WITHIN each cluster
  print(sprintf("Finding DE for Cluster %s",clustIdent))
  clust <- SubsetData(seuratObj, ident.use = clustIdent)
  markersClust <- FindMarkers(object = clust,
                              group.by = groupBy,
                              ident.1 = identGroup1,
                              ident.2 = identGroup2,
                              min.pct = minPCT)
  markersClust1 <- merge(markersClust, geneAnno,by.x = 0,by.y = "gene_name",all.x=T)
  write.table(x = markersClust1,file = sprintf("%s_results.txt",nam),quote = F,row.names = F,sep="\t")
  print(table(markersClust$p_val_adj < 0.05))
  print(markersClust[markersClust$p_val_adj < 0.05,])
  sapply(db,webgX <- function(x) {
    webG_2(counts=NULL,
           resSig = markersClust[markersClust$p_val_adj < 0.05,],
           org = org,
           db = x,
           name = nam,
           gsea.df = markersClust )})
}
