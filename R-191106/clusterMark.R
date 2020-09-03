#' Seurat FindMarkers cluster by cluster and annotate them.
#'
#' This function allows you to use Seurat FindMarkers cluster by cluster and annotate them using an outside reference.
#' @param no defaults.
#' @keywords cluster, scRNA, explore, Seurat
#' @export
#' @examples
#' clusterMark()

#function to find distinguishing features for each cluster and annotate them using outside reference.
clusterMark <- function(obj,clust,anno,nam){
  cluster.markers <- FindMarkers(obj, ident.1 = clust, min.pct = 0.25)
  cluster.markers.ann <- merge(cluster.markers, anno,by.x = 0,by.y = "gene_name",all.x=T)
  write.table(x = cluster.markers.ann,file = sprintf("%s_cluster%s_annotated.txt",nam,clust),sep = "\t")
}
