#' A Seurat Plotting Function
#'
#' This function takes in a list of feature plots from Seurat and removes the legends/axes, and plots them in a grid using cowplot.
#' @param no defaults
#' @keywords cats
#' @export
#' @examples
#' noAxesNoLegendFeaturePlot()



noAxesNoLegendFeaturePlot <- function(p){
  library(cowplot)
  for(i in 1:length(p)) {
    p[[i]] <- p[[i]] + NoLegend() + NoAxes()
  }
  return(cowplot::plot_grid(plotlist = p))
}
