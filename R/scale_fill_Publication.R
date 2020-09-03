#' A Graphing Function
#'
#' This function makes ggplot's theme resemble graphpad.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords ggplot graphing graphpad
#' @export
#' @examples
#' scale_fill_Publication()

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}
