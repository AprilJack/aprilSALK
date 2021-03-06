#' A Graphing Function
#'
#' This function makes ggplot's theme resemble graphpad.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords ggplot graphing graphpad
#' @export
#' @examples
#' scale_colour_Publication()


scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}
