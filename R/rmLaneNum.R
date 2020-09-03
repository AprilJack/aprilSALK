#' A Lane Number Removal Function
#'
#' This function removes the last 5 characters from any string.
#' @param no defaults.
#' @keywords five
#' @export
#' @examples
#' rmLaneNum()


rmLaneNum <- function(test){
  #works to remove last 5 characters from string
  n <- nchar(test)
  n <- n-5
  rt <- substr(test,1,n)
  return(rt)
}
