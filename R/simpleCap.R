#' A String Manipulation Function
#'
#' This function changes case of a word to simple caps. Exp cats -> Cats, or CATS -> Cats.
#' @param no default
#' @keywords string
#' @export
#' @examples
#' simpleCap()



simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste0(toupper(substring(s, 1,1)), tolower(substring(s, 2)))
}
