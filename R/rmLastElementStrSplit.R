#' A String Manipulation Function
#'
#' This function removes the last element of "." or "_" separated string
#' @param test no default
#' @keywords cats
#' @export
#' @examples
#' rmLastElementStrSplit()

rmLastElementStrSplit <- function(test){
  #includes all elments of a stringsplit except the last
  return(paste(strsplit(test,"[._]+")[[1]][1:length(strsplit(test,"[._]+")[[1]])-1],collapse="_"))
}
