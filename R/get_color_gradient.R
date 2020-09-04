#' Generate Color Gradient
#'
#' Generate color gradient for, e.g. gene expression.
#'
#' @param x (numeric) Vector based on which color gradient is generated.
#' 
#' @param col (character) Color vector.
#' 
#' @param breaks (numeric) A set of breakpoints for the colors. Must be the same length of col.
#' 
#' @return (character) Colors.
#'
#' @examples
#' library(pageRank)
#' get_color_gradient(-2:2)
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#' 
#' @import grDevices
#'
#' @export

get_color_gradient <- function (x, col=colorRampPalette(c("Blue", "Grey", "Red"))(100), breaks=seq(-2, 2, length.out=100)){
  color <- NULL
  if (length(col)  == length(breaks)) color <- col[unlist(lapply(x, function(xx, breaks) which.min(abs(xx - breaks)), breaks = breaks))]
  else message("col and breaks should be the same length")
  return(color)}
