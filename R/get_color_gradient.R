#' Generate Color Gradient
#'
#' Generate color gradient for, e.g. gene expression.
#'
#' @param x (numeric) Vector based on which color gradient is generated.
#'
#' @return (character) Colors.
#'
#' @examples
#' get_color_gradient(-2:2)
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @export

get_color_gradient <- function (x, col=colorRampPalette(c("Blue", "Grey", "Red"))(100), breaks=seq(-2, 2, length.out=100)){
  color <- NULL
  if (length(col)  == length(breaks)) color <- col[unlist(lapply(x, function(xx, breaks) which.min(abs(xx - breaks)), breaks = breaks))]
  else message("col and breaks should be the same length")
  return(color)}
