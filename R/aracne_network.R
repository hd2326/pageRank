#' Re-format ARACNe Network.
#'
#' Re-format ARACNe network in regulon object to data.frame with regulator, target and direction columns.
#'
#' @param regulon (regulon) ARACNe network.
#'
#' @return (data.frame) Network, with "reg", "target" and "direction" in column name. For direction, 1/0 denotes positive/negative regulation.
#'
#' @examples
#' library(bcellViper)
#' data(bcellViper)
#' aracne_network(regulon[1:10])
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @export

aracne_network <- function(regulon){
  regulon <- lapply(seq_len(length(regulon)), function(i, regulon) data.frame(reg=rep(names(regulon)[i], length(regulon[[i]]$tfmode)),
                                                                              target=names(regulon[[i]]$tfmode),
                                                                              direction=as.integer(regulon[[i]]$tfmode > 0),
                                                                              stringsAsFactors=FALSE), regulon=regulon)
  regulon <- do.call(rbind, regulon)
  return(regulon)}
