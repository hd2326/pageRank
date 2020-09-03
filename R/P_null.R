#' Build Null Distribution of Probability-Based Distance
#'
#' Build null model for evaluating the significance of interactions by generating random regulator-target pairs.
#'
#' @param expmat (matrix) Gene expression matrix.
#'
#' @param net (data.frame) Network, with "reg" and "target" in column name.
#'
#' @param n (numeric) Number of random pairs.
#'
#' @param sep (numeric) Number of bins for calculating marginal/joint probability.
#'
#' @param method (character) Method for calculating probability-based distance, either PXY-PXPY ("difference") or mutual information ("mi").
#'
#' @return (ecdf) ECDF of null distribution.
#'
#' @examples
#' library(bcellViper)
#' data(bcellViper)
#' dset <- exprs(dset)
#' net <- do.call(rbind, lapply(1:10, function(i, regulon)
#'   data.frame(reg=rep(names(regulon)[i], 10),
#'              target=names(regulon[[i]][[1]])[1:10],
#'              stringsAsFactors = F), regulon=regulon))
#' P_null(dset, net, n=100)
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom stats ecdf
#'
#' @export

P_null <- function(expmat, net, n=10000, sep=5, method=c("difference", "mi")){
  net <- data.frame(reg=sample(unique(net$reg), size=n, replace=T),
                    target=sample(unique(net$target), size=n, replace=T), stringsAsFactors=F)
  bin <- gene_bin(genes=union(net$reg, net$target), expmat=expmat, sep=sep)
  px <- PX(expmat=expmat, bin=bin)
  pxy <- PXY(expmat=expmat, bin=bin, x=net$reg, y=net$target)
  pxpy <- PXPY(px=px, combinations=names(pxy))
  dist <- P_dist(pxy=pxy, pxpy=pxpy, method=method)
  null <- ecdf(unlist(dist))
  return(null)}
