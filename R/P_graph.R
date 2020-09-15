#' Build Probability-Based Network
#'
#' Build probability-based regulator-target interaction network.
#'
#' @param expmat (matrix) Gene expression matrix.
#'
#' @param net (data.frame) Network, with "reg" and "target" in column name.
#'
#' @param sep (numeric) Number of bins for calculating marginal/joint
#' probability.
#'
#' @param method (character) Method for calculating probability-based distance,
#' either PXY-PXPY ("difference") or mutual information ("mi").
#'
#' @param null (ecdf) Null distribution of probability-based distance. Either
#' from random interactions by P_null function, or all interactions in net.
#'
#' @param threshold (numeric) P-value threshold for filtering interactions in
#' net.
#'
#' @return (igraph) Network graph with "pvalue" and "direction", and "pagerank"
#' as edge/vertex attributes.
#'
#' @examples
#' library(bcellViper)
#' data(bcellViper)
#' dset <- exprs(dset)
#' net <- do.call(rbind, lapply(1:10, function(i, regulon){
#'   data.frame(reg=rep(names(regulon)[i], 10),
#'              target=names(regulon[[i]][[1]])[1:10],
#'              direction=rep(1, 10),
#'              stringsAsFactors = FALSE)}, regulon=regulon))
#' P_graph(dset, net, method="difference", null=NULL, threshold=0.05)
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom stats ecdf
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph set_edge_attr
#' @importFrom igraph set_vertex_attr
#' @importFrom igraph page_rank
#'
#' @export

P_graph <- function(expmat, net, sep=5, method=c("difference", "mi"),
                    null=NULL, threshold=1e-3){
    bin <- gene_bin(genes=union(net$reg, net$target), expmat=expmat, sep=sep)
    px <- PX(expmat=expmat, bin=bin)
    pxy <- PXY(expmat=expmat, bin=bin, x=net$reg, y=net$target)
    pxpy <- PXPY(px=px, combinations=names(pxy))
    dist <- P_dist(pxy=pxy, pxpy=pxpy, method=method)
    if (is.null(null)) null <- ecdf(unlist(dist))
    pvalue <- structure(1-null(unlist(dist)), names=names(dist))
    graph <- lapply(names(pvalue)[pvalue <= threshold], function(x){
        unlist(strsplit(x, split = "_"))})
    graph <-  do.call(rbind, graph)
    graph <- graph_from_data_frame(graph[, 2:1], directed=TRUE)
    graph <- set_edge_attr(graph=graph, name="pvalue",
                           value=pvalue[pvalue <= threshold])
    graph <- set_edge_attr(graph=graph, name="direction",
                           value=net$direction[pvalue <= threshold])
    graph <- set_vertex_attr(graph=graph, name="pagerank",
                             value=page_rank(graph)$vector)
    return(graph)}
