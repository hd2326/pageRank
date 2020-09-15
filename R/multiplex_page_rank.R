#' Calculate Multiplex PageRank
#'
#' Calculate multiplex PageRank following defination by Halu et al.
#'
#' @references Halu, Arda, et al. "Multiplex pagerank." PloS one 8.10 (2013).
#'
#' @param graph (igraph) The base graph with pagerank and name as vertex
#' attributes.
#'
#' @param ... (igraph) Supporter graphs with pagerank and name as vertex
#' attributes.
#'
#' @param beta (numeric) Parameters for adjusting supporter graph PageRank
#' values.
#' For the same nodes, PageRank values from different supporter graphs will
#' first be multiplicated.
#' The products will then be exponentiate by beta and gamma, as outgoing edge
#' weights and personalizations of the base graph. Four special multiplex
#' PageRank forms are defined by varying (beta, gamma), including additive
#' (0, 1), multiplicative (1, 0), combined (1, 1) and neutral (0, 0).
#' 
#' @param gamma (numeric) Parameters for adjusting supporter graph PageRank
#' values. For the same nodes, PageRank values from different supporter graphs
#' will first be multiplicated. The products will then be exponentiate by beta
#' and gamma, as outgoing edge weights and personalizations of the base graph.
#' Four special multiplex PageRank forms are defined by varying (beta, gamma),
#' including additive (0, 1), multiplicative (1, 0), combined (1, 1) and
#' neutral (0, 0).
#' 
#' @param damping (numeric) Damping factor.
#'
#' @return (numeric) Multiplex PageRank values.
#'
#' @examples
#' library(igraph)
#' set.seed(1)
#' graph1 <- igraph::erdos.renyi.game(100, 0.01, directed = TRUE)
#' igraph::V(graph1)$name <- 1:100
#' igraph::V(graph1)$pagerank <- igraph::page_rank(graph1)$vector
#' set.seed(2)
#' graph2 <- igraph::erdos.renyi.game(100, 0.01, directed = TRUE)
#' igraph::V(graph2)$name <- 1:100
#' igraph::V(graph2)$pagerank <- igraph::page_rank(graph2)$vector
#' multiplex_page_rank(graph1, graph2)
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom igraph V
#' @importFrom igraph page_rank
#' @importFrom igraph as_data_frame
#'
#' @export

multiplex_page_rank <- function (graph, ..., beta=1, gamma=1, damping=0.85){
    pagerank1 <- structure(V(graph)$pagerank, names=names(V(graph)))
    pagerank2 <- lapply(list(...), function(g) structure(V(g)$pagerank,
                                                         names=names(V(g))))
    pagerank2 <- do.call(rbind, lapply(pagerank2, function(x, names){
        pr <- structure(rep(0, length(names)), names=names)
        pr[intersect(names, names(x))] <- x[intersect(names, names(x))]
        pr},
        names=intersect(names(pagerank1),
                        Reduce(union, lapply(pagerank2, function(x){
                            names(x)})))))
    #pagerank values
    personalized <- structure(rep(0, length(pagerank1)), names=names(pagerank1))
    personalized[colnames(pagerank2)] <- apply(pagerank2, 2, prod)
    personalized[personalized == 0] <- min(personalized[personalized > 0])
    personalized <- personalized^gamma
    #personalization for multiplex pagerank as in halu et al.
    weights <- unlist(lapply(as_data_frame(graph)$to, function(x, pagerank2){
        if (x %in% colnames(pagerank2)) prod(pagerank2[, x])
        else 0}, pagerank2=pagerank2))
    weights[weights == 0] <- min(weights[weights > 0])
    weights <- weights^beta
    #weights for multiplex pagerank as in halu et al.
    pagerank <- page_rank(graph, directed=TRUE, damping=damping,
                          personalized=personalized, weights=weights)$vector
    return(pagerank)}
