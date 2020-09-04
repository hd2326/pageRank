#' Re-calculate PageRank
#'
#' Re-calculate PageRank with updated damping factor, personalized vector and edge weights.
#'
#' @param graph (igraph) The graph to be adjusted.
#'
#' @param damping (numeric) Damping factor.
#'
#' @param personalized (numeric) Personalized vector.
#'
#' @param weights (numeric) Weight vector.
#'
#' @return (igraph) Network with updated "pagerank" as vertex attribute.
#'
#' @examples
#' library(igraph)
#' set.seed(1)
#' graph <- igraph::erdos.renyi.game(100, 0.01, directed = TRUE)
#' V(graph)$name <- 1:100
#' V(graph)$pagerank <- igraph::page_rank(graph, damping=0.85)$vector
#' adjust_graph(graph, damping=0.1)
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom igraph page_rank
#' @importFrom igraph set_vertex_attr
#'
#' @export

adjust_graph <- function(graph, damping=0.85, personalized=NULL, weights=NULL){
  pagerank <- page_rank(graph, damping=damping, personalized=personalized, weights=weights)$vector
  graph <- set_vertex_attr(graph=graph, name="pagerank", value=pagerank)
  return(graph)}
