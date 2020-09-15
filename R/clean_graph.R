#' Clean Graph
#'
#' Remove graph nodes by residing subgraph sizes, vertex names and PageRank
#' values.
#'
#' @param graph (igraph) The graph to be cleaned.
#'
#' @param size (numeric) Subgraph size cutoff.
#'
#' @param vertices (character) Vertices to be kept.
#'
#' @param pagerank (numeric) PageRank cutoff.
#'
#' @return (igraph) Network updated "pagerank" as vertex attribute.
#'
#' @examples
#' library(igraph)
#' set.seed(1)
#' graph <- igraph::erdos.renyi.game(100, 0.01, directed = TRUE)
#' igraph::V(graph)$name <- 1:100
#' igraph::V(graph)$pagerank <- igraph::page_rank(graph)$vector
#' clean_graph(graph, size=5)

#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom igraph V
#' @importFrom igraph delete_vertices
#' @importFrom igraph clusters
#' @importFrom igraph set_vertex_attr
#' @importFrom igraph page_rank
#'
#' @export

clean_graph <- function(graph, size=NULL, vertices=NULL, pagerank=NULL){
    v <- v1 <- v2 <- names(V(graph))
    if (!is.null(vertices)) v1 <- intersect(v, vertices)
    if (!is.null(pagerank)) v2 <- v[V(graph)$pagerank >= pagerank]
    graph <- delete_vertices(graph, setdiff(v, intersect(v1, v2)))
    if (!is.null(size)){
        cluster <- igraph::clusters(graph)
        v <- names(V(graph))[cluster$membership%in%which(cluster$csize < size)]
        graph <- delete_vertices(graph, v)}
    graph <- set_vertex_attr(graph=graph, name="pagerank",
                             value=page_rank(graph)$vector)
    return(graph)}
