#' Calculate Temporal PageRank from Two Graphs
#'
#' Calculate temporal PageRank by changing edges between graph1 and graph2.
#' This is a simplified version of temporal PageRank described by Rozenshtein and Gionis,
#' by only analyzing temporally adjacent graph pairs.
#'
#' @references Rozenshtein, Polina, and Aristides Gionis. "Temporal pagerank." Joint European Conference on Machine Learning and Knowledge Discovery in Databases. Springer, Cham, 2016.
#'
#' @param graph1 (igraph) The 1st graph.
#'
#' @param graph2 (igraph) The 2nd graph.
#'
#' @return (igraph) Network graph1-graph2 with "moi (mode of interaction)" and "pagerank" as edge and vertex attributes.
#'
#' @examples
#' library(pageRank)
#' set.seed(1)
#' graph1 <- igraph::erdos.renyi.game(100, 0.01, directed = TRUE)
#' V(graph1)$name <- 1:100
#' set.seed(2)
#' graph2 <- igraph::erdos.renyi.game(100, 0.01, directed = TRUE)
#' V(graph2)$name <- 1:100
#' diff_graph(graph1, graph2)
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph set_edge_attr
#' @importFrom igraph set_vertex_attr
#' @importFrom igraph page_rank
#'
#' @export

diff_graph <- function(graph1, graph2){
  graph1 <- as.matrix(as_adjacency_matrix(graph1))
  graph2 <- as.matrix(as_adjacency_matrix(graph2))
  vertices <- union(rownames(graph1), rownames(graph2))
  g1 <- g2 <- matrix(0, length(vertices), length(vertices), dimnames = list(vertices, vertices))
  g1[rownames(graph1), colnames(graph1)] <- graph1
  g2[rownames(graph2), colnames(graph2)] <- graph2
  graph <- g1 - g2
  graph <- do.call(rbind, lapply(seq_len(nrow(graph)), function(i, graph){
    data.frame(target=rep(rownames(graph)[i], sum(graph[i, ] != 0)),
               reg=colnames(graph)[graph[i, ] != 0],
               moi=graph[i, graph[i, ] != 0])
  }, graph=graph))
  graph <- graph_from_data_frame(graph, directed=TRUE)
  graph <- set_vertex_attr(graph=graph, name="pagerank", value=page_rank(graph)$vector)
  return(graph)}
