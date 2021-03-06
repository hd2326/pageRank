% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjust_graph.R
\name{adjust_graph}
\alias{adjust_graph}
\title{Re-calculate PageRank}
\usage{
adjust_graph(graph, damping = 0.85, personalized = NULL, weights = NULL)
}
\arguments{
\item{graph}{(igraph) The graph to be adjusted.}

\item{damping}{(numeric) Damping factor.}

\item{personalized}{(numeric) Personalized vector.}

\item{weights}{(numeric) Weight vector.}
}
\value{
(igraph) Network with updated "pagerank" as vertex attribute.
}
\description{
Re-calculate PageRank with updated damping factor, personalized vector and
edge weights.
}
\examples{
library(igraph)
set.seed(1)
graph <- igraph::erdos.renyi.game(100, 0.01, directed = TRUE)
igraph::V(graph)$name <- 1:100
igraph::V(graph)$pagerank <- igraph::page_rank(graph, damping=0.85)$vector
adjust_graph(graph, damping=0.1)

}
\author{
DING, HONGXU (hd2326@columbia.edu)
}
