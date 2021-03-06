% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiplex_page_rank.R
\name{multiplex_page_rank}
\alias{multiplex_page_rank}
\title{Calculate Multiplex PageRank}
\usage{
multiplex_page_rank(graph, ..., beta = 1, gamma = 1, damping = 0.85)
}
\arguments{
\item{graph}{(igraph) The base graph with pagerank and name as vertex
attributes.}

\item{...}{(igraph) Supporter graphs with pagerank and name as vertex
attributes.}

\item{beta}{(numeric) Parameters for adjusting supporter graph PageRank
values.
For the same nodes, PageRank values from different supporter graphs will
first be multiplicated.
The products will then be exponentiate by beta and gamma, as outgoing edge
weights and personalizations of the base graph. Four special multiplex
PageRank forms are defined by varying (beta, gamma), including additive
(0, 1), multiplicative (1, 0), combined (1, 1) and neutral (0, 0).}

\item{gamma}{(numeric) Parameters for adjusting supporter graph PageRank
values. For the same nodes, PageRank values from different supporter graphs
will first be multiplicated. The products will then be exponentiate by beta
and gamma, as outgoing edge weights and personalizations of the base graph.
Four special multiplex PageRank forms are defined by varying (beta, gamma),
including additive (0, 1), multiplicative (1, 0), combined (1, 1) and
neutral (0, 0).}

\item{damping}{(numeric) Damping factor.}
}
\value{
(numeric) Multiplex PageRank values.
}
\description{
Calculate multiplex PageRank following defination by Halu et al.
}
\examples{
library(igraph)
set.seed(1)
graph1 <- igraph::erdos.renyi.game(100, 0.01, directed = TRUE)
igraph::V(graph1)$name <- 1:100
igraph::V(graph1)$pagerank <- igraph::page_rank(graph1)$vector
set.seed(2)
graph2 <- igraph::erdos.renyi.game(100, 0.01, directed = TRUE)
igraph::V(graph2)$name <- 1:100
igraph::V(graph2)$pagerank <- igraph::page_rank(graph2)$vector
multiplex_page_rank(graph1, graph2)

}
\references{
Halu, Arda, et al. "Multiplex pagerank." PloS one 8.10 (2013).
}
\author{
DING, HONGXU (hd2326@columbia.edu)
}
