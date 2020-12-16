#this script is used for complexity analysis

library(igraph)
library(pageRank)

graph <- lapply(2^(5:14), function(n, p){
  lapply(1:11, function(i, n, p){
    set.seed(n+i)
    graph <- igraph::erdos.renyi.game(n, p)
    graph <- igraph::set_vertex_attr(graph, "name", value = 1:n)
    graph <- igraph::set_vertex_attr(graph, "pagerank", value = igraph::page_rank(graph)$vector)
    graph}, n=n, p=p)}, p=0.1)
#random ER graphs
t_multiplex <- lapply(1:10, function(n, graph){
  lapply(1:10, function(i, n, graph){
    t1 <- Sys.time()
    pr <- multiplex_page_rank(graph[[n]][[11]], graph[[n]][[i]])
    t2 <- Sys.time()
    as.numeric(t2-t1, unit = "secs")
  }, n=n, graph=graph)}, graph=graph)
#multiplex PageRank complexity
t_temporal <- lapply(1:10, function(n, graph){
  lapply(1:10, function(i, n, graph){
    t1 <- Sys.time()
    graph <- diff_graph(graph[[n]][[11]], graph[[n]][[i]])
    t2 <- Sys.time()
    as.numeric(t2-t1, unit = "secs")
  }, n=n, graph=graph)}, graph=graph)
#temporal PageRank complexity
save(t_multiplex, t_temporal, file = "complexity_analysis.rda")

pdf("SupFig10.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(5, 6, 4, 1))
bp1 <- do.call(cbind, lapply(t_multiplex, function(x) do.call(rbind, x)))
boxplot(log2(bp1), outline = F, xaxt = "n", cex.axis = 2, xlab = "n, log2", ylab = "seconds, log2", cex.lab = 2, main = "multiplex PageRank", cex.main = 2)
axis(side = 1, at = 1:10, labels = 5:14, cex.axis = 2)
legend("topleft", legend = c("ER(n, p=0.1)", "N(network)=2"), cex = 2, bty = "n")
bp2 <- do.call(cbind, lapply(t_temporal, function(x) unlist(x)))
boxplot(log2(bp2), outline = F, xaxt = "n", cex.axis = 2, xlab = "n, log2", ylab = "seconds, log2", cex.lab = 2, main = "temporal PageRank", cex.main = 2)
axis(side = 1, at = 1:10, labels = 5:14, cex.axis = 2)
legend("topleft", legend = c("ER(n, p=0.1)", "N(network)=2"), cex = 2, bty = "n")
dev.off()
#SupFig10