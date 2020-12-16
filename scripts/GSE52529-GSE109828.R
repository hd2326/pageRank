#this script is used for analyzing the myoblast-muscle developmental process

library(viper)
source("functions.R")

load("./GSE52529/fpkm.rda")
time <- unlist(lapply(strsplit(colnames(fpkm), split = "_"), function(x) x[1]))
regul <- viper::pruneRegulon(get(load("./GSE52529/GSE52529-tf-regulon.rda")), cutoff = 100)
regul <- do.call(rbind, lapply(1:length(regul), function(i, regul) data.frame(reg=rep(names(regul)[i], length(regul[[i]]$tfmode)),
                                                                              target=names(regul[[i]]$tfmode),
                                                                              direction=as.integer(regul[[i]]$tfmode > 0), stringsAsFactors=F), regul=regul))

regul <- regul[regul$reg %in% rownames(fpkm) & regul$target %in% rownames(fpkm), ]
regul <- regul[regul$reg != regul$target, ]
null <- P_null(expmat=fpkm, net=regul, n=10000, sep=5, method="difference")
graph <- structure(lapply(unique(time), function(t, time, fpkm, net, null){
  message(t)
  P_graph(expmat=fpkm[, time==t], sep=5, net=net, method="difference", null=null, threshold=1e-3)
}, time=time, fpkm=fpkm, net=regul, null=NULL), names=unique(time))
save(graph, file = "./GSE52529/graph.rda")
#expression-based networks

load("./GSE52529/fpkm.rda")
time <- unlist(lapply(strsplit(colnames(fpkm), split = "_"), function(x) x[1]))
table <- get(load("./GSE109828/table.rda"))
regulators <- read.table("./source/tf-homo-current-symbol.dat", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
regul <- lapply(table, function(t, regulators) atac_network(table=t, regulators=regulators, species="hg19", upstream=2000, downstream=200), regulators=regulators)
names(regul) <- c("T0", "T24", "T48", "T72")
save(regul, file = "./GSE109828/regul.rda")
graph <- structure(lapply(unique(time), function(t, time, fpkm, net, null){
  message(t)
  net <- net[[t]]
  net <- net[net$reg %in% rownames(fpkm) & net$target %in% rownames(fpkm), ]
  net <- net[net$reg != net$target, ]
  null <- P_null(expmat=fpkm[, time==t], net=net, n=10000, sep=5, method="difference")
  P_graph(expmat=fpkm[, time==t], sep=5, net=net, method="difference", null=null, threshold=1e-3)
}, time=time, fpkm=fpkm, net=regul, null=NULL), names=unique(time))
save(graph, file = "./GSE109828/graph.rda")
#ATAC-Seq networks

graph1 <- get(load("../GSE52529/graph.rda"))
graph2 <- get(load("../GSE109828/graph.rda"))
dgraph1 <- list("T0->T24"=diff_graph(graph1=graph1$T24, graph2=graph1$T0),
                "T24->T48"=diff_graph(graph1=graph1$T48, graph2=graph1$T24),
                "T48->T72"=diff_graph(graph1=graph1$T72, graph2=graph1$T48))
dgraph2 <- list("T0->T24"=diff_graph(graph1=graph2$T24, graph2=graph2$T0),
                "T24->T48"=diff_graph(graph1=graph2$T48, graph2=graph2$T24),
                "T48->T72"=diff_graph(graph1=graph2$T72, graph2=graph2$T48))
load("../GSE52529/fpkm.rda")
time <- unlist(lapply(strsplit(colnames(fpkm), split = "_"), function(x) x[1]))
pdf("Figure2.pdf", width = 7, height = 10)
layout(matrix(c(rep(1, 7), rep(2, 7), rep(3, 7),
                rep(4, 12), rep(5, 9),
                rep(6, 12), rep(7, 9),
                rep(8, 10), rep(10, 1), rep(9, 10)), 4, 21, byrow = T), widths = rep(1, 21), heights = c(7, 8, 2, 6))
dge <- time_expmat(time=time, expmat=fpkm)
dge <- list("T0->T24"=dge[, "T24"]-dge[, "T0"],
            "T24->T48"=dge[, "T48"]-dge[, "T24"],
            "T24->T48"=dge[, "T72"]-dge[, "T48"])
par(mar = c(0, 0, 3, 0))
graph <- clean_graph(graph=graph1$T0, size=5, vertices=NULL, pagerank=NULL)
set.seed(1)
plot(igraph::as.undirected(graph), vertex.color="Grey", vertex.frame.color="Grey",
     vertex.size=V(graph)$pagerank*100+1, vertex.label.cex=V(graph)$pagerank*6+0.5,
     vertex.label.color="Black", edge.color=c("Blue", "Red")[E(graph)$direction+1])
title("T0", cex.main=2)
legend("topleft", legend=c("Pos. Reg.", "Neg. Reg."), fill=c("Red", "Blue"), bty="n", border=NA)
legend("bottomleft", legend=c("Size~PR"), bty="n", border=NA)

graph <- clean_graph(graph=graph1$T24, size=5, vertices=NULL, pagerank=NULL)
set.seed(1)
plot(igraph::as.undirected(graph), vertex.color="Grey", vertex.frame.color="Grey",
     vertex.size=V(graph)$pagerank*100+1, vertex.label.cex=V(graph)$pagerank*6+0.5,
     vertex.label.color="Black", edge.color=c("Blue", "Red")[E(graph)$direction+1])
title("T24", cex.main=2)
legend("topleft", legend=c("Pos. Reg.", "Neg. Reg."), fill=c("Red", "Blue"), bty="n", border=NA)
legend("bottomleft", legend=c("Size~PR"), bty="n", border=NA)

graph <- clean_graph(graph=dgraph1$"T0->T24", size=5, vertices=NULL, pagerank=NULL)
color <- get_color_gradient(x=dge$"T0->T24"[names(V(graph))], col=colorRampPalette(c("Purple", "Grey", "Orange"))(100), breaks=seq(-1, 1, length.out=100))
set.seed(1)
plot(igraph::as.undirected(graph), vertex.color=color, vertex.frame.color=color,
     vertex.size=V(graph)$pagerank*100+1, vertex.label.cex=V(graph)$pagerank*6+0.5,
     vertex.label.color="Black", edge.color=c("Blue", "Grey", "Red")[E(graph)$moi+2])
title("T0->T24", cex.main=2)
legend("topright", legend=c("Gain Reg.", "Loss Reg."), fill=c("Red", "Blue"), bty="n", border=NA)
legend("bottomright", legend=c("Incr. Exp.", "Decr. Exp."), fill=c("Orange", "Purple"), bty="n", border=NA)
legend("bottomleft", legend=c("Size~tPR"), bty="n", border=NA)

par(mar = c(4, 4, 4, 1))
dge <- time_expmat(time=time, expmat=fpkm)
pr <- structure(lapply(names(graph1), function(x, g1, g2) multiplex_page_rank(g1[[x]], g2[[x]], method="combined"), g1=graph1, g2=graph2), names=names(graph1))
gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat) dge[gmat[, x], x], dge=dge, gmat=gmat))
dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=graph1, gmat=gmat))
dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Grey", "Red"))(100), breaks=seq(0, 4, length.out=100), main="Combined PR")  

pr <- structure(lapply(names(dgraph1), function(x, g1, g2) multiplex_page_rank(g1[[x]], g2[[x]], method="combined"), g1=dgraph1, g2=dgraph2), names=names(dgraph1))
gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat){
  xx <- unlist(strsplit(x, split="->"))
  (dge[, xx[1]]-dge[, xx[2]])[gmat[, x]]}, dge=dge, gmat=gmat))
dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=dgraph1, gmat=gmat))
dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Blue", "Grey", "Red"))(100), breaks=seq(-1, 1, length.out=100), main="Combined Temporal PR")

par(mar = c(3, 10, 2, 10))
image(matrix(seq(0, 4, length.out = 100), 100, 1), col = colorRampPalette(c("Grey", "Red"))(100), breaks = seq(0, 4, length.out = 101), axes = F)
axis(side = 1, line = -1, at = seq(0, 1, length.out = 3), tick = F, labels = seq(0, 4, length.out = 3), cex.axis = 1.25)
axis(side = 1, line = 0, at = 0.5, tick = F, labels = "log2(rpm+1)", cex.axis = 1.25)
par(mar = c(3, 6, 2, 6))
image(matrix(seq(-1, 1, length.out = 100), 100, 1), col = colorRampPalette(c("Blue", "Grey", "Red"))(100), breaks = seq(-1, 1, length.out = 101), axes = F)
axis(side = 1, line = -1, at = seq(0, 1, length.out = 3), tick = F, labels = seq(-1, 1, length.out = 3), cex.axis = 1.25)
axis(side = 1, line = 0, at = 0.5, tick = F, labels = "log2(rpm+1)", cex.axis = 1.25)

par(mar = c(3, 4, 0, 2))
table <- structure(lapply(names(graph1), function(x, g1, g2){
  pr <- sort(multiplex_page_rank(g1[[x]], g2[[x]], method="combined"), decreasing = T)[1:20]
  pr1 <- structure(igraph::V(g1[[x]])$"pagerank", names=names(V(g1[[x]])))
  pr2 <- structure(igraph::V(g2[[x]])$"pagerank", names=names(V(g2[[x]])))
  table <- matrix(0, 2, length(pr), dimnames = list(1:2, names(pr)))
  table[1, intersect(names(pr), names(pr1))] <- pr1[intersect(names(pr), names(pr1))]/max(pr1)
  table[2, intersect(names(pr), names(pr2))] <- pr2[intersect(names(pr), names(pr2))]/max(pr2)
  apply(table, 2, prop.table)[2, ]}, g1=graph1, g2=graph2), names=names(graph1))
image(t(do.call(rbind, table)), col = colorRampPalette(c("Grey", "Red"))(100), breaks = seq(0, 1, length.out = 101), axes = F)
for (i in 1:4) text(seq(0, 1, length.out = 20), seq(0, 1, length.out = 4)[i], labels = names(table[[i]]), srt = 90)
axis(side = 1, line = 0, at = seq(0, 1, length.out = 20), labels = paste("#", 1:20, sep=""), tick = F, las = 2)
axis(side = 2, line = 0, at = seq(0, 1, length.out = 4), labels = names(table), tick = F, las = 2)

par(mar = c(6, 5, 0, 1))
table <- structure(lapply(names(dgraph1), function(x, g1, g2){
  pr <- sort(multiplex_page_rank(g1[[x]], g2[[x]], method="combined"), decreasing = T)[1:20]
  pr1 <- structure(igraph::V(g1[[x]])$"pagerank", names=names(V(g1[[x]])))
  pr2 <- structure(igraph::V(g2[[x]])$"pagerank", names=names(V(g2[[x]])))
  table <- matrix(0, 2, length(pr), dimnames = list(1:2, names(pr)))
  table[1, intersect(names(pr), names(pr1))] <- pr1[intersect(names(pr), names(pr1))]/max(pr1)
  table[2, intersect(names(pr), names(pr2))] <- pr2[intersect(names(pr), names(pr2))]/max(pr2)
  apply(table, 2, prop.table)[2, ]}, g1=dgraph1, g2=dgraph2), names=names(dgraph1))
image(t(do.call(rbind, table)), col = colorRampPalette(c("Grey", "Red"))(100), breaks = seq(0, 1, length.out = 101), axes = F)
for (i in 1:3) text(seq(0, 1, length.out = 20), seq(0, 1, length.out = 3)[i], labels = names(table[[i]]), srt = 90)
axis(side = 1, line = 0, at = seq(0, 1, length.out = 20), labels = paste("#", 1:20, sep=""), tick = F, las = 2)
axis(side = 2, line = 0, at = seq(0, 1, length.out = 3), labels = names(table), tick = F)

par(mar = c(5, 0, 5, 0))
image(matrix(seq(0, 1, length.out = 100), 1, 100), col = colorRampPalette(c("Grey", "Red"))(100), breaks = seq(0, 1, length.out = 101), axes = F)
axis(side = 4, line = -1, at = seq(0, 1, length.out = 3), tick = F, labels = seq(0, 1, length.out = 3), cex.axis = 1.25)
axis(side = 4, line = 0, at = 0.5, tick = F, labels = "ATAC-Seq GRN Contribution", cex.axis = 1.25)
dev.off()
#Figure 2

source("./source/functions.R")
graph1 <- get(load("./GSE52529/graph.rda"))
graph2 <- get(load("./GSE109828/graph.rda"))
dgraph1 <- list("T0->T24"=diff_graph(graph1=graph1$T24, graph2=graph1$T0),
                "T24->T48"=diff_graph(graph1=graph1$T48, graph2=graph1$T24),
                "T48->T72"=diff_graph(graph1=graph1$T72, graph2=graph1$T48))
dgraph2 <- list("T0->T24"=diff_graph(graph1=graph2$T24, graph2=graph2$T0),
                "T24->T48"=diff_graph(graph1=graph2$T48, graph2=graph2$T24),
                "T48->T72"=diff_graph(graph1=graph2$T72, graph2=graph2$T48))
load("./GSE52529/fpkm.rda")
time <- unlist(lapply(strsplit(colnames(fpkm), split = "_"), function(x) x[1]))
dge <- time_expmat(time=time, expmat=fpkm)
dge <- list("T0->T24"=dge[, "T24"]-dge[, "T0"],
            "T24->T48"=dge[, "T48"]-dge[, "T24"],
            "T24->T48"=dge[, "T72"]-dge[, "T48"])
pdf("SupFig1.pdf", width = 16, height = 16)
par(mfrow = c(4, 4), mar = c(0, 0, 3, 0))
for (x in names(graph1)){
  g <- clean_graph(graph=graph1[[x]], size=5, vertices=NULL, pagerank=NULL)
  set.seed(1)
  plot(igraph::as.undirected(g), vertex.color="Grey", vertex.frame.color="Grey",
       vertex.size=V(g)$pagerank*200+2, vertex.label.cex=V(g)$pagerank*6+0.5,
       vertex.label.color="Black", edge.color=c("Blue", "Red")[E(g)$direction+1], main=paste("scRNA-Seq", x))}
for (x in names(graph2)){
  g <- clean_graph(graph=graph2[[x]], size=5, vertices=NULL, pagerank=NULL)
  set.seed(1)
  plot(igraph::as.undirected(g), vertex.color="Grey", vertex.frame.color="Grey",
       vertex.size=V(g)$pagerank*200+2, vertex.label.cex=V(g)$pagerank*6+0.5,
       vertex.label.color="Black", edge.color="Grey", main=paste("ATAC-Seq", x))}
for (x in names(dgraph1)){
  g <- clean_graph(graph=dgraph1[[x]], size=5, vertices=NULL, pagerank=NULL)
  color <- get_color_gradient(x=dge[[x]][names(V(g))], col=colorRampPalette(c("Purple", "Grey", "Orange"))(100), breaks=seq(-1, 1, length.out=100))
  set.seed(1)
  plot(igraph::as.undirected(g), vertex.color=color, vertex.frame.color=color,
       vertex.size=V(g)$pagerank*200+2, vertex.label.cex=V(g)$pagerank*6+0.5,
       vertex.label.color="Black", edge.color=c("Blue", "Grey", "Red")[E(g)$moi+2], main=paste("scRNA-Seq", x))}
plot.new()
for (x in names(dgraph2)){
  g <- clean_graph(graph=dgraph2[[x]], size=5, vertices=NULL, pagerank=NULL)
  color <- get_color_gradient(x=dge[[x]][names(V(g))], col=colorRampPalette(c("Purple", "Grey", "Orange"))(100), breaks=seq(-1, 1, length.out=100))
  set.seed(1)
  plot(igraph::as.undirected(g), vertex.color=color, vertex.frame.color=color,
       vertex.size=V(g)$pagerank*200+2, vertex.label.cex=V(g)$pagerank*6+0.5,
       vertex.label.color="Black", edge.color=c("Blue", "Grey", "Red")[E(g)$moi+2], main=paste("ATAC-Seq", x))}
plot.new()
dev.off()
#SupFig1

graph1 <- get(load("./GSE52529/graph.rda"))
graph2 <- get(load("./GSE109828/graph.rda"))
dgraph1 <- list("T0->T24"=diff_graph(graph1=graph1$T24, graph2=graph1$T0),
                "T24->T48"=diff_graph(graph1=graph1$T48, graph2=graph1$T24),
                "T48->T72"=diff_graph(graph1=graph1$T72, graph2=graph1$T48))
dgraph2 <- list("T0->T24"=diff_graph(graph1=graph2$T24, graph2=graph2$T0),
                "T24->T48"=diff_graph(graph1=graph2$T48, graph2=graph2$T24),
                "T48->T72"=diff_graph(graph1=graph2$T72, graph2=graph2$T48))
load("./GSE52529/fpkm.rda")
time <- unlist(lapply(strsplit(colnames(fpkm), split = "_"), function(x) x[1]))
dge <- time_expmat(time=time, expmat=fpkm)
pdf("SupFig3.pdf", width = 10, height = 10)
layout(matrix(1:4, 2, 2, byrow = T), widths = c(4, 3))
par(mar = c(6, 4, 4, 1))
pr <- lapply(graph1, function(g) structure(V(g)$pagerank, names=names(V(g))))
gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat) dge[gmat[, x], x], dge=dge, gmat=gmat))
dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=graph1, gmat=gmat))
dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Grey", "Red"))(100),
            breaks=seq(0, 4, length.out = 100), main="scRNA-Seq PR")

pr <- lapply(dgraph1, function(g) structure(V(g)$pagerank, names=names(V(g))))
gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat){
  xx <- unlist(strsplit(x, split="->"))
  (dge[, xx[1]]-dge[, xx[2]])[gmat[, x]]}, dge=dge, gmat=gmat))
dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=dgraph1, gmat=gmat))
dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Blue", "Grey", "Red"))(100),
            breaks=seq(-1, 1, length.out = 100), main="scRNA-Seq tPR")

pr <- lapply(graph2, function(g) structure(V(g)$pagerank, names=names(V(g))))
gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat) dge[gmat[, x], x], dge=dge, gmat=gmat))
dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=graph2, gmat=gmat))
dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Grey", "Red"))(100),
            breaks=seq(0, 4, length.out = 100), main="ATAC-Seq PR")

pr <- lapply(dgraph2, function(g) structure(V(g)$pagerank, names=names(V(g))))
gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat){
  xx <- unlist(strsplit(x, split="->"))
  (dge[, xx[1]]-dge[, xx[2]])[gmat[, x]]}, dge=dge, gmat=gmat))
dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=dgraph2, gmat=gmat))
dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Blue", "Grey", "Red"))(100),
            breaks=seq(-1, 1, length.out = 100), main="ATAC-Seq tPR")
dev.off()
#SupFig3

source("./source/functions.R")
graph1 <- get(load("./GSE52529/graph.rda"))
graph2 <- get(load("./GSE109828/graph.rda"))
dgraph1 <- list("T0->T24"=diff_graph(graph1=graph1$T24, graph2=graph1$T0),
                "T24->T48"=diff_graph(graph1=graph1$T48, graph2=graph1$T24),
                "T48->T72"=diff_graph(graph1=graph1$T72, graph2=graph1$T48))
dgraph2 <- list("T0->T24"=diff_graph(graph1=graph2$T24, graph2=graph2$T0),
                "T24->T48"=diff_graph(graph1=graph2$T48, graph2=graph2$T24),
                "T48->T72"=diff_graph(graph1=graph2$T72, graph2=graph2$T48))
load("./GSE52529/fpkm.rda")
time <- unlist(lapply(strsplit(colnames(fpkm), split = "_"), function(x) x[1]))
dge <- time_expmat(time=time, expmat=fpkm)
pdf("SupFig4.pdf", width = 14, height = 9)
layout(matrix(1:8, 2, 4, byrow = T), widths = c(4, 3, 4, 3))
par(mar = c(6, 4, 4, 1))
for (method in c("additive", "multiplicative", "combined", "neutral")){
  pr <- structure(lapply(names(graph1), function(x, g1, g2) multiplex_page_rank(g1[[x]], g2[[x]], method=method), g1=graph1, g2=graph2), names=names(graph1))
  gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
  emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat) dge[gmat[, x], x], dge=dge, gmat=gmat))
  dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=graph1, gmat=gmat))
  dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
  bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Grey", "Red"))(100),
              breaks=seq(0, 4, length.out = 100), main=paste(method, "PR"))
  
  pr <- structure(lapply(names(dgraph1), function(x, g1, g2) multiplex_page_rank(g1[[x]], g2[[x]], method=method), g1=dgraph1, g2=dgraph2), names=names(dgraph1))
  gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
  emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat){
    xx <- unlist(strsplit(x, split="->"))
    (dge[, xx[1]]-dge[, xx[2]])[gmat[, x]]}, dge=dge, gmat=gmat))
  dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=dgraph1, gmat=gmat))
  dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
  bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Blue", "Grey", "Red"))(100),
              breaks=seq(-1, 1, length.out = 100), main=paste(method, "tPR"))}
dev.off()
#SupFig4
