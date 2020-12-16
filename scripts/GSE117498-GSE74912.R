#this script is used for analyzing the human hematopoiesis process

library(viper)
source("functions.R")

graph <- lapply(c("HSC", "MPP", "CMP", "GMP", "MEP"), function(x){
  message(x)
  rpm <- get(load(paste("./GSE117498/", x, ".rda", sep="")))
  regul <- viper::pruneRegulon(get(load(paste("./GSE117498/", x, "-tf-regulon.rda", sep=""))), cutoff = 100)
  regul <- do.call(rbind, lapply(1:length(regul), function(i, regul) data.frame(reg=rep(names(regul)[i], length(regul[[i]]$tfmode)),
                                                                                target=names(regul[[i]]$tfmode),
                                                                                direction=as.integer(regul[[i]]$tfmode > 0), stringsAsFactors=F), regul=regul))
  regul <- regul[regul$reg %in% rownames(rpm) & regul$target %in% rownames(rpm), ]
  regul <- regul[regul$reg != regul$target, ]
  null <- P_null(expmat=rpm, net=regul, n=10000, sep=5, method="difference")
  P_graph(expmat=rpm, sep=5, net=regul, method="difference", null=null, threshold=1e-2)})
names(graph) <- c("HSC", "MPP", "CMP", "GMP", "MEP")
save(graph, file = "./GSE117498/graph.rda")
#expression-based networks

graph1 <- get(load("./GSE117498/graph.rda"))
graph2 <- get(load("./GSE74912/graph.rda"))
dgraph1 <- list("HSC->MPP"=diff_graph(graph1=graph1$MPP, graph2=graph1$HSC),
                "MPP->CMP"=diff_graph(graph1=graph1$CMP, graph2=graph1$MPP),
                "CMP->GMP"=diff_graph(graph1=graph1$GMP, graph2=graph1$CMP),
                "CMP->MEP"=diff_graph(graph1=graph1$MEP, graph2=graph1$CMP))
dgraph2 <- list("HSC->MPP"=diff_graph(graph1=graph2$MPP, graph2=graph2$HSC),
                "MPP->CMP"=diff_graph(graph1=graph2$CMP, graph2=graph2$MPP),
                "CMP->GMP"=diff_graph(graph1=graph2$GMP, graph2=graph2$CMP),
                "CMP->MEP"=diff_graph(graph1=graph2$MEP, graph2=graph2$CMP))
dge <- structure(lapply(names(graph1), function(x) rowMeans(get(load(paste("./GSE117498/", x, ".rda", sep=""))))), names=names(graph1))
pdf("SupFig5.pdf", width = 14, height = 9)
layout(matrix(1:8, 2, 4, byrow = T), widths = c(5, 4, 5, 4))
par(mar = c(6, 4, 4, 1))
for (method in c("additive", "multiplicative", "combined", "neutral")){
  pr <- structure(lapply(names(graph1), function(x, g1, g2) multiplex_page_rank(g1[[x]], g2[[x]], method=method), g1=graph1, g2=graph2), names=names(graph1))
  gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
  emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat) dge[[x]][gmat[, x]], dge=dge, gmat=gmat))
  dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=graph1, gmat=gmat))
  dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
  bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Grey", "Red"))(100),
              breaks=seq(0, 10, length.out = 100), main=paste(method, "PR"))
  
  pr <- structure(lapply(names(dgraph1), function(x, g1, g2) multiplex_page_rank(g1[[x]], g2[[x]], method=method), g1=dgraph1, g2=dgraph2), names=names(dgraph1))
  gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
  emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat){
    xx <- unlist(strsplit(x, split="->"))
    (dge[[xx[1]]]-dge[[xx[2]]])[gmat[, x]]}, dge=dge, gmat=gmat))
  dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=dgraph1, gmat=gmat))
  dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
  bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Blue", "Grey", "Red"))(100),
              breaks=seq(-3, 3, length.out = 100), main=paste(method, "tPR"))}
dev.off()
#SupFig5
