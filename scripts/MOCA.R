#this script is used for analyzing the MOCA datasets

library(viper)
source("functions.R")

files <- list.files(pattern = "regulon.rda")
tissues <- unique(unlist(lapply(strsplit(files, split="_"), function(x) x[1])))
for (tissue in tissues){
  graph <- lapply(files[grepl(tissue, files)], function(f){
    message(f)
    tpm <- get(load(gsub("-tf-regulon.rda", ".rda", f)))
    regul <- viper::pruneRegulon(get(load(f)), cutoff = 100)
    regul <- do.call(rbind, lapply(1:length(regul), function(i, regul) data.frame(reg=rep(names(regul)[i], length(regul[[i]]$tfmode)),
                                                                                  target=names(regul[[i]]$tfmode),
                                                                                  direction=as.integer(regul[[i]]$tfmode > 0), stringsAsFactors=F), regul=regul))
    regul <- regul[regul$reg %in% rownames(tpm) & regul$target %in% rownames(tpm), ]
    regul <- regul[regul$reg != regul$target, ]
    null <- P_null(expmat=tpm, net=regul, n=10000, sep=5, method="difference")
    P_graph(expmat=tpm, sep=5, net=regul, method="difference", null=null, threshold=1e-2)})
  names(graph) <- gsub("-tf-regulon.rda", "", files[grepl(tissue, files)])
  save(graph, file = paste(tissue, "_graph.rda", sep = ""))}
#expression-based networks

library(viper)
files <- list.files(pattern = "graph.rda")
nes <- lapply(files, function(f){
  message(f)
  graph <- get(load(f))
  if (sum(grepl("E9.5", names(graph))) == 1) graph <- graph[c(length(graph), 1:(length(graph)-1))]
  names(graph) <- unlist(lapply(strsplit(names(graph), split = "_"), function(x) x[2]))
  
  dgraph <- lapply(2:length(graph), function(i, graph) diff_graph(graph1=graph[[i]], graph2=graph[[i-1]]), graph=graph)
  names(dgraph) <- unlist(lapply(1:(length(graph)-1), function(i, stage) paste(stage[i:(i+1)], collapse="->"), stage=names(graph)))
  tpr <- lapply(dgraph, function(g) structure(V(g)$pagerank, names=names(V(g))))
  
  dge <- do.call(cbind, lapply(names(graph), function(x, f) rowMeans(get(load(gsub("graph", x, f)))), f=f))
  dge <- do.call(cbind, lapply(1:(ncol(dge)-1), function(i, dge) dge[, i+1]-dge[, i], dge=dge))
  colnames(dge) <- names(dgraph)
  
  dregul <- lapply(dgraph, function(g) graph_to_regulon(g))
  dnes <- structure(lapply(1:length(dregul), function(i, regul, mat){
    nes <- aREA(mat[, i], regul[[i]], minsize=5)$nes
    structure(nes[, 1], names=rownames(nes))}, regul=dregul, mat=abs(dge)), names=names(dregul))
  
  structure(lapply(names(dgraph), function(x, tpr, dnes){
    list(tPR=tpr[[x]][intersect(names(tpr[[x]]), names(dnes[[x]]))],
         nES=dnes[[x]][intersect(names(tpr[[x]]), names(dnes[[x]]))])}, tpr=tpr, dnes=dnes), names=names(dgraph))})
names(nes) <- gsub("_graph.rda", "", files)
save(nes, file = "MOCA_nes.rda")
#VIPER analysis

pdf("SupFig2.pdf", width = 40, height = 28)
layout(matrix(1:70, 7, 10, byrow = T), widths = rep(c(5, 4), 5))
par(mar = c(8, 4, 4, 1))
for (f in list.files(pattern = "graph.rda")){
  message(f)
  graph <- get(load(f))
  if (sum(grepl("E9.5", names(graph))) == 1) graph <- graph[c(length(graph), 1:(length(graph)-1))]
  names(graph) <- unlist(lapply(strsplit(names(graph), split = "_"), function(x) x[2]))
  dgraph <- lapply(2:length(graph), function(i, graph) diff_graph(graph1=graph[[i]], graph2=graph[[i-1]]), graph=graph)
  names(dgraph) <- unlist(lapply(1:(length(graph)-1), function(i, stage) paste(stage[i:(i+1)], collapse="->"), stage=names(graph)))
  
  files <- list.files(pattern=gsub("_graph.rda", "", f))
  dge <- do.call(cbind, lapply(files[seq(2, length(files), 2)], function(x) rowMeans(get(load(x)))))
  colnames(dge) <- unlist(lapply(strsplit(gsub(".rda", "", files[seq(2, length(files), 2)]), split="_"), function(x) x[2]))
  dge <- dge[, names(graph)]
  
  pr <- structure(lapply(graph, function(g) structure(V(g)$pagerank, names=names(V(g)))), names=names(graph))
  gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
  emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat) dge[match(gmat[, x], rownames(dge)), x], dge=dge, gmat=gmat))
  emat[is.na(emat)] <- 0
  dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat){
    d <- igraph::degree(graph[[x]], v=gmat[, x][!is.na(gmat[, x])])
    c(d, rep(-0.5, sum(is.na(gmat[, x]))))}, graph=graph, gmat=gmat))
  dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
  bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Grey", "Red"))(100), breaks=seq(0, 8, length.out=100), main=gsub("_graph.rda", "", f))
  
  pr <- structure(lapply(dgraph, function(g) structure(V(g)$pagerank, names=names(V(g)))), names=names(dgraph))
  gmat <- do.call(cbind, lapply(pr, function(x) names(sort(x, decreasing = T))[1:20]))
  emat <- do.call(cbind, lapply(names(pr), function(x, dge, gmat){
    xx <- unlist(strsplit(x, split="->"))
    (dge[, xx[1]]-dge[, xx[2]])[gmat[, x]]}, dge=dge, gmat=gmat))
  dmat <- do.call(cbind, lapply(names(pr), function(x, graph, gmat) igraph::degree(graph[[x]], v=gmat[, x]), graph=dgraph, gmat=gmat))
  dimnames(gmat) <- dimnames(emat) <- dimnames(dmat) <- list(paste("#", 1:20, sep = ""), names(pr))
  bubble_plot(s_mat=dmat, c_mat=emat, n_mat=gmat, col=colorRampPalette(c("Blue", "Grey", "Red"))(100), breaks=seq(-1, 1, length.out=100))}
dev.off()
#SupFig2

pdf("SupFig11S.pdf", width = 35, height = 25)
par(mfrow = c(5, 7), mar = c(5, 5, 3, 1))
for (i in 1:length(nes)){
  plot(NA, NA, xlim = range(lapply(nes[[i]], function(x) x$nES)), ylim = range(lapply(nes[[i]], function(x) x$tPR)),
       xlab = "nES", ylab = "temporal PR", cex.axis = 2, cex.lab = 2, main = names(nes)[i], cex.main = 2)
  for (j in 1:length(nes[[i]])){
    lines(nes[[i]][[j]]$nES, nes[[i]][[j]]$tPR, type = "p", pch = 16, cex = 2, col = j)
    names <- names(sort(nes[[i]][[j]]$tPR, decreasing = T))[1:10]
    text(nes[[i]][[j]]$nES[names], nes[[i]][[j]]$tPR[names], labels = names)}
  r <- unlist(lapply(nes[[i]], function(x) cor(x$nES, x$tPR)))
  legend("topleft", legend = paste(names(r), ", r=", signif(r, 2), sep=""), cex = 2, fill = 1:8, bty = "n", border = NA)}
dev.off()
#SupFig11
