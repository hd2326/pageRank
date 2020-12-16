#this script is used for analyzing the T-cell datasets

library(viper)
library(scatterplot3d)
source("functions.R")

load("./PBMC/T.rda")
regul <- viper::pruneRegulon(get(load("./PBMC/T-tf-regulon.rda")), cutoff = 100)
regul <- do.call(rbind, lapply(1:length(regul), function(i, regul) data.frame(reg=rep(names(regul)[i], length(regul[[i]]$tfmode)),
                                                                              target=names(regul[[i]]$tfmode),
                                                                              direction=as.integer(regul[[i]]$tfmode > 0), stringsAsFactors=F), regul=regul))
genes <- rownames(tpm)[rowSums(tpm > 0)/ncol(tpm) > 0.1]
regul <- regul[regul$reg %in% genes & regul$target %in% genes, ]
regul <- regul[regul$reg != regul$target, ]
null <- P_null(expmat=tpm, net=regul, n=10000, sep=5, method="difference")
graph <- P_graph(expmat=tpm, sep=5, net=regul, method="difference", null=null, threshold=1e-2)
save(graph, file = "./PBMC/graph.rda")
#expression-based networks

load("./PBMC/T.rda")
table <- get(load("./GSE101498/atac_table.rda"))
regulators <- read.table("./source/tf-homo-current-symbol.dat", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
regul <- atac_network(table=table, regulators=regulators, species="hg19", upstream=2000, downstream=200)
save(regul, file = "./GSE101498/atac_regul.rda")
regul <- regul[regul$reg %in% rownames(tpm) & regul$target %in% rownames(tpm), ]
regul <- regul[regul$reg != regul$target, ]
null <- P_null(expmat=tpm, net=regul, n=10000, sep=5, method="difference")
graph <- P_graph(expmat=tpm, sep=5, net=regul, method="difference", null=null, threshold=1e-3)
save(graph, file = "./GSE101498/atac_graph.rda")
#ATAC-Seq networks

load("./PBMC/T.rda")
table <- get(load("./GSE101498/hichip_table.rda"))
regulators <- read.table("./source/tf-homo-current-symbol.dat", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
regul <- conformation_network(table=table, regulators=regulators, species="hg19", upstream=2000, downstream=200)
save(regul, file = "./GSE101498/hichip_regul.rda")
regul <- regul[regul$reg %in% rownames(tpm) & regul$target %in% rownames(tpm), ]
regul <- regul[regul$reg != regul$target, ]
null <- P_null(expmat=tpm, net=regul, n=10000, sep=5, method="difference")
graph <- P_graph(expmat=tpm, sep=5, net=regul, method="difference", null=null, threshold=1e-3)
save(graph, file = "./GSE101498/hichip_graph.rda")
#HiChIP networks

library(scatterplot3d)
regulators <- read.table("../source/tf-homo-current-symbol.dat", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
graph1 <- get(load("../PBMC/graph.rda"))
graph2 <- get(load("../GSE101498/atac_graph.rda"))
graph3 <- get(load("../GSE101498/hichip_graph.rda"))
genes <- Reduce(union, list(names(V(graph1)), names(V(graph2)), names(V(graph3))))
pr1 <- pr2 <- pr3 <- structure(rep(0, length(genes)), names=genes)
pr1[names(V(graph1))] <- V(graph1)$pagerank
pr2[names(V(graph2))] <- V(graph2)$pagerank
pr3[names(V(graph3))] <- V(graph3)$pagerank
pr <- sort(multiplex_page_rank(graph1, graph2, graph3, method="combined"), decreasing = T)
pdf("SupFig6.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
cex <- structure(rep(15*min(pr)+0.5, length(genes)), names=genes)
cex[names(pr)] <- 15*pr+0.5
plot <- scatterplot3d(pr1, pr2, pr3, xlab="scRNA-Seq", ylab="ATAC-Seq", zlab="HiChIP",
                      color=(genes%in%regulators)+1, cex.symbols=cex, cex.lab=1.5, main="Combined PR", cex.main=2)
genes <- names(pr)[pr >= quantile(pr, 1-20/length(pr))]
text(plot$xyz.convert(cbind(pr1, pr2, pr3)[genes, ]), labels=genes, pos=3, col=4)
legend("topleft", legend = c("TF", "Target"), fill = 2:1, cex = 1.5, bty = "n", border = NA)
legend("bottomright", legend=c("Size~PR"), cex=1.5, bty="n", border=NA)

par(mar = c(8, 2, 2, 2))
table <- rbind(pr1[names(pr)[1:20]]/max(pr1), pr2[names(pr)[1:20]]/max(pr2), pr3[names(pr)[1:20]]/max(pr3))
barplot(apply(table, 2, prop.table), yaxt = "n", ylim = c(0, 1.5), col = c("Red", "Green", "Blue"), names.arg = names(pr)[1:20], las = 2)
axis(side = 1, line = 5, at = 12, labels = "Rank", cex.axis = 2, tick = F)
legend("topright", legend=c("scRNA-Seq", "ATAC-Seq", "HiChIP"), fill=c("Red", "Green", "Blue"), cex=1.5, bty="n", border=NA)
dev.off()
#SupFig6

graph1 <- get(load("./PBMC/graph.rda"))
graph2 <- get(load("./GSE101498/atac_graph.rda"))
graph3 <- get(load("./GSE101498/hichip_graph.rda"))
genes <- Reduce(union, list(names(V(graph1)), names(V(graph2)), names(V(graph3))))
pr <- list()
for (x in c("additive", "multiplicative", "combined", "neutral")){
  pr1 <- pr2 <- pr3 <- pr4 <- pr5 <- pr6 <- structure(rep(0, length(genes)), names=genes)
  pr1[names(V(graph1))] <- multiplex_page_rank(graph1, graph2, graph3, method=x)
  pr2[names(V(graph2))] <- multiplex_page_rank(graph2, graph3, graph1, method=x)
  pr3[names(V(graph3))] <- multiplex_page_rank(graph3, graph1, graph2, method=x)
  pr4[names(V(graph1))] <- V(graph1)$pagerank
  pr5[names(V(graph2))] <- V(graph2)$pagerank
  pr6[names(V(graph3))] <- V(graph3)$pagerank
  pr[[x]] <- list(pr1, pr2, pr3, pr4, pr5, pr6)}
pdf("SupFig7.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
scatterplot3d(unlist(lapply(pr, function(x) x[[1]]))^0.5, 
              unlist(lapply(pr, function(x) x[[2]]))^0.5,
              unlist(lapply(pr, function(x) x[[3]]))^0.5,
              color=unlist(lapply(2:5, function(i, genes) rep(i, length(genes)), genes=genes)),
              xlab="scRNA-Seq, sqrt(PR)", ylab="ATAC-Seq, sqrt(PR)", zlab="HiChIP, sqrt(PR)",
              cex.lab=1.5, main="Base Network", cex.main=2)
legend("topleft", legend=c("additive", "multiplicative", "combined", "neutral"), fill=2:5, bty="n", border = NA)

par(mar = c(5, 5, 3, 2))
plot(NA, NA, xlim=range(pr)^0.5, ylim=range(pr)^0.5, xlab="regular, sqrt(PR)", ylab="multiplex, sqrt(PR)", cex.lab=1.5, main="PageRank Method", cex.main=2)
for (i in 1:4) for (j in 1:3) points(pr[[i]][[j+3]]^0.5, pr[[i]][[j]]^0.5, col=i+1, pch=j)
abline(0, 1, NULL, NULL, col = 1, lty = 2)
legend("bottomright", legend = c("scRNA-Seq", "ATAC-Seq", "HiChIP"), pch = 1:3, bty = "n", border = NA)
dev.off()
#SupFig7
