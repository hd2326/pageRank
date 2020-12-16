#this script contains all necessary functions

library(igraph)
library(annotate)
library(GenomicFeatures)
library(GenomicRanges)
library(TFBSTools)
library(JASPAR2018)
library(motifmatchr)

gene_bin <- function(genes, expmat, sep=5){
  expmat <- expmat[intersect(genes, rownames(expmat)), ]
  table <- t(apply(expmat, 1, function(x, sep) seq(min(x), max(x), length.out = (sep+1)), sep=sep))
  colnames(table) <- paste("bin", 1:(sep+1), sep = "_")
  return(table)}

PX <- function(expmat, bin){
  pb <- txtProgressBar(min = 1, max = nrow(bin), style = 3)
  table <- do.call(rbind, lapply(1:nrow(bin), function(i, expmat, bin, pb){
    setTxtProgressBar(pb, i)
    x <- expmat[rownames(bin)[i], ]
    b <- bin[i, ]
    table <- structure(rep(0, length(b)-1), names=paste("bin", 1:(length(b)-1), sep = "_"))
    for (p in 1:length(b)-1) table[p] <- sum(x >= b[p] & x <= b[p+1])
    table/sum(table)}, expmat=expmat, bin=bin, pb=pb))
  rownames(table) <- rownames(bin)
  return(table)}

PXY <- function(expmat, bin, x, y){
  table <- list()
  pb <- txtProgressBar(min = 1, max = length(x), style = 3)
  for (i in 1:length(x)){
    setTxtProgressBar(pb, i)
    xx <- expmat[x[i], ]
    yy <- expmat[y[i], ]
    bx <- bin[x[i], ]
    by <- bin[y[i], ]
    t <- matrix(0, (length(bx)-1), (length(by)-1))
    for (p in 1:(length(bx)-1)) for (q in 1:(length(by)-1)) t[p, q] <- sum(xx>=bx[p] & xx<=bx[p+1] & yy>=by[q] & yy<=by[q+1])
    dimnames(t) <- list(paste(x[i], 1:(length(bx)-1), sep = "_"), paste(y[i], 1:(length(by)-1), sep = "_"))
    table[[paste(x[i], y[i], sep = "_")]] <- t/sum(t)}
  return(table)}

PXPY <- function(px, combinations){
  pb <- txtProgressBar(min = 1, max = length(combinations), style = 3)
  table <- lapply(1:length(combinations), function(i, combinations, px, pb){
    setTxtProgressBar(pb, i)
    x <- unlist(lapply(strsplit(combinations[i], split = "_"), function(x) x[1]))
    y <- unlist(lapply(strsplit(combinations[i], split = "_"), function(x) x[2]))
    pxx <- px[x, ]
    pxy <- px[y, ]
    table <- outer(pxx, pxy, "*")
    dimnames(table) <- list(paste(x, 1:length(pxx), sep = "_"), paste(y, 1:length(pxy), sep = "_"))
    table}, combinations=combinations, px=px, pb=pb)
  names(table) <- combinations
  return(table)}  
    
P_dist <- function(pxy, pxpy, method=c("difference", "mi")){
  pb <- txtProgressBar(min = 1, max = length(pxy), style = 3)
  dist <- lapply(1:length(pxy), function(i, pxy, pxpy, method, pb){
    setTxtProgressBar(pb, i)
    if (method == "difference") dist <- pxy[[i]]-pxpy[[i]]
    else if (method == "mi") dist <- pxy[[i]]*log2(pxy[[i]]/pxpy[[i]])
    else (message("method should be difference or mi"))
    dist[dist < 0 | !is.finite(dist)] <- 0
    sum(dist)}, pxy=pxy, pxpy=pxpy, method=method, pb=pb)
  names(dist) <- names(pxy)
  return(dist)} 

P_null <- function(expmat, net, n=10000, sep=5, method=c("difference", "mi")){
  net <- data.frame(reg=sample(unique(net$reg), size=n, replace=T),
                    target=sample(unique(net$target), size=n, replace=T), stringsAsFactors=F)
  bin <- gene_bin(genes=union(net$reg, net$target), expmat=expmat, sep=sep)
  px <- PX(expmat=expmat, bin=bin)
  pxy <- PXY(expmat=expmat, bin=bin, x=net$reg, y=net$target)
  pxpy <- PXPY(px=px, combinations=names(pxy))
  dist <- P_dist(pxy=pxy, pxpy=pxpy, method=method)
  null <- ecdf(unlist(dist))
  return(null)}

P_graph <- function(expmat, net, sep=5, method=c("difference", "mi"), null=NULL, threshold=1e-3){
  bin <- gene_bin(genes=union(net$reg, net$target), expmat=expmat, sep=sep)
  px <- PX(expmat=expmat, bin=bin)
  pxy <- PXY(expmat=expmat, bin=bin, x=net$reg, y=net$target)
  pxpy <- PXPY(px=px, combinations=names(pxy))
  dist <- P_dist(pxy=pxy, pxpy=pxpy, method=method)
  if (is.null(null)) null <- ecdf(unlist(dist))
  pvalue <- structure(1-null(unlist(dist)), names=names(dist))
  graph <- do.call(rbind, lapply(names(pvalue)[pvalue <= threshold], function(x) unlist(strsplit(x, split = "_"))))
  graph <- igraph::graph_from_data_frame(graph[, 2:1], directed=T) %>%
    set_edge_attr(name="pvalue", value=pvalue[pvalue <= threshold]) %>%
    set_edge_attr(name="direction", value=net$direction[pvalue <= threshold])
  graph <- set_vertex_attr(graph=graph, name="pagerank", value=igraph::page_rank(graph)$vector)
  return(graph)}

time_expmat <- function(time, expmat){
  table <- lapply(unique(time), function(t, time, expmat) rowMeans(expmat[, time == t]), time=time, expmat=expmat)
  table <- structure(do.call(cbind, table), dimnames=list(rownames(expmat), unique(time)))
  return(table)}

atac_network <- function(table, regulators, species=c("hg19", "mm9", "hg38", "mm10"), upstream=2000, downstream=200){
  table <- makeGRangesFromDataFrame(table)
  #convert peak table to GRanges
  if (species == "hg19"){
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    promoter <- promoters(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), upstream=upstream, downstream=downstream)
    names(promoter) <- getSYMBOL(names(promoter), data="org.Hs.eg")
  }else if (species == "mm9"){
    library(org.Mm.eg.db)
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    promoter <- promoters(genes(TxDb.Mmusculus.UCSC.mm9.knownGene), upstream=upstream, downstream=downstream)
    names(promoter) <- getSYMBOL(names(promoter), data="org.Mm.eg")
  }else if (species == "hg38"){
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    promoter <- promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), upstream=upstream, downstream=downstream)
    names(promoter) <- getSYMBOL(names(promoter), data="org.Hs.eg")
  }else if (species == "mm10"){
    library(org.Mm.eg.db)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    promoter <- promoters(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), upstream=upstream, downstream=downstream)
    names(promoter) <- getSYMBOL(names(promoter), data="org.Mm.eg")
  }else message("only hg19, hg38, mm9 and mm10 are supported")
  promoter <- promoter[!is.na(names(promoter))]
  #get promoter region of genes
  target <- as.data.frame(findOverlaps(table, promoter))
  target$subjectHits <- names(promoter)[target$subjectHits]
  #map between peaks and promoters by GRanges overlapping searching
  id <- unique(names(table)[target$queryHits])
  #peaks can be associated with promoters
  if (species == "hg19"){
    library(BSgenome.Hsapiens.UCSC.hg19)
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Homo sapiens"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg <- matchMotifs(PFMatrixList, table[id], genome="BSgenome.Hsapiens.UCSC.hg19")
  }else if (species == "mm9"){
    library(BSgenome.Mmusculus.UCSC.mm9)
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Mus musculus"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg <- matchMotifs(PFMatrixList, table[id], genome="BSgenome.Mmusculus.UCSC.mm9")
  }else if (species == "hg38"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Homo sapiens"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg <- matchMotifs(PFMatrixList, table[id], genome="BSgenome.Hsapiens.UCSC.hg38")
  }else if (species == "mm10"){
    library(BSgenome.Mmusculus.UCSC.mm10)
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Mus musculus"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg <- matchMotifs(PFMatrixList, table[id], genome="BSgenome.Mmusculus.UCSC.mm10")}
  reg <- as.matrix(motifMatches(reg))
  colnames(reg) <- structure(unlist(lapply(PFMatrixList, function(x) x@name)), names=NULL)
  #map between peaks and regulators by motif searching
  target <- lapply(unique(target$queryHits), function(i, target) target$subjectHits[target$queryHits == i], target=target)
  reg <- apply(reg, 1, function(x) names(x)[x])
  names(target) <- names(reg) <- id
  regul <- lapply(id, function(x, target, reg) expand.grid(target=target[[x]], reg=reg[[x]], stringsAsFactors=F), target=target, reg=reg)
  regul <- do.call(rbind, regul)
  #regulator-target network
  return(regul)}

diff_graph <- function(graph1, graph2){
  graph1 <- as.matrix(igraph::as_adjacency_matrix(graph1))
  graph2 <- as.matrix(igraph::as_adjacency_matrix(graph2))
  vertices <- union(rownames(graph1), rownames(graph2))
  g1 <- g2 <- matrix(0, length(vertices), length(vertices), dimnames = list(vertices, vertices))
  g1[rownames(graph1), colnames(graph1)] <- graph1
  g2[rownames(graph2), colnames(graph2)] <- graph2
  graph <- g1 - g2
  graph <- do.call(rbind, lapply(1:nrow(graph), function(i, graph){
    data.frame(target=rep(rownames(graph)[i], sum(graph[i, ] != 0)),
               reg=colnames(graph)[graph[i, ] != 0],
               moi=graph[i, graph[i, ] != 0])
  }, graph=graph))
  graph <- igraph::graph_from_data_frame(graph, directed=T) %>%
    set_edge_attr(name="moi", value=graph$moi)
  graph <- set_vertex_attr(graph=graph, name="pagerank", value=igraph::page_rank(graph)$vector)
  return(graph)}

clean_graph <- function(graph, size=NULL, vertices=NULL, pagerank=NULL){
  v <- v1 <- v2 <- names(V(graph))
  if (!is.null(vertices)) v1 <- intersect(v, vertices)
  if (!is.null(pagerank)) v2 <- v[V(graph)$pagerank >= pagerank]
  graph <- igraph::delete_vertices(graph, setdiff(v, intersect(v1, v2)))
  if (!is.null(size)){
    cluster <- igraph::clusters(graph)
    v <- names(V(graph))[cluster$membership %in% which(cluster$csize < size)]
    graph <- igraph::delete_vertices(graph, v)}
  graph <- set_vertex_attr(graph=graph, name="pagerank", value=igraph::page_rank(graph)$vector)
  return(graph)}

adjust_graph <- function(graph, damping=0.85, personalized=NULL, weights=NULL){
  pagerank <- igraph::page_rank(graph, damping=damping, personalized=personalized, weights=weights)$vector
  graph <- set_vertex_attr(graph=graph, name="pagerank", value=pagerank)
  return(graph)}

multiplex_page_rank <- function (graph, ..., method=c("additive", "multiplicative", "combined", "neutral"), damping=0.85){
  pagerank1 <- structure(V(graph)$pagerank, names=names(V(graph)))
  pagerank2 <- lapply(list(...), function(g) structure(V(g)$pagerank, names=names(V(g))))
  pagerank2 <- do.call(rbind, lapply(pagerank2, function(x, names){
    pr <- structure(rep(0, length(names)), names=names)
    pr[intersect(names, names(x))] <- x[intersect(names, names(x))]
    pr}, names=intersect(names(pagerank1), Reduce(union, lapply(pagerank2, function(x) names(x))))))
  #pagerank values
  personalized <- structure(rep(0, length(pagerank1)), names=names(pagerank1))
  personalized[colnames(pagerank2)] <- apply(pagerank2, 2, prod)
  personalized[personalized == 0] <- min(personalized[personalized > 0])
  weights <- unlist(lapply(igraph::as_data_frame(graph)$to, function(x, pagerank2){
    if (x %in% colnames(pagerank2)) prod(pagerank2[, x])
    else 0}, pagerank2=pagerank2))
  weights[weights == 0] <- min(weights[weights > 0])
  #personalized and weights for multiplex pagerank as in halu et al.
  if (method == "additive") pagerank <- igraph::page_rank(graph, directed=T, damping=damping, personalized=personalized, weights=NULL)$vector
  else if (method == "multiplicative") pagerank <- igraph::page_rank(graph, directed=T, damping=damping, personalized=NULL, weights=weights)$vector
  else if (method == "combined") pagerank <- igraph::page_rank(graph, directed=T, damping=damping, personalized=personalized, weights=weights)$vector
  else if (method == "neutral") pagerank <- igraph::page_rank(graph, directed=T, damping=damping, personalized=NULL, weights=NULL)$vector
  else message("method should be additive, multiplicative, combined or neutral")
  return(pagerank)}

get_color_gradient <- function (x, col=colorRampPalette(c("Blue", "Grey", "Red"))(100), breaks=seq(-2, 2, length.out=100)){
  color <- NULL
  if (length(col)  == length(breaks)) color <- col[unlist(lapply(x, function(xx, breaks) which.min(abs(xx - breaks)), breaks = breaks))]
  else message("col and breaks should be the same length")
  return(color)}

bubble_plot <- function(s_mat, c_mat, n_mat, col=colorRampPalette(c("Blue", "Grey", "Red"))(100),
                        breaks=seq(-2, 2, length.out=100), main=NULL){
  row_names <- Reduce(intersect, list(rownames(s_mat), rownames(c_mat), rownames(n_mat)))
  col_names <- Reduce(intersect, list(colnames(s_mat), colnames(c_mat), colnames(n_mat)))
  s_mat <- s_mat[row_names, col_names]
  c_mat <- c_mat[row_names, col_names]
  n_mat <- n_mat[row_names, col_names]
  s_mat <- 2*t(t(s_mat)/apply(s_mat, 2, max))+0.5
  c_mat <- structure(t(apply(c_mat, 1, function(x) get_color_gradient(x, col=col, breaks=breaks))), dimnames=dimnames(c_mat))
  plot(NA, NA, xlim = c(0.5, length(col_names)+1), ylim = c(0.5, length(row_names)+0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = main, cex.main = 2)
  for (i in 1:ncol(s_mat)) points(rep(i, nrow(s_mat)), 1:nrow(s_mat), cex = s_mat[, i], col = c_mat[, i], pch = 16)
  for (i in 1:ncol(n_mat)) text(rep(i+0.5, nrow(n_mat)), 1:nrow(n_mat), labels = n_mat[, i])
  axis(side = 1, at = 1:length(col_names), labels = col_names, tick = F, las = 2)
  axis(side = 2, at = 1:length(row_names), labels = row_names, tick = F, las = 2)}

graph_to_regulon <- function(graph){
  graph <- igraph::as_data_frame(graph)
  regul <- structure(lapply(unique(graph$to), function(x, graph){
    target <- graph$from[graph$to == x]
    list(tfmode=structure(rep(1, length(target)), names=target), likelihood=rep(1, length(target)))
  }, graph=graph), names=unique(graph$to))
  return(regul)}
