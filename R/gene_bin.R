#' Bin Gene Expression Space
#'
#' Bin gene expression space for marginal/joint probability calculation.
#'
#' @param genes (character) Genes to be analyzed.
#'
#' @param expmat (matrix) Gene expression matrix.
#'
#' @param sep (numeric) Number of bins.
#'
#' @return (matrix) Border values of gene expression bins.
#'
#' @author DING, HONGXU (hd2326@columbia.edu)

gene_bin <- function(genes, expmat, sep=5){
  expmat <- expmat[intersect(genes, rownames(expmat)), ]
  table <- t(apply(expmat, 1, function(x, sep) seq(min(x), max(x), length.out = (sep+1)), sep=sep))
  colnames(table) <- paste("bin", 1:(sep+1), sep = "_")
  return(table)}
