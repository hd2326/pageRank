#' Marginal Probability Calculation
#'
#' Calculate marginal probability.
#'
#' @param expmat (matrix) Gene expression matrix.
#'
#' @param bin (matrix) Results of gene_bin function.
#'
#' @return (matrix) Marginal probability matrix.
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' 
#' @keywords internal

PX <- function(expmat, bin){
    pb <- txtProgressBar(min = 1, max = nrow(bin), style = 3)
    table <- do.call(rbind, lapply(seq_len(nrow(bin)), function(i, expmat,
                                                                bin, pb){
        setTxtProgressBar(pb, i)
        x <- expmat[rownames(bin)[i], ]
        b <- bin[i, ]
        table <- structure(rep(0, length(b)-1),
                           names=paste("bin", seq_len(length(b)-1), sep = "_"))
        for (p in seq_len(length(b)-1)) table[p] <- sum(x >= b[p] & x <= b[p+1])
        table/sum(table)}, expmat=expmat, bin=bin, pb=pb))
    rownames(table) <- rownames(bin)
    return(table)}
