#' Joint Probability Calculation
#'
#' Calculate joint probability.
#'
#' @param expmat (matrix) Gene expression matrix.
#'
#' @param bin (matrix) Results of gene_bin function.
#'
#' @param x (character) The first variable.
#'
#' @param y (character) The pairing second variable.
#'
#' @return (list) Joint probability.
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' 
#' @keywords internal

PXY <- function(expmat, bin, x, y){
    table <- list()
    pb <- txtProgressBar(min = 1, max = length(x), style = 3)
    for (i in seq_len(length(x))){
        setTxtProgressBar(pb, i)
        xx <- expmat[x[i], ]
        yy <- expmat[y[i], ]
        bx <- bin[x[i], ]
        by <- bin[y[i], ]
        t <- matrix(0, (length(bx)-1), (length(by)-1))
        for (p in seq_len(length(bx)-1))
            for (q in seq_len(length(by)-1))
                t[p, q] <- sum(xx>=bx[p]&xx<=bx[p+1]&yy>=by[q]&yy<=by[q+1])
        dimnames(t) <- list(paste(x[i], seq_len(length(bx)-1), sep = "_"),
                            paste(y[i], seq_len(length(by)-1), sep = "_"))
        table[[paste(x[i], y[i], sep = "_")]] <- t/sum(t)}
  return(table)}
