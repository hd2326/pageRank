#' Generate Timewise Average Gene Expression
#'
#' Generate timewise average gene expression.
#'
#' @param time (character) Time-annotation of samples.
#'
#' @param expmat (matrix) Gene expression matrix.
#'
#' @return (matrix) Time-wise average gene expression.
#'
#' @examples
#' expmat <- matrix(rnorm(90), 10, 9, dimnames=list(LETTERS[1:10], 1:9))
#' time <- c(rep("T1", 3), rep("T2", 3), rep("T3", 3))
#' time_expmat(time, expmat)
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @export

time_expmat <- function(time, expmat){
    table <- lapply(unique(time), function(t, time, expmat)
        rowMeans(expmat[, time == t]), time=time, expmat=expmat)
    table <- structure(do.call(cbind, table),
                       dimnames=list(rownames(expmat), unique(time)))
    return(table)}
