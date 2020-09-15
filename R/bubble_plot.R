#' Make Bubbleplot
#'
#' Make bubbleplot.
#'
#' @param s_mat (matrix) Matrix denotes the size of bubbles.
#'
#' @param c_mat (matrix) Matrix denotes the color of bubbles.
#'
#' @param n_mat (matrix) Matrix denotes the name of bubbles.
#'
#' @param col (character) Colors.
#'
#' @param breaks (numeric) Breakpoints of colors.
#'
#' @param main (character) Title.
#' 
#' @return (NULL) A bubbloplot.
#'
#' @examples
#' s_mat <- c_mat <- n_mat <- matrix(1:12, 3, 4, dimnames=list(1:3, 1:4))
#' bubble_plot(s_mat, c_mat, n_mat, breaks=seq(1, 12, length.out=100), main="")
#' 
#' @import graphics
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @export

bubble_plot <- function(s_mat, c_mat, n_mat,
                        col=colorRampPalette(c("Blue", "Grey", "Red"))(100),
                        breaks=seq(-2, 2, length.out=100), main=NULL){
    row_names <- Reduce(intersect, list(rownames(s_mat),
                                        rownames(c_mat),
                                        rownames(n_mat)))
    col_names <- Reduce(intersect, list(colnames(s_mat),
                                        colnames(c_mat),
                                        colnames(n_mat)))
    s_mat <- s_mat[row_names, col_names]
    c_mat <- c_mat[row_names, col_names]
    n_mat <- n_mat[row_names, col_names]
    s_mat <- 2*t(t(s_mat)/apply(s_mat, 2, max))+0.5
    c_mat <- structure(t(apply(c_mat, 1, function(x) get_color_gradient(x, col=col, breaks=breaks))), dimnames=dimnames(c_mat))
    plot(NA, NA, xlim = c(0.5, length(col_names)+1), ylim = c(0.5, length(row_names)+0.5),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = main, cex.main = 2)
    for (i in seq_len(ncol(s_mat))) points(rep(i, nrow(s_mat)),
                                           seq_len(nrow(s_mat)),
                                           cex = s_mat[, i],
                                           col = c_mat[, i], pch = 16)
    for (i in seq_len(ncol(n_mat))) text(rep(i+0.5, nrow(n_mat)),
                                         seq_len(nrow(n_mat)),
                                         labels = n_mat[, i])
    axis(side = 1, at = seq_len(length(col_names)), labels = col_names,
         tick = FALSE, las = 2)
    axis(side = 2, at = seq_len(length(row_names)), labels = row_names,
         tick = FALSE, las = 2)}
