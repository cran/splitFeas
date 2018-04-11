#' Compute soft-max
#' 
#' \code{softmax} returns the soft maximum of a collection of reals.
#' 
#' @param x input
#' @param a scaling factor
#' @import matrixStats
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10
#' x <- rnorm(n)
#' softmax(x)
softmax <- function(x, a=100){
  return(logSumExp(x*a)/a)
}

