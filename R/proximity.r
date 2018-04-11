#' Proximity function
#' 
#' \code{proximity} computes the proximity function.
#' 
#' @param x Current iterate
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param plist1 list of projection functions for first set of constraints; each takes a single point and returns its projection
#' @param plist2 list of projection functions for second set of constraints; each takes a single point and returns its projection
#' @param h Function handle for output mapping
#' @import compiler
#' @export
proximity <- cmpfun(function(x, v, w, plist1, plist2, h) {
  n <- length(plist1)
  p <- length(x)
  obj <- 0
  for (i in 1:n) {
    z <- as.matrix(x - plist1[[i]](x))
    obj <- obj + v[i]*norm(z,'f')**2
  }
  n <- length(plist2)
  for (i in 1:n) {
    hx <- h(x)
    z <- as.matrix(hx - plist2[[i]](hx))
    obj <- obj + w[i]*norm(z,'f')**2
  }
  return(0.5*obj)
}) 
