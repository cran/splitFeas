#' MM algorithm for nonlinear multiple-sets split feasibility problem
#' 
#' \code{nmsfp_mm} uses quasi-Newton updates to solve the nonlinear multiple-sets split feasibility problem.
#' 
#' @param x0 Initial iterate
#' @param f objective function
#' @param df gradient of objective function
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param plist1 list of projection functions for first set of constraints; each takes a single point and returns its projection
#' @param plist2 list of projection functions for second set of constraints; each takes a single point and returns its projection
#' @param h Function handle for output mapping
#' @param hgrad Handle for output mapping Jacobian
#' @param tol Stopping tolerance
#' @param max_iter Maximum number of iterations
#' @import compiler
#' @export
#' @seealso \code{mmqn_step}
nmsfp_mm <- cmpfun(function(x0, v, w, plist1, plist2, f, df, h, hgrad, tol=1e-10, max_iter=1e3) {
  x <- x0
  xhist <- matrix(NA,length(x0),max_iter+1)
  xhist[,1] <- x
  loss <- double(max_iter+1)
  loss[1] <- f(x)
  for (iter in 1:max_iter) {
    x <- mmqn_step(x, v, w, plist1, plist2, f, df, h, hgrad)
    xhist[,iter+1] <- x
    loss[iter+1] <- f(x)
    if (loss[iter+1] < tol) break
  }
  loss <- loss[1:(iter+1)]
  xhist <- xhist[,1:(iter+1),drop=FALSE]
  return(list(x=xhist,loss=loss))
})

#' MM algorithm (accelerated) for nonlinear multiple-sets split feasibility problem
#' 
#' \code{nmsfp_mmqn} uses quasi-Newton updates to solve the nonlinear multiple-sets split feasibility problem.
#' 
#' @param x0 Initial iterate
#' @param f objective function
#' @param df gradient of objective function
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param plist1 list of projection functions for first set of constraints; each takes a single point and returns its projection
#' @param plist2 list of projection functions for second set of constraints; each takes a single point and returns its projection
#' @param h Function handle for output mapping
#' @param hgrad Handle for output mapping Jacobian
#' @param qn number of secants
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @import compiler
#' @export
#' @seealso \code{mmqn_step}, \code{qnamm}
nmsfp_mmqn <- cmpfun(function(x0, v, w, plist1, plist2, f, df, h, hgrad, qn=5, tol=1e-10, max_iter=1e3) {
  x <- x0
  fx_mm <- function(x) {return(mmqn_step(x, v, w, plist1, plist2, f, df, h, hgrad))}
  sol <- qnamm(x, fx_mm, qn, f, max_iter=max_iter, tol=tol) 
  return(list(x=sol$Xhist,loss=sol$objective))
})
