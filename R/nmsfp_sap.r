#' Self-adaptive projection-type method algorithm for nonlinear multiple-sets split feasibility problem
#' 
#' \code{nmsfp_sap} performs the self-adaptive projection-type method of Li et al.
#' 
#' @param x0 Initial iterate
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param plist1 list of projection functions for first set of constraints; each takes a single point and returns its projection
#' @param plist2 list of projection functions for second set of constraints; each takes a single point and returns its projection
#' @param proj handle to projection operation.
#' @param f objective function
#' @param df gradient of objective function
#' @param h Function handle for output mapping
#' @param hgrad Handle for output mapping Jacobian
#' @param delta step-size parameter
#' @param mu step-size parameter
#' @param beta0 initial 
#' @param tol Tolerance
#' @param max_iter Maximum number of iterations
#' @export
nmsfp_sap <- cmpfun(function(x0,v,w,plist1,plist2,proj,f,df,h,hgrad,delta=0.3,mu=0.7,beta0=1,tol=1e-10,max_iter=1e3) {
  x <- x0
  xhist <- matrix(NA,2,max_iter+1)
  xhist[,1] <- x
  loss <- double(max_iter+1)
  loss[1] <- f(x)
  tau <- double(max_iter)
  tau[1:50] <- 0.5
  for (iter in 1:max_iter) {
    gamma <- beta0
    sol <- nmsfp_sap_one_step(x,delta,mu,tau[iter],gamma,df,proj)
    x <- sol$x
    gamma <- sol$gamma
    xhist[,iter+1] <- x
    loss[iter+1] <- f(x)
    if (loss[iter+1] < tol) break
  }
  loss <- loss[1:(iter+1)]
  xhist <- xhist[,1:(iter+1),drop=FALSE]
  return(list(x=xhist,loss=loss))
})
