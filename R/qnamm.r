#' Quasi-Newton acceleration of MM algorithm
#' 
#' \code{qnamm} performs Quasi-Newton acceleration of an MM algorithm.
#' 
#' @param x initial iterate
#' @param fx_mm MM algorithm map
#' @param qn number of secants
#' @param fx_obj handle to objective function
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @param ... Additional arguments to pass to \code{fx_mm}
#' @import compiler corpcor
#' @export
#' @references H Zhou, D Alexander, and K Lange. (2011) A quasi-Newton acceleration method for high-dimensional optimization algorithms, Statistics and Computing, 21(2):261-273.
qnamm <- cmpfun(function(x, fx_mm, qn, fx_obj, max_iter=100, tol=1e-6, ...) {
  n <- length(x)
  U <- matrix(0,n,qn)
  V <- matrix(0,n,qn)
  objval <- Inf
  objective <- double(max_iter)
  Xhist <- matrix(NA,n,qn+max_iter)
  #
  #   accumulate the first QN differences for Quasi-Newton acceleration  
  #  
  for (i in 1:qn) {
    Xhist[,i] <- x
    x_old <- x    
    x <- fx_mm(x, ...)
    U[,i] <- x - x_old
  }
  V[,1:(qn-1)] <- U[,2:qn]
  x_old <- x    
  x <- fx_mm(x, ...)
  V[,qn] <- x - x_old
  old_secant <- 1
  C <- t(U)%*%(U-V)
  nacc <- 0
  nrej <- 0
  for (i in 1:max_iter) {
    Xhist[,qn+i] <- x
    objval_old <- objval
    x_old <- x
    x <- fx_mm(x, ...)      
    #
    #   do one more MM step to accumulate secant pairs  
    #
    U[,old_secant] <- x - x_old
    x_old <- x
    x <- fx_mm(x, ...)
    V[,old_secant] <- x - x_old
    C[old_secant,] <- t(U[,old_secant,drop=FALSE]) %*% (U-V)
    C[,old_secant] <- t(U) %*% (U[,old_secant,drop=FALSE] - V[,old_secant,drop=FALSE])
    new_secant <- old_secant
    old_secant <- (old_secant %% qn) + 1
    objval_MM <- fx_obj(x, ...)   
    #  
    #   quasi-Newton jump
    #     
    #      x_qn <- x_old + V %*% solve(C, t(U)%*%U[,new_secant,drop=FALSE])
    x_qn <- x_old + V %*% pseudoinverse(C) %*% (t(U)%*%U[,new_secant,drop=FALSE])
    x_qn <- fx_mm(x_qn, ...)
    objval_QN <- fx_obj(x_qn, ...)
    #
    #     choose MM vs QN jump
    #    
    if (objval_QN < objval_MM) {
      x <- x_qn;
      objval <- objval_QN;
      nacc <- nacc + 1
    } else {
      objval <- objval_MM;
      nrej <- nrej + 1
    }
    objective[i] <- objval
    #    
    # stopping rule
    #    
    print(norm(as.matrix(x-x_old),'f')/(norm(as.matrix(x_old),'f')+1))
    if (norm(as.matrix(x-x_old),'f') < tol*(norm(as.matrix(x_old),'f')+1)) break
  }
  print(paste("Accepted:", nacc))
  print(paste("Rejected:", nrej))
  return(list(x=x, objective=objective[1:i], iter=i+qn, Xhist=Xhist[,1:(i+qn),drop=FALSE]))
})
