#' MM-quasi-Newton step
#' 
#' \code{mmqn_step} computes a single step.
#' 
#' @param x Current iterate
#' @param f objective function
#' @param df gradient of objective function
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param plist1 list of projection functions for first set of constraints; each takes a single point and returns its projection
#' @param plist2 list of projection functions for second set of constraints; each takes a single point and returns its projection
#' @param h Function handle for output mapping
#' @param hgrad Handle for output mapping Jacobian
#' @param woodbury Boolean: TRUE to use the Woodbury inversion formula
#' @import compiler
#' @export
mmqn_step <- cmpfun(function(x, v, w, plist1, plist2, f, df, h, hgrad, woodbury=TRUE) {
  df <- dg(x, x, v, w, plist1, plist2, h, hgrad)
  if (woodbury==TRUE) {
    dx <- -wood_inv_solve(x, v, w, hgrad, df)
  } else {
    H <- ddg(x, v, w, hgrad)
    dx <- -solve(H, df)
  }
  t <- backtrack(x, dx, f, df)
  return(x + t*dx)
})

#' Backtracking Line Search
#' 
#' Given a descent direction \code{backtrack} computes a step size that ensures sufficient decrease in an objective.
#' 
#' @param x Current iterate
#' @param dx Descent direction
#' @param f objective function
#' @param df gradient of objective function
#' @param alpha sufficient decrease parameter
#' @param beta sufficient decrease parameter
#' @import compiler
#' @export
backtrack <- cmpfun(function(x, dx, f, df, alpha=0.01, beta=0.8) {
  t <- 1
  g <- df(x)
  u <- alpha*sum(g*dx)
  k <- 1
  repeat {
    if (f(x + t*dx) <= f(x) + t*u) break
    t <- beta*t
    print(paste0("backtrack ",k))
    k <- k + 1
  }
  return(t)
})

#' Compute the approximate Hessian of the majorization.
#' 
#' \code{ddg} computes the Hessian of the majorization of the proximity function.
#' 
#' @param x non-anchor point
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param hgrad Handle for output mapping Jacobian
#' @import compiler
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10
#' p <- 2
#' x <- matrix(rnorm(p),p,1)
#' v <- 1
#' w <- 1
#' A <- matrix(rnorm(n*p),n,p)
#' hgrad <- function(x) {return(t(A))}
#' sol <- ddg(x,v,w,hgrad)
ddg <- cmpfun(function(x, v, w, hgrad) {
  p <- length(x)
  dh <- as.matrix(hgrad(x))
  return(diag(sum(v),p,p) + sum(w)*dh%*%t(dh))
})

#' Compute the gradient of the majorization.
#' 
#' \code{dg} computes the gradient of the majorization of the proximity function.
#' 
#' @param x non-anchor point
#' @param xa Anchor point
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param plist1 list of projection functions for first set of constraints; each takes a single point and returns its projection
#' @param plist2 list of projection functions for second set of constraints; each takes a single point and returns its projection
#' @param h Function handle for output mapping
#' @param hgrad Handle for output mapping Jacobian
#' @import compiler
#' @export
dg <- cmpfun(function(x, xa, v, w, plist1, plist2, h, hgrad) {
  p <- length(x)
  SPC <- double(p)
  Sum2 <- sum(w)*h(x)
  
  SPQ <- double(length(Sum2))
  n1 <- length(plist1)
  for (i in 1:n1) {
    PCi <- as.matrix(plist1[[i]](xa))
    SPC <- SPC + v[i]*PCi
  }
  Sum1 <- sum(v)*x - SPC
  n2 <- length(plist2)
  for (j in 1:n2) {
    PQj <- as.matrix(plist2[[j]](h(xa)))
    SPQ <- SPQ + w[j]*PQj
  }
  Sum2 <- Sum2 - SPQ
  return(Sum1 + hgrad(x)%*%Sum2)
})

#' Compute the inverse approximate Hessian of the majorization using the Woodbury inversion formula.
#
#' \code{wood_inv_solve} computes the inverse of the Hessian term of the majorization of the proximity function
#' using the Woodbury formula. The function \code{mmqn_step} invokes \code{wood_inv_solve} instead of {ddg} if the argument
#' \code{woodbury=TRUE}. This should be used when p << n.
#'
#' @param x non-anchor point
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param hgrad Handle for output mapping Jacobian
#' @param df Right hand side
#' @import compiler
#' @export
wood_inv_solve <- cmpfun(function(x, v, w, hgrad, df) {
  n <- length(x)
  vv <- sum(v); ww <- sum(w)
  dh <- as.matrix(hgrad(x))
  p <- dim(dh)[2]  
  return( (1/vv)*df - (ww/vv^2)*dh %*% solve(diag(p) + (ww/vv)*t(dh) %*% dh , t(dh)%*%df ) )
})
