#' One step of self-adaptive projection-type method for the NMSFP
#' 
#' \code{nmsfp_sap_one_step} performs the self-adaptive projection-type method of Li et al.
#' 
#' @param x current iterate
#' @param delta step-size parameter
#' @param mu step-size parameter
#' @param tau step-size parameter
#' @param gamma step-size parameter
#' @param df handle to gradient of objective function
#' @param proj handle to projection operation.
#' @import compiler
#' @export
#' @references Li, Han, and Zhang. (2013) A self-adaptive projection-type method for nonlinear 
#' multiple-sets split feasibility problem, Inverse Problems in Science and Engineering,
#' 21(1):155-170.
nmsfp_sap_one_step <- cmpfun(function(x, delta, mu, tau, gamma, df, proj) {
  beta <- gamma
  gradf <- df(x)
  repeat {
    xp <- proj(x - beta*gradf)
    gradfp <- df(xp)
    s <- as.matrix(gradf - gradfp)
    dx <- as.matrix(x - xp)
    lhs <- beta*t(s)%*%s
    rhs <- t(dx)%*%s
    if (lhs <= (2-delta)*rhs) break
    beta <- mu*beta
  }
  if (lhs <= 0.5*rhs) {
    gamma <- (1 + tau)*beta
  } else {
    gamma <- beta
  }
  return(list(x=xp,gamma=gamma))
})
