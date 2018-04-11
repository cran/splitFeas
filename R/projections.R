#' Projection onto a ball
#' 
#' \code{project_ball} computes the Euclidean projection of a point onto a ball.
#' 
#' @param x Point to project
#' @param center Center of the sphere
#' @param r Radius of the sphere
#' @import compiler
#' @export
#' @examples
#' set.seed(12345)
#' p <- 3
#' center <- rnorm(p)
#' r <- runif(1)
#' x <- rnorm(p)
#' y <- project_ball(x,center,r)
project_ball <- cmpfun(function(x, center, r) {
  d <- x - center
  dl <- norm(as.matrix(d),'f')
  if (dl <= r) {
    return(x)
  } else {
    d <- d / dl
    return(center + r*d)
  }
})

#' Project onto a cube
#' 
#' \code{project_cube} computes the Euclidean projection of a point onto a cube.
#' 
#' @param x Point to project
#' @param center Center of the square
#' @param r Half the length of a side
#' @import compiler
#' @export
#' @examples
#' set.seed(12345)
#' p <- 3
#' center <- matrix(rnorm(p),p,1)
#' r <- runif(1)
#' x <- matrix(rnorm(p),p,1)
#' y <- project_cube(x,center,r)
project_cube <- cmpfun(function(x, center, r) {
  y1 <- min(max(center[1]-r, x[1]), center[1]+r)
  y23 <- project_square(x[-1], center[-1], r)
  y <- matrix(c(y1, y23), ncol=1)
  return(y)
})

#' Projection onto a halfspace
#' 
#' \code{project_halfspace} computes the Euclidean projection of a point onto a closed half-space.
#' The function returns the projection onto the set
#       C = {y : t(a)%*%y <= b}
#' 
#' @param x Point to project
#' @param a is the normal vector
#' @param b is the threshold
#' @import compiler
#' @export
#' @examples
#' set.seed(12345)
#' p <- 3
#' a <- matrix(rnorm(p),p,1)
#' a <- a/norm(a,'f')
#' b <- runif(1)
#' x <- matrix(rnorm(p),p,1)
#' y <- project_halfspace(x,a,b)
project_halfspace <- cmpfun(function(x, a, b) {
  x <- as.matrix(x)
  y <- x
  a <- as.matrix(a)
  if (t(a)%*%x > b) {
    y <- x - a%*%((t(a)%*%x - b)/norm(a,'f')**2)
  }
  return(y)
})


#' Project onto a square
#' 
#' \code{project_square} computes the Euclidean projection of a point onto a square.
#' 
#' @param x Point to project
#' @param center Center of the square
#' @param r Half the length of a side
#' @import compiler
#' @export
#' @examples
#' set.seed(12345)
#' p <- 2
#' center <- matrix(rnorm(p),p,1)
#' r <- runif(1)
#' x <- matrix(rnorm(p),p,1)
#' y <- project_square(x,center,r)
project_square <- cmpfun(function(x, center, r) {
  d <- length(center)
  if (norm(as.matrix(x-center),'I') <= r) y <- x
  if (x[1] - center[1] > r & x[2] - center[2] > r) y <- matrix(c(center[1] + r, center[2] + r), d, 1)
  if (x[1] - center[1] < -r & x[2] - center[2] > r) y <- matrix(c(center[1] - r, center[2] + r), d, 1)
  if (x[1] - center[1] < -r & x[2] - center[2] < -r) y <- matrix(c(center[1] -r ,center[2] - r), d, 1)
  if (x[1] - center[1] > r & x[2] -center[2] < -r) y <- matrix(c(center[1] + r, center[2] - r), d, 1)
  if (abs(x[1] - center[1]) <= r & x[2] - center[2] > r) y <- matrix(c(x[1], center[2]+r), d, 1)
  if (abs(x[1] - center[1]) <= r & x[2] - center[2] < -r) y <- matrix(c(x[1], center[2]-r), d, 1)
  if (x[1] - center[1] > r & abs(x[2] - center[2]) <= r) y <- matrix(c(center[1]+r, x[2]), d, 1)
  if (x[1] - center[1] < -r & abs(x[2] - center[2]) <= r) y <- matrix(c(center[1]-r, x[2]), d, 1)
  return(y)
})
