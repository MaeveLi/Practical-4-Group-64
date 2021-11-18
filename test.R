## test

rb <- function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by 'bfgs'
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
}

rbf <- function(theta,k=10) {
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  f
}
rbg <- function(theta,k=10) {
  z <- theta[1]; x <- theta[2]
  c(2*k*(z-x^2),
    -4*k*x*(z-x^2) -2*(1-x))
}
optim(c(-1,2),rbf,rbg,method="BFGS",hessian=TRUE)

bfgs(c(-1,2),rb)
