#### Practical 2 Group 64: Xinyu HOU(s2145863), Maeve LI(Minqing LI)(s2167017), Di WU(s2176435)
#### Github repo: https://github.com/MaeveLi/Practical-4-Group-64

## Overview: the aim is to write a function named "bfgs" that implements the BFGS 
## quasi-Newton minimization method for general purpose optimization.

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  ## a function implementing the BFGS quasi-Newton minimization method.
  
  ## inputs: 
  ## theta: a set of initial values for parameters to be optimized in f.
  ## f: the objective function to be minimized.
  ## tol: the convergence tolerance.
  ## fscale: a rough estimate of the magnitude of f at the optimum which is used in convergence testing.
  ## maxit: the maximum number of BFGS iterations.
  
  ## output:
  ## a named list that contains:
  ## f: the scalar value of the objective function at the minimum.
  ## theta: the vector of values of the parameters at the minimum.
  ## iter: the number of iterations taken to reach the minimum.
  ## g: the gradient vector at the minimum (so the user can judge closeness to numerical zero).
  ## H: the approximate Hessian matrix (obtained by finite differencing) at the minimum.
  
  n = length(theta) ## length of theta
  B = diag(n) ## initialize unit matrix
  x = theta
  i = 0
  
  g <- function(theta,f,...){
    ## A function that computes the gradient when given a set of parameter values (theta) and a function f.
    if (is.null(attr(f(theta,...),"gradient"))){ 
      ## if the function does not have an attribute "gradient", compute the gradient
      ## using the finite differencing method
      th0 <- fd <- theta
      f0 <- f(th0,...)
      eps <- 1e-7 ## finite difference interval
      for (i in 1:n) { ## loop over parameters
        th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps
        f1 <- f(th1,...) ## compute resulting f 
        fd[i] <- (f1 - f0)/eps ## approximate -dl/dth[i]
      }
    }
    else{
      ## if the gradient is provided in the function we directly use it
      fd=attr(f(theta,...),"gradient")
    }
    return(fd)
  }
  
  if (is.infinite(f(theta)) || is.infinite(g(theta,f,...))) ## check if objective is infinite
    stop("the objective or derivatives are not finite at initial theta")
  
  for (k in 1:maxit){
    lambda = 1 ## initialize the first step length
    fx = f(theta=x) ## compute the current objective function value
    if (max(abs(g(x,f,...))) < (abs(fx)+fscale)*tol) ## to judge whether the gradient vector is close enough to zero
      break
    dx = - B %*% g(x,f,...) ## compute the step direction
    
    while(TRUE){
      ## a suitable lambda is found using this loop
      x1 = x + lambda * dx
      fx1 = f(theta=x1)
      if (is.infinite(fx1))
        lambda = 0.5*lambda ## reduce the step length if a step leads to a non-finite value
      else if (is.na(fx1))
        lambda = 0.5*lambda ## reduce the step length if a step leads to NA
      else if (fx1 > fx){ 
        lambda = 0.5*lambda ## reduce the step length if the objective function value isn't reduced
      }
      else if (t(g(x+lambda*dx,f,...)) %*% dx < 0.9*t(g(x,f,...)) %*% dx){
        lambda = 1.05*lambda ## increase the step length if the step length doesn't satisfy the second Wolfe condition
      }
      else
        break
    }
    
    x1 = x + lambda * dx  ## compute the new x based on the search direction and step length
    fx1 = f(theta=x1) ## compute the new f value
    s = x1 - x
    y = g(x1,f,...) - g(x,f,...)
    
    ## update the B matrix
    p = 1 / (t(s) %*% y)
    p = drop(p)
    ## rewrite the B formula to ensure that its cost is O(n^2)
    B = B - p*(B%*%y)%*%t(s) - p*s%*%(t(y)%*%B) + p^2*s%*%(t(y)%*%B%*%y)%*%t(s) + p*(s%*%t(s))
    
    i = i+1
    x = x1
  }
  
  if (i==maxit & max(abs(g(x,f,...))) >= (abs(fx)+fscale)*tol)
    stop("maxit is reached without convergence")
  
  ## Compute the Hessian matrix by finite differencing the final gradient
  H <- matrix(0,n,n) ## initialize matrix
  for (j in 1:n) { ## loop over parameters
    eps=1e-7
    th1 <- x; th1[j] <- th1[j] + eps ## increase x[j] by eps
    h1 <- g(th1,f,...) ## compute the resulting gradient of th1
    H[j,] <- (h1 - g(x,f))/eps ## approximate second derivatives
  }
  H <- 0.5 * (t(H) + H) ## ensure H's symmetricity
  
  ## return a list
  l = list(f=fx,theta=x,iter=i,g=g(x,f),H=H)
  return(l)
}
