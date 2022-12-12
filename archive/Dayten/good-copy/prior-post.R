dlaplace <-
  function(x, mean=0, sd=1) {
    # Laplace (double exponential) density function with mean equal to \code{mean} and standard deviation equal to \code{sd}. 
    # 'x': Vector of quantiles.
    # 'mean': Population mean.
    # 'sd': Population standard deviation.
    # example:    dlaplace( seq( 20, 80, length.out=11 ), 50, 10 )
    if (!is.numeric(x))  stop("'x' must be numeric.")
    if (!is.numeric(mean))  stop("'mean' must be numeric.")
    if (!is.numeric(sd))  stop("'sd' must be numeric.")
    if (sd<0)  stop("'sd' cannot be negative.")
    if (sd==0)  return( dnorm(x, mean, 0) )
    exp(-abs(x-mean)*sqrt(2)/sd)/(sd*sqrt(2))
  }

# priors for shape parameters
p1 <<- function (x, var) {
  dunif(x, var-3, var+3) #neg values ok
}

p2 <<- function (x, var) {
  dnorm(x, var, 3)
}  
  
p3 <<- function(x, var) {
  dlaplace(x, var, 3)
}

# our likelihood function
lh <- function(mu, sigma, pi, pj) {
  LL_obs <- mapply(dlogitnorm, mu=mu, sigma=sigma, MoreArgs = list(x=og_data), SIMPLIFY = T)
  LL_mat <- matrix(apply(LL_obs, 2, prod), nrow=dim(mu)[1], ncol=dim(sigma)[1])
  LL_mat * pi(mu, var=mu0) * pj(sigma, var=sigma0)
}

# intg <- integrate(function(x) p2(x, 3.16, T), lower=-Inf, upper=Inf)
# curve((1/(intg$value))*p2(x, 3.16, T), from=-10, to=20)
# 
# integrate(function(x) (1/(intg$value))*p2(x, 3.16, T), lower=-Inf, upper=Inf)
