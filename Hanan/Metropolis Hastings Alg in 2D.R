##### An example of MCMC with a MVN posterior
source("seed.R")
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}

#### METROPOLIS-HASTINGS ALGORITHM 

## proposal function needs to be symmetric, as per M-H Alg (`Detailed Balance' property).
proposal <- function(x, d=8) x + runif(length(x), -d/2, d/2)

## What's handy about M-H is we don't need a transition matrix P
## State-transitions function
step <- function(x, posterior, jump) {
  ## Pick new point
  xp <- jump(x)
  ## Acceptance probability:
  alpha <- min(1, ifelse(posterior(x)==0, 1e-5, posterior(xp)/posterior(x)) )
  #alpha <- min(1, (posterior(xp) / posterior(x)))
  ## Accept new point with probability alpha:
  if (runif(1) < alpha)
    x <- xp
  ## Returning the point:
  x
}

## MCMC wrapper (input is initial state, posterior function, the proposal dist'n, and # of steps)
run <- function(x, posterior, jump, nsteps) {
  res <- matrix(NA, nsteps, length(x))
  for (i in seq_len(nsteps))
    res[i,] <- x <- step(x, posterior, jump)
  drop(res)
}


## posterior (2D) weighted sum of two distn's
#p=0.75;
mu=c(5,1); sd=c(2.5,4)
#a0=3; b0=2
posterior <- function(x){
  #p*dnorm(x, mu, sd)+(1-p)*dgamma(x, a0, b0) 
  dmvnorm(x, mean = mu, sigma = matrix(c(rev(sd), (sd)), byrow=T, nrow=2))
}
#par(mfrow=c(1,1))
#curve(posterior(x), col='red', -5, 15, n = 301, las=1)


#call the MCMC alg functions
res.long <- run(c(5,1), posterior, proposal, 1000)


### Run diagnostics for convergence/etc.



#plots!
{x <- seq(-1, 11, length=71)
y <- seq(-6, 9, length=61)
xy <- expand.grid(x=x, y=y)
z <- matrix(apply(as.matrix(xy), 1, posterior), length(x), length(y))

image(x, y, z, xlim=range(x, res.long[,1]), ylim=range(x, res.long[,2]))
contour(x, y, z, add=TRUE)
lines(res.long[,1], res.long[,2], col="#00000088", lwd=0.75)
}
# hist(res.long[,1], 100, freq=FALSE, main="", ylim=c(0, 2.5*max(hist(res.long[,1])$density)), las=1,
#        xlab="x", ylab="Probability density", col="grey")


###################### END

## Things to consider:
# multi-dimensional posterior distributions
# burn-in iterations?
# lagged states (due to autocorrelations) & effective sample size
# tuning the proposal density variance (?)
# other samplers, like Gibb's for example?
# MCMC Diagnostics!!!
# model comparisons? (if we're doing models??)
# GOAL: want to compare the effect of different priors for a given posterior/model using MH-MCMC (?)


