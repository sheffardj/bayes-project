##### An example of MCMC with a JANKY posterior
source("seed.R")

#### METROPOLIS-HASTINGS ALGORITHM 


## proposal function needs to be symmetric, as per M-H Alg (`Detailed Balance' property).
# best if we just use the normal distribution, 
# but we can play around with the variance hyperparameter. 
proposal <- function(x) rnorm(1, x, 1)

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


## posterior (1D) weighted sum of two distn's
p=0.75;
mu=5; sd=2.5
a0=3; b0=2
posterior <- function(x){
  p*dnorm(x, mu, sd)+(1-p)*dgamma(x, a0, b0) 
}
#par(mfrow=c(1,1))
curve(posterior(x), col='red', -5, 15, n = 301, las=1)


#call the MCMC alg functions
res.long <- run(1, posterior, proposal, 10000)


### Run diagnostics for convergence/etc.



#plots!
hist(res.long, 100, freq=FALSE, main="", ylim=c(0, 1.5*max(hist(res.long)$density)), las=1,
       xlab="x", ylab="Probability density", col="grey")
#z <- integrate(posterior, -Inf, Inf)$value
#curve(posterior(x) / z, add=TRUE, col="red", n=200)
curve(posterior(x), add=TRUE, col="blue", n=200) 
  

{layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
plot(res.long, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1, lwd=0.25)
usr <- par("usr")
xx <- seq(usr[3], usr[4], length=301)
plot(posterior(xx), xx, type="l", yaxs="i", axes=FALSE, xlab="")
par(mfrow=c(1,1))
}

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


