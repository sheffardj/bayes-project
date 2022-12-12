source(paste0(getwd(),'/posteriors/multi_modal_1d.R'))

#some devtool here.

## guess at a  prior
prior <- function(x) rnorm(1, x, 7)

if(1){ #### EXTRACT THIS STEP TO NEW FILE AND SOURCE()
  ## mcmc step?
  step <- function(x, posterior, prior) {
    ## Pick new point
    xp <- prior(x)
    ## Acceptance probability:
    alpha <- min(1, posterior(xp) / posterior(x))
    ## Accept new point with probability alpha:
    if (runif(1) < alpha)
      x <- xp
    ## Returning the point:
    x
  }
  
  ## mcmc step?
  run <- function(x, posterior, prior, nsteps) {
    res <- matrix(NA, nsteps, length(x))
    for (i in seq_len(nsteps))
      res[i,] <- x <- step(x, posterior, prior)
    drop(res)
  }
}

if(0){ # NOT RUN! change  0 to 1 to force run
  ## mcmc trajectory
  res <- run(-10, posterior, prior, 1000)
  layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
  par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
  plot(res, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
  usr <- par("usr")
  xx <- seq(usr[3], usr[4], length=301)
  plot(posterior(xx), xx, type="l", yaxs="i", axes=FALSE, xlab="")
  par(mfrow=c(1,1))
}

if(0){ # NOT RUN!
  ## plot of toy example
  hist(res, 50, freq=FALSE, main="", ylim=c(0, 1.5*max(hist(res)$density)), las=1,
       xlab="x", ylab="Probability density")
  z <- integrate(posterior, -Inf, Inf)$value
  curve(posterior(x) / z, add=TRUE, col="red", n=200)
}

# PROPER RUN
res.long <- run(-10, posterior, prior, 50000)
hist(res.long, 100, freq=FALSE, main="", ylim=c(0, 1.5*max(hist(res)$density)), las=1,
     xlab="x", ylab="Probability density", col="grey")
z <- integrate(posterior, -Inf, Inf)$value
curve(posterior(x) / z, add=TRUE, col="red", n=200)


if(0){ # NOT RUN
  ## examples of using different priors with more extreme beliefs
  res.fast <- run(-10, posterior, function(x) rnorm(1, x,  33), 1000)
  res.slow <- run(-10, posterior, function(x) rnorm(1, x,  .3), 1000)
  layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
  par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
  plot(res, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1,
       col="grey")
  lines(res.fast, col="red")
  lines(res.slow, col="blue")
  plot(posterior(xx), xx, type="l", yaxs="i", axes=FALSE)
  
  par(mfrow=c(1, 3), mar=c(4, 2, 3.5, .5))
  acf(res.slow, las=1, main="Small steps")
  acf(res, las=1, main="Intermediate")
  acf(res.fast, las=1, main="Large steps")
}
