source(paste0(getwd(),'/posteriors/multi_modal_2d.R'))
## inherit from core_1d.R
if(1){
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
  
  run <- function(x, posterior, prior, nsteps) {
    res <- matrix(NA, nsteps, length(x))
    for (i in seq_len(nsteps))
      res[i,] <- x <- step(x, posterior, prior)
    drop(res)
  }
}

## example MCMC sampling of 2d multi-modal
q <- function(x, d=8) {x + runif(length(x), -d/2, d/2)}

x0 <- c(-4, -4)
samples <- run(x0, f, q, 1000)

image(x, y, z, xlim=range(x, samples[,1]), ylim=range(x, samples[,2]))
contour(x, y, z, add=TRUE)
lines(samples[,1], samples[,2], col="#00000088")

## drawing a f-ton with MCMC
samples <- run(x0, f, q, 100000)

## take a loook
smoothScatter(samples)
contour(x, y, z, add=TRUE)

## marginalized dist in 'x'
hist(samples[,1], freq=FALSE, main="", xlab="x",
     ylab="Probability density")

## marginalized dist in 'y'
hist(samples[,2], freq=FALSE, main="", xlab="y",
     ylab="Probability density", breaks=80)

## get density in 'x'
m <- function(x1) {
  g <- Vectorize(function(x2) f(c(x1, x2)))
  integrate(g, -Inf, Inf)$value
}

xx <- seq(min(samples[,1]), max(samples[,1]), length=201)
yy <- sapply(xx, m)
z <- integrate(splinefun(xx, yy), min(xx), max(xx))$value

hist(samples[,1], freq=FALSE, main="", las=1, xlab="x",
     ylab="Probability density of 'x'", ylim=c(0, 0.25))
lines(xx, yy/z, col="red")


## get density in 'y'
m2 <- function(x2) {
  g <- Vectorize(function(x1) f(c(x1, x2)))
  integrate(g, -Inf, Inf)$value
}

xx <- seq(min(samples[,2]), max(samples[,2]), length=201)
yy <- sapply(xx, m2)
z <- integrate(splinefun(xx, yy), min(xx), max(xx))$value

hist(samples[,2], freq=FALSE, main="", las=1, xlab="y",
     ylab="Probability density of 'y'", ylim=c(0, 0.25))
lines(xx, yy/z, col="red")
