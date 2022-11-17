library(logitnorm)

# one-param estimation
dat <- rlogitnorm(10000, 0.1, 2.1) # mu= 0.1 sigma = 2.1
hist(dat, breaks=50, probability = T)
sigma <- 2.1
hist(rbeta(10000, 0.5, 0.5), breaks=50, probability = T)
runif(10000, -1, 1)
hist(dbeta(dat, 0.5, 0.5)*dnorm(dat, 0, 1), breaks=1000, probability = T)


hist(dat, breaks=50, probability = T)
curve(dlogitnorm(x, mu=0.1, sigma=2.1), from=0, to=1, add=T, col='blue')

mu0 <- min(dbeta(dat, 0.5, 0.5)*dnorm(dat, 0, 1))
# hist(dbeta(dat, 0.5, 0.5)*dnorm(dat, 0, 1), breaks=100)
curve(dlogitnorm(x, mu=mu0, sigma=2.1), from=0, to=1, add=T, col='red')

mu1 <- min(dbeta(dat, 0.5, 0.5)*dunif(dat, -1, 1))
# hist(dbeta(dat, 0.5, 0.5)*dunif(dat, -1, 1), breaks=100, prob=T)
curve(dlogitnorm(x, mu=mu1, sigma=2.1), from=0, to=1, add=T, col='green')

mu2 <- min(dbeta(dat, .5,.5))
# hist(dbeta(dat, .5,.5), breaks = 50, probability = T)
curve(dlogitnorm(x, mu=mu2, sigma=2.1), from=0, to=1, add=T, col='orange')

