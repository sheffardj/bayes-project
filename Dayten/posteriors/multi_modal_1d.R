p <- 0.4
mu <- c(-5, 1)
sd <- c(5, 2)

posterior <- function(x){
  p * dnorm(x, mu[1], sd[1]) + (1-p) * dnorm(x, mu[2], sd[2])
}

curve(posterior(x), col='red', -15, 15, n = 301, las=1)
