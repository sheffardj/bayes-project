## "given" some data
xx <- rlogitnorm(10000,0.1,2.5)

## we decide to fit a logit-normal dist'n

## here is a custom logit-normal pdf
## but it doesn't draw the legs !
f <- function(x, mu, sigma) {
  scale <- (1 / (sigma * sqrt(2 * pi)))
  expo <- exp(-((logit(x) - mu) ^ 2) / (2 * sigma ^ 2))
  jacobian <- (1 / (x * (1 - x)))
  return(scale * expo * jacobian)
}


est_mu <- mean(rnorm(xx, 0.1, 2.5) * dlogitnorm(xx))

hist(xx, breaks=50, probability = T, ylim=c(0,10))
curve(dlogitnorm(x, 0.1, 2.5), from=0, to=1, add=T, col='blue')
curve(f(x, 0.1, 2.5), from=0, to=1, add=T, col='blue')
curve(f(x, est_mu, 2.5), from=0, to=1, add=T, col='red')

