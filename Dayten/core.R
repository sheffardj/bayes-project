xx <- rlogitnorm(10000,0.1,2.5)
pp <- exp(xx) / (1+sum(exp(xx)))


f <- function(x, mu, sigma) {
  scale <- (1 / (sigma * sqrt(2 * pi)))
  jacobian <- (1 / (x * (1 - x)))
  expo <- exp(-((logit(x) - mu) ^ 2) / (2 * sigma ^ 2))
  return(scale * expo * jacobian)
}

logit(runif(10000, 0, 1)) %>% hist # takes (0,1) to normal
f(runif(10000,0,1), 0, 1) %>% hist # takes (0,1) to logit-normal

f(pp, 0.1, 2.5) %>% hist # takes (0,1) to logit-normal
pp %>% hist


yy <- rnorm(10000, 0.1, 2.1) #normal
yy %>% hist

pp <- exp(xx) / (1+sum(exp(xx))) #logit transform
f(pp, 0.1, 2.1) %>% hist # takes (0,1) to logit-normal
pp %>% hist






