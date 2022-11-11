library(dplyr)

# GROUND TRUTH
N <- 5
theta <- 0.8
# y <- rbinom(N, 1, theta) COULD RUN THIS
# but book uses:
y <- c(1,0,1,1,1)


# SECTION 4.6.2.1
## BERNOULLI LIKELIHOOD == ## Binomial Likelihood
lik_bern <- function(y, theta) {
  vals <- c()
  for(parm in theta){
    bernoulli <- choose(N, sum(y))*(parm ^ sum(y)) * ((1 - parm) ^ (N - sum(y)))
    vals <- c(vals, bernoulli)
  }
  vals
}

theta_seq <- seq(0, 1, by=0.01)
lik_bern(y, theta_seq) %>% plot(theta_seq, ., type='l')
#notice the spike around 0.8 - the true theta value


## SECTION 4.6.2.5
# NEED TO USE BETA DISTRIBUTIONS FOR REST:

## SET UP PRIOR (Observed 2 heads and 2 tails before seeing data)
a <- 2 #heads
b <- 2 #tails

x_seq <- seq(0, 1, by=0.01)
y_sum <- y %>% sum()

update <- function (a, b){
  prior <- dbeta(x_seq, a, b)
  lik <- dbeta(x_seq, y_sum + 1, N - y_sum + 1)
  posterior <-  dbeta(x_seq, y_sum + a, N - y_sum + b)
  
  plot(x_seq, prior, type = 'l', col = 'red', ylab = 'y', ylim = c(0, 3))
  lines(x_seq, lik, lty = 4, col='green')
  lines(x_seq, posterior, lty = 2, col = 'blue')
}

update(a, b)


## UNIFORM PRIOR
a <- 1 #heads
b <- 1 #tails
update(a, b)