library(dplyr)

## edited updating function
update <- function (a, b, i = 0, trial){
  posterior <- dbeta(x_seq, a, b)
  y_max <- max(posterior)
  plot(x_seq, posterior, type = 'l', col = (trial + 2), ylab = 'y', ylim = c(0, y_max))
  mtext(paste("trial", i), side=3, cex=0.8)
  mtext(ifelse(trial == 1, "Success", "Fail"), side=4, cex=0.6)
}

## Set up parms
N <- 89 # 10 bernoulli trials
theta <- 0.6 # biased to success
# will iterate for "y"
# and use uninformative prior
trials <- c()

a <- 1
b <- 1
par(mfrow=c(5,6), mar=c(2,2,2,2))
update(a,b, 0, 1)
## in this loop, we perform a trial and 
## first we set a & b based on uninformative prior
## then we update a & b based on new knowledge (and the prior)
for(i in 1:65){
  trial <- rbinom(1, 1, theta) #success or fail
  trials <- c(trials, trial) # append current trial
  num_success <- trials %>% sum
  num_fails <- i - num_success
  
  a <- ifelse(length(trials) == 0, 1, num_success + 1)
  b <- ifelse(length(trials) == 0, 1, num_fails + 1)
  
  if(i %in% 1:11 | i %% 3 == 0){
    update(a, b, i, trial)
  }
}
