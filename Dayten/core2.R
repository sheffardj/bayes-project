# Set the numer of tosses.  
n <- 100
# Set the number of heads obtained.  
h <- 73  
# Define our likelihood function.   
# Since our model is a binomial model, we can use:  
likelihood <- function(h,n,p){  
  lh <- dbinom(h,n,p)  
  lh  
}  
# Set the starting value of p  
p <- runif(1,0,1)  
# Create an empty data.frame to store the accepted p values for each iteration.  
# Remember: "the posterior probability is just an updated version of the prior"  
posterior <- data.frame()  
# Set the lenght of the loop (Marcov Chain, number of iterations).  
nrep <- 5000  
# Start the loop (MCMC)  
for (i in 1:nrep) {  
  # Obtain a new proposal value for p  
  p_prime <- p + runif(1, -0.05,0.05)  
  # Avoid values out of the range 0 - 1  
  if (p_prime < 0) {p_prime <- abs(p_prime)}  
  if (p_prime > 1) {p_prime <- 2 - p_prime}  
  # Compute the acceptance proability using our likelihood function and the  
  # beta(1,1) distribution as our prior probability.  
  R <- likelihood(h,n,p_prime)/likelihood(h,n,p) * (dbeta(p_prime,1,1)/dbeta(p,1,1))  
  # Accept or reject the new value of p  
  if (R > 1) {R <- 1}  
  random <- runif (1,0,1)  
  if (random < R) {  
    p <- p_prime  
  }  
  # Store the likelihood of the accepted p and its value  
  posterior[i,1] <- log(likelihood(h, n, p))  
  posterior[i,2] <- p  
  print(i)  
}  