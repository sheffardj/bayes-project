##### PREAMBLE #####
set.seed(8776)
library(logitnorm)
library(sn)
library(dplyr)
library(philentropy)
library(bayestestR)
library(ggplot2)
library(ggpubr)

##### PRIORS #####

un <- function (d, v) {
  dunif(d, v-1,v+1 )
}
nor <- function (d, v){
  dnorm(d, 0, 1)
}
lp <- function(x, mean, sd = 1) {
  # Laplace (double exponential) density function with mean equal to \code{mean} and standard deviation equal to \code{sd}. 
  # 'x': Vector of quantiles.
  # 'mean': Population mean.
  # 'sd': Population standard deviation.
  # example:    dlaplace( seq( 20, 80, length.out=11 ), 50, 10 )
  if (!is.numeric(x))  stop("'x' must be numeric.")
  if (!is.numeric(mean))  stop("'mean' must be numeric.")
  if (!is.numeric(sd))  stop("'sd' must be numeric.")
  if (sd<0)  stop("'sd' cannot be negative.")
  if (sd==0)  return( dnorm(x, mean, 0) )
  exp(-abs(x-mean)*sqrt(2)/sd)/(sd*sqrt(2))
}

##### PARAMETERS #####
n <- 30 # number of samples
data <- matrix(nrow = n, ncol = 9) # storage for original data
posteriors <- matrix(nrow = 2 * n + 1, ncol = 54) # storage for found posteriors

base_cases <- expand.grid(mu=0:2, sigma=c(3.16, 1.78, 0.32))
base_cases$case_id <- 1:9

cases <- expand.grid( # lay out all the cases
  mu = 0:2,
  sigma = c(3.16, 1.78, 0.32),
  prior = c('un', 'nor', 'lp'),
  fixed = c("mu", "sigma")
)
cases$best_est <- NA

# Enforce identical `original data` for all cases with similar `mu` and `sigma`
for(i in 1:9){
  data[,i] <- rlogitnorm(n, base_cases[i,'mu'], base_cases[i,'sigma'])
}

##### CASE FUNCTION #####
## Includes Simpson's integration and Posterior estimations
## Write current Case data to Global Environment (`dat` and `posterior`)
## Also updates `cases$best_est` in Global Env.
call_case <- function(fixed, fixed_val, free_val, prior, grid_size=200){

  
  if(fixed == "mu"){
    case_id <- base_cases %>% 
      filter(mu == fixed_val & sigma == free_val) %>% 
      pull(case_id)
    dat <- data[, case_id]
  } else {
    case_id <- base_cases %>% 
      filter(mu == free_val & sigma == fixed_val) %>% 
      pull(case_id)
    dat <- data[, case_id]
  }

  # Toggle variable estimate whether `fixed` is mu or sigma 
  # If mu fixed, get sample VARIANCE. If not (so sigma fixed), find sample MEAN.
  var0 <- ifelse(fixed=="mu", var(dat), mean(dat))
  
  # set up simpsons params
  a <- ifelse(fixed=="mu",max(var0 - 1, 0.01), var0 - 1) # toggled on `fixed`
  b <- var0 + 1
  nn = 2*n+1 
  h = (b - a) / (nn - 1) #length of subdivision
  xx <- seq(a, b, length=nn)
  
  # obtain likelihood
  if(fixed == "mu"){
    # Likelihood with mu fixed
    LL_obs <- mapply(dlogitnorm, mu=fixed_val, sigma=xx, MoreArgs = list(x=dat), SIMPLIFY = T)
  } else {
    # Likelihood with sigma fixed
    LL_obs <- mapply(dlogitnorm, mu=xx, sigma=fixed_val, MoreArgs = list(x=dat), SIMPLIFY = T)
  }
  LL <- abs(apply(LL_obs, 2, prod))
  
  # obtain
  y <- LL*prior(xx, var0)
  # y %>% plot(type='l')
  
  # compute the denominator with Simpsons rule
  dnom <-
    h * (y[1] + y[nn] + 4 * sum(y[(2:(nn - 1))[2:(nn - 1) %% 2 == 0]]) + 2 *
           sum(y[(2:(nn - 1))[2:(nn - 1) %% 2 == 1]])) / 3
  
  # compute the posterior and assign to Global Env,.
  posterior <- c(y/dnom)
  
  # finds the best estimate and stores to global env.
  best_est_index <- which.max(posterior)
  best_est <- xx[best_est_index]
  
  return(list(dat = dat, posterior = posterior, best_est=best_est))
}

##### MAIN LOOP #####
for(ii in 1:dim(cases)[1]){
  case <- cases[ii, ] # pull out case
  
  # verbosely assign case variables (fixed and free)
  fixed <- as.character(case$fixed)
  fixed_val <- ifelse(case$fixed == 'mu', case$mu, case$sigma)
  free_val <- ifelse(case$fixed == 'mu', case$sigma, case$mu)
  
  # pull the prior from the environment from case assignment (i.e., get("un"))
  prior <- get(as.character(case$prior))
  
  results <- call_case(fixed, fixed_val, free_val, prior)
  cases[ii,"best_est"] <- results$best_est
  posteriors[,ii] <- results$posterior
  if(ii %% 6 == 0) cat(paste0(round(100*ii/dim(cases))[1], '% '))
}




##### PLOTS #####
### MU FIXED
par(mfrow=c(3,3))
for(ii in c(0,1,2)){
  for(jj in c(3.16, 1.78, 0.32)){
    group_priors <- as.data.frame(cases) %>% 
      filter(fixed=='mu' & mu==ii & sigma == jj)
    
    case_id <- base_cases %>% 
      filter(mu == ii & sigma == jj) %>% 
      pull(case_id)
    dat <- data[, case_id]
    
    best_un <- group_priors[1,'best_est']
    best_nor <- group_priors[2, 'best_est']
    best_lp <- group_priors[3, 'best_est']
    
    sig.hat <- round(mean(best_un, best_nor, best_lp), 3) # just a ref for plotting
    
    # this gets a nice `ymax` since any curve or histogram prob can be the max
    hist <- hist(dat, plot=F)
    dens_og <- dlogitnorm(seq(0,1,length=401), ii, jj)
    dens_un <- dlogitnorm(seq(0,1,length=401), ii, best_un)
    dens_nor <- dlogitnorm(seq(0,1,length=401), ii, best_nor)
    dens_lp <- dlogitnorm(seq(0,1,length=401), ii, best_lp)
    ymax <- max(c(hist$density, dens_og, dens_un, dens_nor, dens_lp))
    
    hist(dat, probability = T, ylim=c(0,ymax+1), xlim=c(0,1),
         main=bquote(mu~" = "~.(ii)~" "~sigma~" = "~.(jj)~", "~bar(sigma)~" = "~.(sig.hat)))
    curve(dlogitnorm(x, ii, jj), col = 'red', add = T)
    curve(dlogitnorm(x, ii, best_un), col ='blue', lty='dashed', add = T)
    curve(dlogitnorm(x, ii, best_nor), col ='green', lty='dashed', add = T)
    curve(dlogitnorm(x, ii, best_lp), col ='orange', lty='dashed', add = T)
  }
}


### SIGMA FIXED
par(mfrow=c(3,3))
for(ii in c(0,1,2)){
  for(jj in c(3.16, 1.78, 0.32)){
    group_priors <- as.data.frame(cases) %>% 
      filter(fixed=='sigma' & mu==ii & sigma == jj)
    
    case_id <- base_cases %>% 
      filter(mu == ii & sigma == jj) %>% 
      pull(case_id)
    dat <- data[, case_id]
    
    best_un <- group_priors[1,'best_est']
    best_nor <- group_priors[2, 'best_est']
    best_lp <- group_priors[3, 'best_est']
    
    mu.hat <- mean(best_un,best_nor,best_lp) # just a ref for plotting
    
    # this gets a nice `ymax` since any curve or histogram prob can be the max
    hist <- hist(dat, plot=F)
    dens_og <- dlogitnorm(seq(0,1,length=401), ii, jj)
    dens_un <- dlogitnorm(seq(0,1,length=401), best_un, jj)
    dens_nor <- dlogitnorm(seq(0,1,length=401), best_nor, jj)
    dens_lp <- dlogitnorm(seq(0,1,length=401), best_lp, jj)
    ymax <- max(c(hist$density, dens_og, dens_un, dens_nor, dens_lp))
    
    hist(dat, probability = T, ylim=c(0,ymax+1),
         main=bquote(atop(mu~" = "~.(ii)~" "~sigma~" = "~.(jj)~", "~
                            bar(mu)~" = "~.(mu.hat))))
    curve(dlogitnorm(x, ii, jj), col = 'red', add = T)
    curve(dlogitnorm(x, best_un, jj), col ='blue', lty='dashed', add = T)
    curve(dlogitnorm(x, best_nor, jj), col ='green', lty='dashed', add = T)
    curve(dlogitnorm(x, best_lp, jj), col ='orange', lty='dashed', add = T)
  }
}
