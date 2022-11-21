# data that we're "given"
zz <- rlogitnorm(1e4, 0.1, 2.1)
zz %>% hist(., probability=T)

# we assume it is logit-normal, so we convert data to normal and get sample params
mu0 <- logit(zz) %>% mean()
sigma0 <- sqrt(logit(zz) %>% var())

# but we will need vector estimates for the prior for mu, assuming ~ Normal
sample_mu <- rnorm(1e4, mu0, 0.1)
sample_mu %>% hist(., probability=T)

# and again for sigma, assuming ~ Uniform
sample_sigma <- runif(1e4, sigma0 - 1, sigma0 + 1)
sample_sigma %>% hist(., probability=T)

# Take a look the prior distribution with the estimated mu and sigma vectors
prior_zz <- rlogitnorm(1e4, sample_mu, sample_sigma)
prior_zz %>% hist(., probability=T)

# data plotted against reality (blue) we want to match
curve(dlogitnorm(x, 0.1, 2.1), from=0, to=1, add=T, col='blue')

# grid approximation
mu.list <- seq(from = mu0 - 0.1,
               to = mu0 + 0.1,
               length.out = 100)

sigma.list <- seq(from = sigma0 - 0.1,
                  to = sigma0 + 0.1,
                  length.out = 100)

# prepare a list of possible mu and sigma
post <- expand.grid(mu = mu.list, sigma = sigma.list)

##### This part is a big step!! Read carefully, follow closely.
# Remember that in large samples, all unique samples are VERY unlikely.
# So, we have to be careful to do everything on the log scale or else
# the rounding error from R will make posterior probabilities ZERO!!!
#
# Thus, for each pair of ('mu', 'sigma') from the 'post' matrix
# we want to add together all the log-likelihoods of the data conditioned on 
# each mu[i] and sigma[i] pair!
post$LL <- sapply(1:nrow(post), function(i) {
  sum(dlogitnorm(zz, post$mu[i], post$sigma[i], log = TRUE))
})

##### Another important step!!
# Normally, we would want to multiply the prior by the likelihood right but 
# we're in the in the log scale, so we add the priors (also in the log scale!) 
# to the likelihood. (this is equivalent to multiplying but less rounding error)
post$prod <- post$LL + dnorm(post$mu, mu0, 0.1, log = TRUE) +
  dunif(post$sigma, sigma0 - 1, sigma0 + 1, log = TRUE)

# Now we want to convert back up to the probability scale, but the naive attempt
# of 'exp(post$prod)' would send most probabilities to ZERO!
# So we can instead scale all the log-products by the maximum log-product
# These are called relative posterior probabilities and NOT exact probabilities.
# Perhaps this is where integrating can help us achieve a better result?
post$prob <- exp(post$prod - max(post$prod))

# These functions come from the 'rethinking' package
contour_xyz(post$mu, post$sigma, post$prob)
image_xyz(post$mu, post$sigma, post$prob)

# an alternative visualization
sample.rows <-sample(1:nrow(post),size=1e4,replace=TRUE,
                          prob=post$prob )
sample.mu <-post$mu[sample.rows]
sample.sigma <-post$sigma[sample.rows]
plot(
  sample.mu,
  sample.sigma,
  cex = 0.5,
  pch = 16,
  col = col.alpha(rangi2, 0.1)
)





####
simps_2d <- function(f, ax,bx,nx, ay=ax,by=bx,ny=nx){
  #works for rectangular areas of integration
  #assumes function 'f' takes two arguments
  
  n1=2*nx+1 
  n2=2*ny+1      #number of points, including midway points
  
  h1=(bx-ax)/(n1-1)
  h2=(by-ay)/(n2-1) #length of subdivisions
  
  #create the Simpson matrix:
  s1=c(1,rep(2,n1-2)^(1:(n1-2)%%2 +1) , 1) 
  s2=c(1,rep(2,n2-2)^(1:(n2-2)%%2 +1) , 1)
  S=outer(s1,s2)
  
  #create the variable matrices
  xx=matrix(seq(ax,bx,length=n1), nrow=n1, ncol=n2)
  yy=matrix(seq(ay,by,length=n2), nrow=n1, ncol=n2, byrow=T)
  
  h1*h2*sum(S*f(xx,yy))/9 #compute the integral
}


sig <- 0.1
s_init <- 0
s_final <- 1

post2 <- function(mu, sigma){
  dlogitnorm(zz, mu, sigma)* dnorm(zz, mu, sig)* dunif(zz, s_init, s_final)
}

scale.factor <- simps_2d(post2, 0, 1, 100)

mu.est <- post$mu %>% mean()
sigma.est <- post$sigma %>% mean()

proportional <- dlogitnorm(zz, mu.est, sigma.est) * dnorm(zz, mu.est, sig) * dunif(zz, s_init, s_final)

prop <- post2(mu0, sigma0) / simps_2d(post2, 0, 1, 100, 1, 3, 100)




