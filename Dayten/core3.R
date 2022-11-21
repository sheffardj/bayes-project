library(tidyr)

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

# simpsons params
nx = ny = 1000 #number of subdivisions
n1 = 2 * nx + 1 # length of mu sequence
n2 = 2 * ny + 1 # length of sigma sequence
ax = mu0 - 2 # range for mu
bx = mu0 + 2
ay = sigma0 - 2 # range for sigma
by = sigma0 + 2
h1 = (bx - ax) / (n1 - 1)  #length of subdivisions
h2 = (by - ay) / (n2 - 1)

#create the Simpson matrix:
s1 = c(1, rep(2, n1 - 2) ^ (1:(n1 - 2) %% 2 + 1) , 1)
s2 = c(1, rep(2, n2 - 2) ^ (1:(n2 - 2) %% 2 + 1) , 1)
S = outer(s1, s2)

# set out the sequences (used elsewhere)
x_seq <- seq(ax,bx,length=n1)
y_seq <- seq(ay,by,length=n2)

# create the variable matrices
xx=matrix(x_seq, nrow=n1, ncol=n2)
yy=matrix(y_seq, nrow=n1, ncol=n2, byrow=T)

# our likelihood function
lh <- function(mu, sigma) {
  dlogitnorm(zz, mu, sigma) * dnorm(mu, mu0, 1) *
    dunif(sigma, 0, sigma0 + 1)
}

# running simpsons
scale <- h1 * h2 * sum(S * lh(xx, yy)) / 9 #compute the integral
lh.est <- (lh(xx, yy)/scale) # return posterior
lh.est %>% image()




#convert matrix to paired vectors
rownames(lh.est) <- x_seq
colnames(lh.est) <- y_seq
M = lh.est %>% as.data.frame()
M["rowid"] = row.names(M)
arr <- gather(M, colid, value, -rowid)

arr[,1:2] %>% image()

arr <- apply(arr, 2, as.numeric) %>% as.data.frame()
# These functions come from the 'rethinking' package
contour_xyz(arr$rowid, arr$colid, arr$value)
image_xyz(arr$rowid, arr$colid, arr$value)

# best pair (close but no cigar)
mu.post <- arr[which.max(arr[,3]),1]
sig.post <- arr[which.max(arr[,3]),2]

curve(dlogitnorm(x, mu.post, sig.post), from=0, to=1, add=T, col='red')

arr %>% arrange(desc(value)) %>% head(10)
# # Looking at LOG likelihood
llh <- function(mu, sigma) {
  dlogitnorm(zz, mu, sigma, log = TRUE) + dnorm(mu, mu0, 0.1, log = TRUE) +
    dunif(sigma, 0, sigma0 + 1, log = TRUE)
}

scale <- h1 * h2 * sum(S * llh(xx, yy)) / 9 #compute the integral
exp((llh(xx, yy) / (-scale))) %>% image()

