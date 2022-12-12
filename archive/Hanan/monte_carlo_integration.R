library(tidyr)

source("seed.R")
# data that we're "given"
zz <- rlogitnorm(1000, 0.1, 2.1)
zz %>% hist(., probability=T, ylim=c(0,2))

# we assume it is logit-normal, so we convert data to normal and get sample params
mu0 <- logit(zz) %>% mean()
sigma0 <- sqrt(logit(zz) %>% var())

# Take a look the prior distribution with the estimated mu and sigma vectors
#prior_zz <- rlogitnorm(1e4, sample_mu, sample_sigma)
#prior_zz %>% hist(., probability=T)

# data plotted against reality (blue) we want to match
curve(dlogitnorm(x, 0.1, 2.1), from=0, to=1, add=T, col='blue', n=301)

# monte carlo integration set up

n1 = 201 # number of observations
n2 = n1 # length of sigma sequence
ax = mu0 - 0.1 # range for mu
bx = mu0 + 0.1
ay = sigma0 - 0.1 # range for sigma
by = sigma0 + 0.1
D = (bx - ax)*(by - ay)  #area of integration

# generate random observations
x_seq <- sort(runif(n1, ax, bx))
y_seq <- sort(runif(n2, ay, by))

# create the variable matrices
xx=matrix(x_seq, nrow=n1, ncol=n2)
yy=matrix(y_seq, nrow=n1, ncol=n2, byrow=T)

# our likelihood function
lh <- function(mu, sigma) {
  LL_obs <- mapply(dlogitnorm, mu=mu, sigma=sigma, MoreArgs = list(x=zz), SIMPLIFY = T)
  LL_mat <- matrix(apply(LL_obs, 2, prod), nrow=dim(mu)[1], ncol=dim(sigma)[2])
  LL_mat * dnorm(mu, 1, 0.1) *
    dunif(sigma, 1, 3)
}

LH.MAT <- lh(xx,yy)
scale <- D * mean(LH.MAT) #compute the integral
lh.est <- LH.MAT/scale # return posterior
lh.est %>% image()

#convert matrix to paired vectors
rownames(lh.est) <- x_seq
colnames(lh.est) <- y_seq
M = lh.est %>% as.data.frame()
M["rowid"] = row.names(M)
arr <- gather(M, colid, value, -rowid)

arr <- apply(arr, 2, as.numeric) %>% as.data.frame()
# These functions come from the 'rethinking' package
contour_xyz(arr$rowid, arr$colid, arr$value)
image_xyz(arr$rowid, arr$colid, arr$value)

# best pair (close but no cigar)
mu.post <- arr[which.max(arr[,3]),1]
sig.post <- arr[which.max(arr[,3]),2]

zz %>% hist(., probability=T,ylim=c(0,2))
curve(dlogitnorm(x, 0.1, 2.1), from=0, to=1, add=T, col='blue', n=301)
curve(dlogitnorm(x, mu.post, sig.post), from=0, to=1, add=T, col='red', n=301)
