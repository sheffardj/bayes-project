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
curve(dlogitnorm(x, 0.1, 2.1), from=0, to=1, add=T, col='blue',n=301)

# simpsons params
nx = ny = 100 #number of subdivisions
n1 = 2 * nx + 1 # length of mu sequence
n2 = 2 * ny + 1 # length of sigma sequence
ax = mu0 - 0.1 # range for mu
bx = mu0 + 0.1
ay = sigma0 - 0.1 # range for sigma
by = sigma0 + 0.1
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
  LL_obs <- mapply(dlogitnorm, mu=mu, sigma=sigma, MoreArgs = list(x=zz), SIMPLIFY = T)
  LL_mat <- matrix(apply(LL_obs, 2, prod), nrow=dim(mu)[1], ncol=dim(sigma)[2])
  LL_mat * dnorm(mu, 1, 0.1) *
    dunif(sigma, 1, 3)
}

# ll_mat <- 1;
# for (i in 1:length(zz)) {
#   ll_mat = ll_mat*dnorm(zz[i], xx, yy)
# }
# ll_mat <- ll_mat*dnorm(xx,1,0.1)*dunif(yy,0,3)
# ll_scaled <- ll_mat / (h1 * h2 * sum(S * ll_mat) / 9)
# persp(x_seq, y_seq,z = ll_scaled)

# lh <- function(param){
#   ll_prod <- rep(1,length(param))
#   for (i in 1:length(param)) {
#     ll_prod = ll_prod*(dlogitnorm(zz, param[i,1], param[i,2])) * dnorm(param[i,1], 1, 0.1) *
#       dunif(param[i,2], 0, 3)
#   }
#   ll_prod
# }

# running simpsons
LH.MAT <- lh(xx, yy)
scale <- h1 * h2 * sum(S * LH.MAT) / 9 #compute the integral
lh.est <- LH.MAT/scale # return posterior
#lh.est <- lh(expand.grid(x_seq,y_seq))/scale # return posterior
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

zz %>% hist(., probability=T, ylim = c(0,2))
curve(dlogitnorm(x, mu.post, sig.post), from=0, to=1, add=T, col='red', n=301)

# # Looking at LOG likelihood
# llh <- function(mu, sigma) {
#   dlogitnorm(zz, mu, sigma, log = TRUE) + dnorm(zz, mu, 0.1, log = TRUE) +
#     dunif(zz, 0, sigma + 1, log = TRUE)
# }
# 
# llh.est <- dlogitnorm(zz, mu.est, sigma.est, log = TRUE) +
#   dnorm(zz, mu.est, 0.1, log = TRUE) +
#   dunif(zz, 0, 3, log = TRUE)
# 
# scale <- h1 * h2 * sum(S * llh(xx, yy)) / 9 #compute the integral
# exp((llh(xx, yy) / (-scale))) %>% image()
# 
