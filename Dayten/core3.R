library(tidyr)
library(rethinking)
source("seed.R")
# data that we're "given"
dev.off()
zz <- rlogitnorm(1000, 0.1, 2.1)
zz %>% hist(., probability=T, ylim=c(0,1.9))

# data plotted against reality (blue) we want to match
curve(dlogitnorm(x, 0.1, 2.1), from=0, to=1, add=T, col='blue')

# we assume it is logit-normal, so we convert data to normal and get sample params
mu0 <- logit(zz) %>% mean()
sigma0 <- sqrt(logit(zz) %>% var())


# simpsons params
nx = ny = 100 #number of subdivisions
n1 = 2 * nx + 1 # length of mu sequence
n2 = 2 * ny + 1 # length of sigma sequence
ax = mu0 - 0.5 # range for mu
bx = mu0 + 0.5
ay = sigma0 - 0.5 # range for sigma
by = sigma0 + 0.5
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

# our likelihood functions
lh1 <- function(mu, sigma) {
  LL_obs <- mapply(dlogitnorm, mu=mu, sigma=sigma, MoreArgs = list(x=zz), SIMPLIFY = T)
  LL_mat <- matrix(apply(LL_obs, 2, prod), nrow=dim(mu)[1], ncol=dim(sigma)[2])
  LL_mat * dnorm(mu, mu0, 1) *
    dunif(sigma, 0, 3)
}

lh2 <- function(mu, sigma) {
  LL_obs <- mapply(dlogitnorm, mu=mu, sigma=sigma, MoreArgs = list(x=zz), SIMPLIFY = T)
  LL_mat <- matrix(apply(LL_obs, 2, prod), nrow=dim(mu)[1], ncol=dim(sigma)[2])
  LL_mat * dexp(mu, 1/mu0) *
    dnorm(sigma, sigma0, 3)
}


estimate <- function(lh){
  func <- as.character(substitute(lh))
  print(paste('func is', func))
  
  # running simpsons
  LH.MAT <- lh(xx, yy)
  scale <- h1 * h2 * sum(S * LH.MAT) / 9 #compute the integral
  lh.est <- LH.MAT/scale # return posterior
  
  
  #convert matrix to paired vectors
  rownames(lh.est) <- x_seq
  colnames(lh.est) <- y_seq
  M = lh.est %>% as.data.frame()
  M["rowid"] = row.names(M)
  arr <- gather(M, colid, value, -rowid)
  
  arr <- apply(arr, 2, as.numeric) %>% as.data.frame()
  # These functions come from the 'rethinking' package
  layout(mat = matrix(c(1, 3, 2, 3), 
                      nrow = 2, 
                      ncol = 2),
         heights = c(1, 1),    # Heights of the two rows
         widths = c(1, 1))     # Widths of the two columns
  
  
  contour_xyz(arr$rowid, arr$colid, arr$value, main=paste("Likelihood", func))
  image_xyz(arr$rowid, arr$colid, arr$value)
  
  # best pair (close but no cigar)
  mu.post <- arr[which.max(arr[,3]),1]
  sig.post <- arr[which.max(arr[,3]),2]
  
  zz %>% hist(., probability=T, ylim=c(0,2.5))
  curve(dlogitnorm(x, 0.1, 2.1), from=0, to=1, add=T, col='blue')
  curve(dlogitnorm(x, mu.post, sig.post), from=0, to=1, add=T, col='red')
  print(paste(mu.post, sig.post))
}


estimate(lh1)
estimate(lh2)
