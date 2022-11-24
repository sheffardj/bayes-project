library(tidyr)
library(rethinking)
source("seed.R")
# data that we're "given"
n <- 1000
zz <- rlogitnorm(n, 0.1, 2.1)
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

lh3 <- function(mu, sigma) {
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
  arr <- arr %>% arrange(rowid, colid)
  
  
  # plot layout
  layout(mat = matrix(c(1, 3, 2, 3),
                      nrow = 2,
                      ncol = 2),
         heights = c(1, 1),
         widths = c(1, 1))

  # These functions come from the 'rethinking' package
  contour_xyz(arr$rowid, arr$colid, arr$value, main=paste("Likelihood", func))
  image_xyz(arr$rowid, arr$colid, arr$value)
  
  # best pair (close but no cigar)
  mu.post <- arr[which.max(arr[,3]),1]
  sig.post <- arr[which.max(arr[,3]),2]
  
  zz %>% hist(., probability=T, ylim=c(0,2.5))
  curve(dlogitnorm(x, 0.1, 2.1), from=0, to=1, add=T, col='blue')
  curve(dlogitnorm(x, mu.post, sig.post), from=0, to=1, add=T, col='red')
  print(paste(mu.post, sig.post))
  
  
  assign("estimates", rbind(estimates, c(
    func = func,
    mu.post = mu.post,
    sig.post =  sig.post
  )), envir = .GlobalEnv)
  # fig <- plot_ly(z = ~ lh.est) %>% add_surface(contours = list(
  #   z = list(
  #     show = TRUE,
  #     usecolormap = TRUE,
  #     highlightcolor = "#ff0000",
  #     project = list(z = TRUE)
  #   )
  # ),
  # opacity = 0.6) %>% layout(scene = list(camera = list(eye = list(
  #   x = 1.87, y = 0.88, z = -0.64
  # ))))

  # htmlwidgets::saveWidget(widget = fig, file=paste0(func,".html"), selfcontained=TRUE)
}

# setup a storage matrix to store many Bayesian estimates
estimates <- matrix(NA, nrow=0, ncol=3) %>% as.data.frame()
colnames(estimates) <- c('lh', 'mu.post', 'sigma.post')

# estimate posterior, create plots, store estimates
estimate(lh1)
estimate(lh2)

# need to convert to numeric values
estimates <- estimates %>% mutate_at(c('mu.post', 'sigma.post'), as.numeric)

# > estimates
# lh    mu.post sigma.post
# 1 lh1 0.11821886   2.059281
# 2 lh2 0.08321886   2.059281

# we generate samples from our estimates
sample1 <- rlogitnorm(n, estimates[1,2],  estimates[1,3])
sample2 <- rlogitnorm(n, estimates[2,2],  estimates[2,3])

# store in matrix to compare with original logitnorm dat
dat1 <- rbind(zz/sum(zz),sample1/sum(sample1))
dat2 <- rbind(zz/sum(zz),sample2/sum(sample2))

# install.packages('philentropy')
library(philentropy)

KL(dat1) # kullback-leibler = 0.871625 
KL(dat2) # kullback-leibler = 0.8939065 

zz %>% hist(., probability=T, ylim=c(0,2.5))
curve(dlogitnorm(x, 0.1, 2.1), from=0, to=1, add=T, col='blue')
curve(dlogitnorm(x, estimates[1,2],  estimates[1,3]), from=0, to=1, add=T, col='red')
curve(dlogitnorm(x, estimates[2,2],  estimates[2,3]), from=0, to=1, add=T, col='green')
