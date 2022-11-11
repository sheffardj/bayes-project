# https://link.springer.com/content/pdf/10.1007/978-0-387-71385-4.pdf
install.packages('LearnBayes')
library(LearnBayes)

#### CHAPTER 2 ####
###### DISCRETE PRIOR (Section 2.1 - 2.3, page 29, book page 19) ######
# sleep study ( >= 8 hours is a success)
p = seq(0.05, 0.95, by = 0.1) # initial beliefs what values acceptable
prior = c(2, 4, 8, 8, 4, 2, 1, 1, 1, 1) # weights for initial beliefs
prior = prior/sum(prior) # get prior proportionally
plot(p, prior, type = "h", ylab="Prior Probability")

data = c(11, 16) # 11 slept >= 8 hours, 16 slept < 8 hours
post = pdisc(p, prior, data) # generates the posterior
cbind(p, prior, post)
plot(p, post, type = "h", ylab="Posterior Probability")

###### CONTINUOUS PRIOR (section 2.4, page 32/22 pdf/book) ######
# beta => conjugate - prior and posterior have same functional form
p = seq(0, 1, length = 500)
a = 3.4 # a,b found by trial & error (90th percentiles given by .3 and .5)
b = 7.4
s = 11 # 11 slept >= 8hrs
f = 16 # 16 did not

prior=dbeta(p,a,b)
like=dbeta(p,s+1,f+1)
post=dbeta(p,a+s,b+f)

plot(p,post,type="l",ylab="Density",lty=2,lwd=3)
lines(p,like,lty=1,lwd=3)
lines(p,prior,lty=3,lwd=3)
legend(.7,4,c("Prior","Likelihood","Posterior"), lty=c(3,1,2),lwd=c(3,3,3))

# probability 50% of sleepers get >= 8 hours
1 - pbeta(0.5, a + s, b + f)
#0.0684256943

# 95% interval estimate
qbeta(c(0.025, 0.975), a + s, b + f)
#       2.5%       97.5% 
#[1] 0.235099648 0.538742334

###### SIMULATING FOR POSTERIOR (section 2.4, page  32/22) ######
ps = rbeta(1000, a + s, b + f)
hist(ps,xlab="p",main="")
sum(ps >= 0.5)/1000 #very close
#[1] 0.075

quantile(ps, c(0.025, 0.975))
#       2.5%       97.5% 
#0.244865531 0.539359586 




###### GRID APPROXIMATION (section 2.5, page 36/26) ######
midpt = seq(0.05, 0.95, by = 0.1) # midpoints of desired interval/grid
prior = c(2, 4, 8, 8, 4, 2, 1, 1, 1, 1) # weights at midpt
prior = prior/sum(prior) # get the prior proportionally

p = seq(0, 1, length = 500) # for a smooth function
plot(p, histprior(p, midpt, prior), type = "l",
       ylab = "Prior density", ylim = c(0, .25))

like = dbeta(p, s + 1, f + 1)
post = like * histprior(p, midpt, prior)
plot(p, post, type = "l",ylab="Posterior density") # updated likelihood * prior

post = post/sum(post)
ps = sample(p, replace = TRUE, prob = post)
hist(ps, xlab="p")

###### PREDICTIONS (section 2.6, page 39/29) ######
# Using discrete - BINOMIAL
p = seq(0.05, 0.95, by = .1) # considered values
prior = c(2, 4, 8, 8, 4, 2, 1, 1, 1, 1) # weights
prior = prior / sum(prior) # proportional
m = 20 # future sample size
ys = 0:20 # number of successes of interest
pred = pdiscp(p, prior, m, ys) #
cbind(0:20, pred) # 5 and 6 MOST LIKELY

# Using continuous BETA
ab=c(3.4, 7.4)
m=20; ys=0:20
pred=pbetap(ab, m, ys)
plot(1:21, pred, type = "h", ylab="") # 6 and 7 MOST LIKELY

# Simulation

p = rbeta(1000, 3.4, 7.4)
y = rbinom(1000, 20, p)
table(y)
freq = table(y)
ys = c(0:max(y))
predprob = freq / sum(freq)
plot(ys, predprob, type = "h", xlab = "y") # 4 and 7 MOST LIKELY

dist=cbind(ys,predprob)
dist
covprob=.95
discint(dist, covprob) #1-13 successes is the 95% interval


#### CHAPTER 3 ####