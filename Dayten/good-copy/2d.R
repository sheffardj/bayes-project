start <- Sys.time()
source("seed.R")
source(paste0(getwd(),"/good-copy/prior-post.R")) # load in priors and posterior
source(paste0(getwd(),"/good-copy/case.R")) # import function to setup simpsons
source(paste0(getwd(),"/good-copy/estimate.R")) # import our estimation function

# number of points to simulate from logitnorm data
n <- 30

# all pairs of parameters
parms <-cbind(c(rep(0,3), rep(1,3), rep(2, 3)), rep(c(3.16, 1.78, 0.32),3))
colnames(parms) <- c("mu", "sigma")

# ids to track which priors pair we're using
run_ids <- cbind(c(rep(1,3), rep(2,3), rep(3, 3)),rep(1:3,3))
colnames(run_ids) <- c('pi', 'pj')

# setup a storage matrix to store many Bayesian estimates
estimates <- matrix(NA, nrow=0, ncol=4) %>% as.data.frame()
colnames(estimates) <- c('case', "run_id", 'mu.post', 'sigma.post')

grid_size <- 30
# place to store posteriors (grid values) for heatmaps
arr <- matrix(nrow = (2 * grid_size + 1)^2, ncol=0) 

# estimate posterior, create plots, store estimates
for(case_num in 1:dim(parms)[1]){
  print(paste("case", case_num))
  
  # setup simpsons integration matrix per param pair
  case(n, case_num, parms[case_num, 1], parms[case_num, 2], grid_size) 
  
  # run each pair of priors
  for(ii in 1:dim(run_ids)[1]){
    estimate(case_num, ii,
             get(paste0("p", run_ids[ii, 1])), # to call each pair of priors
             get(paste0("p", run_ids[ii, 2])))
  }
}
beep(sound = 4) # alert to alg completion

 # need to convert to numeric values
colnames(estimates) <- c('case', "run_id", 'mu.post', 'sigma.post')
# estimated values not numeric due to an earlier transformation
estimates <- estimates %>% mutate_at(c('mu.post', 'sigma.post'), as.numeric)
estimates$pi <- run_ids[,1] # identify pairs of priors pi
estimates$pj <- run_ids[,2] # and pj

# just aligns the underlying mu/sigma pairs
estimates <- estimates %>% 
  mutate(mu = c(rep(0,27), rep(1,27), rep(2,27))) %>% 
  mutate(sigma = rep(c(rep(3.16,9), rep(1.78, 9), rep(0.32, 9)),3))

estimates

estimates$kl <- NA
estimates$map_mu <- NA
estimates$map_sig <- NA
estimates$mu_low <- NA
estimates$mu_hi <- NA
estimates$sig_low <- NA
estimates$sig_hi <- NA
estimates$pd_mu <- NA
estimates$ps_mu <- NA
estimates$pd_sig <- NA
estimates$ps_sig <- NA


end <- Sys.time()
end - start



if(0){
### BELOW JUST SOME PLOTS
par(mfrow=c(3,3))
for(jj in 1:dim(parms)[1]){
  datj <- get(paste0('dat',jj))
  est.mu <- logit(datj) %>% mean()
  est.sig <- sqrt(logit(datj) %>% var())
                    
  hist(get(paste0("dat",jj)), probability = T,
       main=
         bquote(mu[est]~" = "~.(est.mu)~" "~sigma[est]~" ="~.(est.sig)),
       xlab="",
       ylim=c(0,12))
  est <- estimates %>% filter(case == jj)
  for(ii in 1:dim(run_ids)[1]){
    curve(
      dlogitnorm(x, mu = est[ii, "mu.post"], sigma = est[ii, "sigma.post"]),
      from = 0,
      to = 1,
      add
      = T,
      col = ii,
      lty=2,
      n = 300
    )
  }
  curve(
    dlogitnorm(x, mu = parms[jj, "mu"], sigma = parms[jj, "sigma"]),
    from = 0,
    to = 1,
    add = T,
    col = 'red',
    n = 300,
  )
}

par(mfrow=c(3,3))
for(ii in 1:9){
  tmp <- arr[,1:3]
   image_xyz(tmp$mu, tmp$sigma, tmp$prob, main=paste("Likelihood", ii))
}
}

for(ii in 1:dim(estimates)[1]){
  mu.post <- estimates[ii,"mu.post"]
  sigma.post <- estimates[ii,"sigma.post"]
  
  # generate samples from our estimates
  samp <- rlogitnorm(n, mu.post,  sigma.post)
  
  # store in matrix to compare with original logitnorm dat
  tmp <- rbind(og_data/sum(og_data),samp/sum(samp))
  
  KL(tmp) # kullback-leibler = 0.8939065 
  estimates[ii,"kl"] <- KL(tmp)
  arr0 <- arr[, (3 * ii - 2):(3 * ii)]
  
  sample.rows <-sample(1:nrow(arr0),size=1e3,replace=TRUE,
                       prob=arr0$prob )
  sample.mu <-arr0$mu[sample.rows]
  sample.sigma <-arr0$sigma[sample.rows]
  
  post <- tibble(sample.mu,
                 sample.sigma)
  
  CI <- describe_posterior(
    post,
    centrality = "MAP",
    test = c("p_direction", "p_significance")
  )
  
  estimates[ii,'map_mu'] <- CI$MAP[1]
  estimates[ii,'map_sig'] <- CI$MAP[2]
  estimates[ii,'mu_low'] <- CI$CI_low[1]
  estimates[ii,'mu_hi'] <-  CI$CI_high[2]
  estimates[ii,'sig_low'] <-  CI$CI_low[2]
  estimates[ii,'sig_hi'] <-  CI$CI_high[2]
  estimates[ii,'pd_mu'] <- CI$pd[1]
  estimates[ii,'ps_mu'] <- CI$ps[1]
  estimates[ii,'pd_sig'] <- CI$pd[2]
  estimates[ii,'ps_sig'] <- CI$ps[2]
}

estimates %>% as.data.frame() %>%  arrange(kl)
beep(sound = 2) # alert to alg completion

fwrite(estimates, file='estimates.csv')
