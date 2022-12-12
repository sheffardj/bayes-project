start <- Sys.time()
source("seed.R")
source(paste0(getwd(),"/good-copy/prior-post.R")) # load in priors and posterior
source(paste0(getwd(),"/good-copy/case.R")) # import function to setup simpsons
source(paste0(getwd(),"/good-copy/estimate.R")) # import our estimation function
custPalette <- c("#002754", "#005493",  "#C63527", "#F5AA1C")

# number of points to simulate from logitnorm data
n <- 30

# all pairs of parameters
parms <- cbind(c(rep(0,3), rep(1,3), rep(2, 3)), rep(c(3.16, 1.78, 0.32),3))
colnames(parms) <- c("mu", "sigma")

# ids to track which priors pair we're using
run_ids <- cbind(c(rep(1,3), rep(2,3), rep(3, 3)),rep(1:3,3))
colnames(run_ids) <- c('pi', 'pj')

# setup a storage matrix to store many Bayesian estimates
estimates <- matrix(NA, nrow=0, ncol=4) %>% as.data.frame()
colnames(estimates) <- c('case', "run_id", 'mu.post', 'sigma.post')

grid_size <- 100
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


# RESULTS PLOT
if(0){
  png(file="2d-results.png",
      width=1800, height=1200, units="px", res=190, bg='white')
  par(mfrow=c(3,3), mar=c(3,4.1,4,2))
  for(jj in 1:dim(parms)[1]){
    datj <- get(paste0('dat',jj))
    est.mu <- logit(datj) %>% mean()
    est.sig <- sqrt(logit(datj) %>% var())
                      
    max <- max(dlogitnorm(seq(0, 1, by = 0.01), 
                          mu = parms[jj, "mu"], 
                          sigma = parms[jj, "sigma"]))
    print(max)
    hist(get(paste0("dat",jj)), probability = T,
         main=
           bquote(hat(mu)~" = "~.(round(est.mu,3))~" "~hat(sigma)~" ="~.(round(est.sig,3))),
         xlab="",
         breaks=5,
         col='lightblue',
         ylim=c(0, max + ifelse(jj %in% c(2,3,8,9), 1, 0) + ifelse( jj == 9, 2, 0)))
    est <- estimates %>% filter(case == jj)
    for(ii in 1:dim(run_ids)[1]){
      curve(
        dlogitnorm(x, mu = est[ii, "mu.post"], sigma = est[ii, "sigma.post"]),
        from = 0,
        to = 1,
        add
        = T,
        col = rev(hcl.colors(12, 'Viridis'))[ii],
        lty=c(6,4,2,3)[ii %% 4 + 1],
        n = 300
      )
    }
    curve(
      dlogitnorm(x, mu = parms[jj, "mu"], sigma = parms[jj, "sigma"]),
      from = 0,
      to = 1,
      add = T,
      col = custPalette[3],
      n = 300,
    )
  }
  dev.off()
}
# POSTERIOR HEATMAP FOR ONE CASE
if(0){
  png(file="2d-posterior.png",
      width=1800, height=1200, units="px", bg='white')
  labels <- c("Uniform", "Normal", "Laplace")
  par(mfrow=c(3,3), mar=c(4,4,2,0.5))
  for(ii in 10:18){
    priors <- run_ids[ii - 9,]
    tmp <- arr[,(ii*3-2):(3*ii)] %>% arrange(mu, sigma)
    image_xyz(x=tmp$mu, y=tmp$sigma, z=tmp$prob,
              xlab=bquote(mu["est"]), ylab=bquote(sigma["est"]), main=paste("Priors:",labels[priors[1]], "x", labels[priors[2]]))
  }
  dev.off()
}

# 2D ESTIMATES GT TABLE
if(0){
for(ii in 1:dim(estimates)[1]){
  mu.post <- estimates[ii,"mu.post"]
  sigma.post <- estimates[ii,"sigma.post"]
  
  # generate samples from our estimates
  samp <- rlogitnorm(n, mu.post,  sigma.post)
  
  # store in matrix to compare with original logitnorm dat
  tmp <- rbind(og_data/sum(og_data),samp/sum(samp))
  
  # store KL divergence in estimates table
  estimates[ii,"kl"] <- KL(tmp)
  
  # start CREDIBLE intervals
  arr0 <- arr[, (3 * ii - 2):(3 * ii)]
  
  #sample from mu and sigma by probability
  sample.rows <-sample(1:nrow(arr0),size=1e3,replace=TRUE,
                       prob=arr0$prob )
  sample.mu <-arr0$mu[sample.rows]
  sample.sigma <-arr0$sigma[sample.rows]
  
  # arrange sample data for credible interval computation
  post <- tibble(sample.mu, sample.sigma)
  
  # retrieve CI and set relative values
  CI <- describe_posterior(
    post,
    centrality = "MAP",
    test = c("p_direction", "p_significance"),
    ci=0.997
  )
  
  estimates[ii,'map_mu'] <- CI$MAP[1]
  estimates[ii,'map_sig'] <- CI$MAP[2]
  estimates[ii,'mu_low'] <- CI$CI_low[1]
  estimates[ii,'mu_hi'] <-  CI$CI_high[1]
  estimates[ii,'sig_low'] <-  CI$CI_low[2]
  estimates[ii,'sig_hi'] <-  CI$CI_high[2]
  estimates[ii,'pd_mu'] <- CI$pd[1]
  estimates[ii,'ps_mu'] <- CI$ps[1]
  estimates[ii,'pd_sig'] <- CI$pd[2]
  estimates[ii,'ps_sig'] <- CI$ps[2]
}

estimates %>% as.data.frame() %>%  arrange(kl)
beep(sound = 2) # alert to alg completion

# fwrite(estimates, file='estimates.csv')
# estimates<-fread("estimates.csv")

estimates %>% 
  group_by(case) %>% 
  slice_min(order_by = kl, n=1) -> best_cases

best_cases %>% 
  select(case, mu, sigma, pi, pj, mu.post, sigma.post, kl,
         map_mu, map_sig, mu_low, mu_hi, sig_low, sig_hi) -> best_cases

best_cases[,c(3,6:7,9:14)] <- apply(best_cases[,c(3,6:7,9:14)], 2, function(x) round(x,2))
best_cases[,8] <- apply(best_cases[,8], 2, function(x) round(x,4))
best_cases[,4:5] <- apply(best_cases[,4:5], 2, function(x) c("U", "N", "L")[x])

gt(ungroup(best_cases)) %>% 
  tab_header(title = "Best Estimate Per Case") %>% 
  tab_stubhead(label = "Case") %>%
  cols_label(
    case= "Case",
    mu= html("&mu;"),
    sigma= html("&sigma;"),
    mu.post = html("&mu;<sub>post</sub>"),
    sigma.post = html("&sigma;<sub>post</sub>"),
    kl="KL",
    map_mu = html("&mu;<sub>MAP</sub>"),
    map_sig = html("&sigma;<sub>MAP</sub>"),
    mu_low=html("&mu;<sub>low</sub>"),
    mu_hi=html("&mu;<sub>high</sub>"),
    sig_low=html("&sigma;<sub>low</sub>"),
    sig_hi=html("&sigma;<sub>high</sub>"),
  ) -> base_table

# add spanners
base_table <- base_table %>% 
  tab_spanner(label = md("*Configuration p(x)*"),
              columns = c(case, mu, sigma, pi, pj)) %>% 
  tab_spanner(label = md("*Posteriors q(x)*"),
              columns = c(mu.post, sigma.post)) %>% 
  tab_spanner(label = md("*Score*"),
              columns = c(kl))%>% 
  tab_spanner(label = md("*.997 Credibility Interval*"),
              columns = c(map_mu, map_sig, mu_low, mu_hi, sig_low, sig_hi))

  
  
  
# add colors  
base_table <- base_table  %>%
  data_color(
    columns = c(mu.post, sigma.post),
    colors = scales::col_numeric(
      palette = paletteer::paletteer_d(
        palette = "ggsci::red_material"
      ) %>%
        as.character(),
      domain = NULL
    )
  ) %>%
  data_color(
    columns = c(kl),
    colors = scales::col_numeric(
      palette = paletteer::paletteer_d(
        palette = "ggsci::green_material"
      ) %>%
        as.character(),
      domain = NULL
    )
  ) %>% 
  data_color(
    columns = c(map_mu, map_sig, mu_low, mu_hi, sig_low, sig_hi),
    colors = scales::col_numeric(
      palette = paletteer::paletteer_d(
        palette = "ggsci::yellow_material"
      ) %>%
        as.character(),
      domain = NULL
    )
  ) %>% 
  data_color(
    columns = c(case, mu, sigma),
    colors = scales::col_numeric(
      palette = paletteer::paletteer_d(
        palette = "ggsci::blue_material"
      ) %>%
        as.character(),
      domain = NULL
    )
  ) %>% 
  data_color(
    columns = c(pi,pj),
    colors = scales::col_factor(
      palette = paletteer::paletteer_d(
        palette = "ggsci::purple_material"
      ) %>%
        as.character(),
      domain = NULL
    )
  ) 
base_table %>% 
  cols_width(everything() ~ px(60)) -> table_plt
table_plt
}