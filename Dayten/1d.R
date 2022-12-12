##### PREAMBLE #####
set.seed(8776)
library(logitnorm)
library(sn)
library(dplyr)
library(philentropy)
library(bayestestR)
library(ggplot2)
library(ggpubr)
custPalette <- c("#002754", "#005493",  "#C63527", "#F5AA1C")

##### PRIORS #####

un <- function (d, v) {
  dunif(d, v-3,v+3 )
}
nor <- function (d, v){
  dnorm(d, 0, 3)
}
lp <- function(x, mean, sd = 3) {
  # Laplace (double exponential) density function with mean equal to...
  # ... \code{mean} and standard deviation equal to \code{sd}. 
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
post_length <- 3 * 300 + 1 #based on simpsons 3/8
data <- matrix(nrow = n, ncol = 9) # storage for original data
posteriors <- matrix(nrow = post_length, ncol = 54) # storage for found posteriors
xx_seqs <- matrix(nrow = post_length, ncol = 54)
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
  var0 <- ifelse(fixed=="mu", sd(logit(dat)), mean(logit(dat)))
  
  # set up simpsons params
  a <- ifelse(fixed=="mu", max(var0 - 1, 0.01), var0 - 1) # toggled on `fixed`
  b <- var0 + 1
  nn = post_length
  h = (b - a) / (nn - 1) #length of subdivision
  xx <- seq(a, b, length=nn)
  
  # obtain likelihood
  if(fixed == "mu"){
    # Likelihood with mu fixed
    LL_obs <- mapply(dlogitnorm, mu=fixed_val, sigma=xx, 
                     MoreArgs = list(x=dat), SIMPLIFY = T)
  } else {
    # Likelihood with sigma fixed
    LL_obs <- mapply(dlogitnorm, mu=xx, sigma=fixed_val, 
                     MoreArgs = list(x=dat), SIMPLIFY = T)
  }
  LL <- abs(apply(LL_obs, 2, prod))
  
  # obtain
  y <- LL*prior(xx, var0)
  
  # compute the denominator with Simpsons rule
  dnom <-
    h * (y[1] + y[nn] + 3 * sum(y[(2:(nn - 1))[2:(nn - 1) %% 3 != 1]])  + 2 *
           sum(y[(2:(nn - 1))[2:(nn - 1) %% 3 == 1]])) / 8
  
  # compute the posterior and assign to Global Env,.
  posterior <- c(y/dnom)
  
  # finds the best estimate and stores to global env.
  best_est_index <- which.max(posterior)
  best_est <- xx[best_est_index]
  
  return(list(dat = dat, posterior = posterior, best_est=best_est, xx_seq = xx))
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
  xx_seqs[,ii] <- results$xx_seq
  if(ii %% 6 == 0) cat(paste0(round(100*ii/dim(cases))[1], '% '))
}



##### PLOTS #####
if(0){
### MU FIXED
par(mfrow=c(3,3), mar=c(3,4,3,1))
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
    
    sig.hat <- round(mean(best_un, best_nor, best_lp), 3) # just a plotting ref
    colors <- rev(hcl.colors(12, 'Viridis'))[c(1,5,10)]
    
    # this gets a nice `ymax` since any curve or histogram prob can be the max
    hist <- hist(dat, plot=F)
    
    dens_og  <- dlogitnorm(seq(0, 1, length = 60), ii, jj)
    dens_un  <- dlogitnorm(seq(0, 1, length = 60), ii, best_un)
    dens_nor <- dlogitnorm(seq(0, 1, length = 60), ii, best_nor)
    dens_lp  <- dlogitnorm(seq(0, 1, length = 60), ii, best_lp)
    ymax <- max(c(hist$density, dens_og, dens_un, dens_nor, dens_lp))
    print(ymax)
    
    hist(dat, probability = T, ylim=c(0,ymax+1), xlim=c(0,1),
         main=bquote(mu~" = "~.(ii)~" "~sigma~" = "~.(jj)~", "~
                       bar(sigma)~" = "~.(sig.hat)),
         col='lightblue')
    curve(dlogitnorm(x, ii, jj), col = custPalette[3], add = T)
    curve(dlogitnorm(x, ii, best_un),  col =colors[1], lty=6, add = T)
    curve(dlogitnorm(x, ii, best_nor), col =colors[2], lty=4, add = T)
    curve(dlogitnorm(x, ii, best_lp),  col =colors[3], lty=3, add = T)
  }
}


### SIGMA FIXED
par(mfrow=c(3,3), mar=c(3,4,3,1))
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
    
    mu.hat <- round(mean(best_un,best_nor,best_lp), 3) # just a ref for plotting
    colors <- rev(hcl.colors(12, 'Viridis'))[c(1,5,10)]
    
    # this gets a nice `ymax` since any curve or histogram prob can be the max
    hist <- hist(dat, plot=F)
    dens_og  <- dlogitnorm(seq(0,1,length=60), ii, jj)
    dens_un  <- dlogitnorm(seq(0,1,length=60), best_un, jj)
    dens_nor <- dlogitnorm(seq(0,1,length=60), best_nor, jj)
    dens_lp  <- dlogitnorm(seq(0,1,length=60), best_lp, jj)
    ymax <- max(c(hist$density, dens_og, dens_un, dens_nor, dens_lp))
    print(ymax)
    
    hist(dat, probability = T, ylim=c(0, ymax + 1), xlim=c(0,1),
         main=bquote(atop(mu~" = "~.(ii)~" "~sigma~" = "~.(jj)~", "~
                            bar(mu)~" = "~.(mu.hat))),
         col='lightblue')
    curve(dlogitnorm(x, ii, jj), col = custPalette[3], add = T)
    curve(dlogitnorm(x, best_un, jj),  col =colors[1], lty=6, add = T)
    curve(dlogitnorm(x, best_nor, jj), col =colors[2], lty=4, add = T)
    curve(dlogitnorm(x, best_lp, jj),  col =colors[3], lty=3, add = T)
  }
}
}


##### DIVERGENCE #####
cases$kl <- NA
cases$map_est <- NA
cases$est_low <- NA
cases$est_hi <- NA
cases$pd <- NA
cases$ps <- NA



for(ii in 1:dim(cases)[1]){
  est <- cases[ii,"best_est"] 
  
  case_id <- base_cases %>% 
    filter(mu == cases[ii,"mu"] & sigma == cases[ii,"sigma"]) %>% 
    pull(case_id)
  dat <- data[, case_id]
  
  if(cases[ii,"fixed"] == 'mu'){
    samp <- rlogitnorm(n, cases[ii,"mu"],  est)
  } else {
    samp <- rlogitnorm(n, est,  cases[ii,"sigma"])
  }
  # generate samples from our estimates
  
  # store in matrix to compare with original logitnorm dat
  tmp <- rbind(dat/sum(dat),samp/sum(samp))
  
  # store KL divergence in estimates table
  cases[ii,"kl"] <- KL(tmp)
  
  # start CREDIBLE intervals
  post <- posteriors[,ii]
  xx <- xx_seqs[,ii]
  #sample from mu and sigma by probability
  sample.rows <-sample(1:length(post),size=1e3,replace=TRUE,prob=post)
  sample.est <- xx[sample.rows]
  
  # arrange sample data for credible interval computation
  
  # retrieve CI and set relative values
  CI <- describe_posterior(
    sample.est,
    centrality = "MAP",
    test = c("p_direction", "p_significance")
  )
  
  cases[ii,'map_est'] <- CI$MAP
  cases[ii,'est_low'] <- CI$CI_low
  cases[ii,'est_hi'] <- CI$CI_high
  cases[ii,'pd'] <- CI$pd
  cases[ii,'ps'] <- CI$ps
}

cases %>% str

cases <- cases %>% 
  mutate(prior = case_when(prior == 'un' ~ 1,
                           prior == 'nor' ~ 2,
                           TRUE ~ 3))

##### TABLE FOR MU FIXED #####
if(0){
cases %>% 
  filter(fixed=='mu') %>% 
  group_by(mu, sigma) %>% 
  slice_min(order_by = kl, n=1) -> best_cases
best_cases[,c(5,7:11)] <- apply(best_cases[,c(5,7:11)], 2, function(x) round(x,3))
best_cases[,6] <- apply(best_cases[,6], 2, function(x) round(x,4))
best_cases[,3] <- apply(best_cases[,3], 2, function(x) c("U", "N", "L")[x])
  
best_cases$case <- 1:9
best_cases <- best_cases[c(12, 1:3, 5:9)]


gt(ungroup(best_cases)) %>% 
  tab_header(title = html("Best Estimate Per Case in 1d (&mu; fixed)")) %>% 
  tab_stubhead(label = "Case") %>%
  cols_label(
    case= "Case",
    mu= html("&mu;<sub>fixed</sub>"),
    sigma= html("&sigma;"),
    prior= html("Prior"),
    best_est = html("&sigma;<sub>est</sub>"),
    kl="KL",
    map_est = html("&sigma;<sub>MAP</sub>"),
    est_low= html("CI<sub>low</sub>"),
    est_hi=html("CI<sub>high</sub>"),
  ) -> base_table

# add spanners
base_table <- base_table %>% 
  tab_spanner(label = md("*Configuration*"),
              columns = c(case, mu, sigma, prior)) %>% 
  tab_spanner(label = md("*Post*"),
              columns = c(best_est)) %>% 
  tab_spanner(label = md("*Score*"),
              columns = c(kl))%>% 
  tab_spanner(label = md("*Credibility Interval*"),
              columns = c(map_est, est_low, est_hi))


# add colors  
base_table <- base_table  %>%
  data_color(
    columns = c(best_est),
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
    columns = c(map_est, est_low, est_hi),
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
    columns = c(prior),
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


if(0){
  cases %>% 
    filter(fixed=='sigma') %>% 
    group_by(mu, sigma) %>% 
    slice_min(order_by = kl, n=1) -> best_cases
  best_cases[,c(5,7:11)] <- apply(best_cases[,c(5,7:11)], 2, function(x) round(x,3))
  best_cases[,6] <- apply(best_cases[,6], 2, function(x) round(x,4))
  best_cases[,3] <- apply(best_cases[,3], 2, function(x) c("U", "N", "L")[x])
  
  best_cases$case <- 1:9
  best_cases <- best_cases[c(12, 1:3, 5:9)]
  
  
  gt(ungroup(best_cases)) %>% 
    tab_header(title = html("Best Estimate Per Case in 1d (&sigma; fixed)")) %>% 
    tab_stubhead(label = "Case") %>%
    cols_label(
      case= "Case",
      mu= html("&mu;<sub>fixed</sub>"),
      sigma= html("&sigma;<sub>fixed</sub>"),
      prior= html("Prior"),
      best_est = html("&mu;<sub>est</sub>"),
      kl="KL",
      map_est = html("&mu;<sub>MAP</sub>"),
      est_low= html("CI<sub>low</sub>"),
      est_hi=html("CI<sub>high</sub>"),
    ) -> base_table
  
  # add spanners
  base_table <- base_table %>% 
    tab_spanner(label = md("*Configuration*"),
                columns = c(case, mu, sigma, prior)) %>% 
    tab_spanner(label = md("*Post*"),
                columns = c(best_est)) %>% 
    tab_spanner(label = md("*Score*"),
                columns = c(kl))%>% 
    tab_spanner(label = md("*Credibility Interval*"),
                columns = c(map_est, est_low, est_hi))
  
  
  # add colors  
  base_table <- base_table  %>%
    data_color(
      columns = c(best_est),
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
      columns = c(map_est, est_low, est_hi),
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
      columns = c(prior),
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