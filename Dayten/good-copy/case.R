source(paste0(getwd(), "/good-copy/simpsons2d.R"))

# just sets up the environment per case
case <- function(n, case_num, parm1, parm2, grid_size) {
  print(paste(n, case_num, parm1, parm2))
  
  # this cases data from a pair of mu/sigma parameters
  data <-  rlogitnorm(n, parm1, parm2)
  assign('og_data', data, env = .GlobalEnv)
  assign(paste0('dat', case_num), data, env = .GlobalEnv)
  
  # plot & curve
  # data %>% hist(., probability=T, ylim=c(0,1.9))
  # curve(dlogitnorm(x, 0, 3.6), from=0, to=1, add=T, col='blue')
  
  # assuming we have logitnorm data, convert normal with logit transformation
  # to estimate mu and sigma
  mu0 <<- (logit(og_data) %>% mean())
  sigma0 <<- (sqrt(logit(og_data) %>% var()))
  
  #write h1, h2, S, xx, and yy to global.env
  simpsons2d(mu0, sigma0, grid_size)
}

# TEST
# for(case_num in 1:dim(parms)[1]){
#   print(paste("case", case_num))
#
#   # setup simpsons integration matrix per param pair
#   case(n, case_num, parms[case_num, 1], parms[case_num, 2]) }
