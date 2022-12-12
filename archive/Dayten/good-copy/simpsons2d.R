# sets up all the simpsons parameters to run our posterior estimation
simpsons2d <- function(mu0, sigma0, grid_size){
  # simpsons params
  nx = ny = grid_size
  n1 = 2 * nx + 1 # length of mu sequence
  n2 = 2 * ny + 1 # length of sigma sequence
  
  ax = mu0 - 0.5# range for mu +/- one
  bx = mu0 + 0.5
  ay = max(sigma0 - 0.55, 0.01) # range for sigma cant be <= 0, breaks estimation
  by = sigma0 + 0.75 # and 2 past current estimate
  
  h1 <<- (bx - ax) / (n1 - 1)  #length of subdivisions
  h2 <<- (by - ay) / (n2 - 1)
  
  #create the Simpson matrix:
  s1 = c(1, rep(2, n1 - 2) ^ (1:(n1 - 2) %% 2 + 1) , 1)
  s2 = c(1, rep(2, n2 - 2) ^ (1:(n2 - 2) %% 2 + 1) , 1)
  S <<- outer(s1, s2)
  
  # set out the sequences (used elsewhere)
  x_seq <<- seq(ax,bx,length=n1)
  y_seq <<- seq(ay,by,length=n2)
  
  # create the variable matrices
  assign("xx", matrix(x_seq, nrow=n1, ncol=n2), envir = .GlobalEnv)
  assign("yy", matrix(y_seq, nrow=n1, ncol=n2, byrow=T), envir = .GlobalEnv)
}

