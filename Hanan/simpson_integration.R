## numerical integration: 1d
## Simpson's rules (1/3) and (3/8), with uniform subdivision
## Thurs. Nov 17, 2022

##Simpson's 1/3 rule
simps_mtd <- function(f, a, b, n){
  # function f
  # [a,b] interval of integration
  # n = number of BIG subdivisions (will be multiplied by two to include middle ones)
  nn = 2*n+1 
  h=(b-a)/n #length of subdivision
  xx = seq(a,b,length=nn) 
  y = f(xx)
  #pass in vector of integrand values (y) and length of subdivision (h)
  n=length(y)
  h*( y[1] + y[nn] + 4*sum(y[(2:(nn-1))[2:(nn-1)%%2==0]]) + 2*sum(y[(2:(nn-1))[2:(nn-1)%%2==1]]) )/6
}

## Simpson's 3/8 rule
simps_three8 <- function(f, a, b, n){
  # function f
  # [a,b] interval of integration
  # n = number of BIG subdivisions (will be multiplied by three to include middle ones)
  nn = 3*n+1
  h = (b - a)/n
  xx <- seq(a, b, length=nn)
  y = f(xx)
  
  h*( y[1] + y[nn] + 3*sum(y[(2:(nn-1))[2:(nn-1)%%3!=1]])  + 2*sum(y[(2:(nn-1))[2:(nn-1)%%3==1]]) )/8
}

###### A test function:
omega=10
f <- function(x) cos(omega*x) #integrand

a=0; b=1 #lower and upper bounds of integral
n = 1000 #number of subdivisions


sin(omega)/omega #true value
simps_mtd(f,a,b,n)-sin(omega)/omega      #1/3 rule error
simps_three8(f,a,b,n)-sin(omega)/omega   #3/8 rule error

########## Consider the following~~~
########
###### something more relevant to our project:

# suppose our likelihood is X ~ Normal(2,3) and our prior information is improper (i.e. =1):
# so our posterior is also Normal(2,3). Lets integrate this:

mu=2; sd=3
simps_three8(function(x) dnorm(x, mu, sd), 0, 1, n)
## check against true cdf:
(pnorm(1, mu, sd)-pnorm(0, mu, sd))

# now suppose our posterior is slightly more janky:
posterior <- function(x){
  dnorm(x, mu, sd)*dcauchy(x,location=0,scale=1)
}

# we use Simpson's rule to integrate it:
simps_three8(function(x) dnorm(x, mu, sd)*dcauchy(x), 0, 1, n)
#compare it against the base-R function:
integrate(posterior,0,1,subdivisions = n)

##Next: Gaussian quadrature (1d?) & monte-carlo integration (dw it is not related to markov chains lol) for 2d+



