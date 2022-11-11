library(dplyr)

x1 <- rnorm(300, -1, 1.25)
x2 <- rnorm(600, 5, 2.5)

c(x1, x2) %>% hist(breaks=35)
