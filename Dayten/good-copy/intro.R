# source(paste0(getwd(),'/Dayten/init.R'))
set.seed(90210)
custPalette <- c("#002754", "#005493",  "#C63527", "#F5AA1C")

parms <- matrix(
  c(0, 3.16, 0, 1.78, 0, .32,
    1, 3.16, 1, 1.78, 1, 0.32,
    2, 3.16, 2, 1.78, 2, 0.32),
  nrow=9, ncol=2, byrow = T
)
colnames(parms) <- c("mu", "sigma")

n <- 100000

if(1){ # DONT RUN THIS IF YOU DONT NEED TO
  dat  <- matrix(nrow = 0, ncol= 2)
  
  for (ii in 1:dim(parms)[1]) {
    new_dat <-
      cbind(rlogitnorm(n, parms[ii, 1], parms[ii, 2]), paste0("Case", ii))
    dat <- rbind(dat, new_dat)
  }
  
  colnames(dat) <- c("data", "group")
  dat <- dat %>% as.data.frame() %>% mutate_at("data", as.numeric)
  
  
  # visualize per "Case"
  plots_list <- list()
  for (ii in 1:dim(parms)[1]) {
    # retrieve this iter's data
    tmp <- dat %>% filter(group == paste0('Case',ii))
    
    # Freedman Diaconis bin-widt
   
    bins <- c(rep(c(2,0.175,0.2),2), 2, 2, 0.2)*sqrt(n)
    # tmp %>% ggplot(aes(x=dat, y=xat)) +  geom_point(aes(y=dat))
    #plot
    p <- tmp %>%
      ggplot(aes(x = data)) +
      geom_histogram(
        aes(y = ..density..),
        bins = bins[ii],
        fill = custPalette[2],
        alpha = 0.8
      ) +
      stat_function(
        fun = dlogitnorm,
        geom = "line",
        color=custPalette[4],
        args = list(mu = parms[ii, 1], sigma = parms[ii, 2])
      ) + 
      ggtitle(bquote(mu ~"=" ~ .(parms[ii, 1]) ~" "~ sigma ~ "=" ~ .(parms[ii, 2]))) + 
      theme(plot.title = element_text(hjust = 0.5))
      # geom_density(color = custPalette[4]) 
      # p
    plots_list[[ii]] <- p + theme_minimal() 
  }
}
library(gridExtra)

ggarrange(do.call("grid.arrange", c(plots_list, ncol=3)), nrow=1)
# ggsave(
#   "cases.png",
#   last_plot(),
#   device = "png",
#   width = 5,
#   height = 4,
#   units = "in",
#   dpi = 300,
#   bg="white"
# )
# 
