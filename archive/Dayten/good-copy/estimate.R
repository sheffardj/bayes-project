estimate <- function(case_num, ii, pri, prj){
  # set up prior
  LH.MAT <- lh(xx, yy, pri, prj)
  scale <- h1 * h2 * sum(S * LH.MAT) / 9 #compute the integral
  lh.est <- LH.MAT/scale # return posterior
  
  #convert matrix to paired vectors
  rownames(lh.est) <- x_seq
  colnames(lh.est) <- y_seq
  M = lh.est %>% as.data.frame()
  M["rowid"] = row.names(M)
  arr0 <- gather(M, colid, value, -rowid) %>%  arrange(rowid, colid)
  
  arr0 <- apply(arr0, 2, as.numeric) %>% as.data.frame()
  arr0 <- arr0 %>% arrange(rowid, colid)
  colnames(arr0) <- c("mu", "sigma", "prob")

  
  # plot layout
  # layout(mat = matrix(c(1, 3, 2, 3),
  #                     nrow = 2,
  #                     ncol = 2),
  #        heights = c(1, 1),
  #        widths = c(1, 1))

  # These functions come from the 'rethinking' package
  # contour_xyz(arr0$mu, arr0$sigma, arr0$prob, main=paste("Likelihood"))
  # image_xyz(arr0$mu, arr0$sigma, arr0$prob)
  # og_data %>% hist(., probability=T, ylim=c(0,1.9))
  # curve(dlogitnorm(x, 0, 3.16), from=0, to=1, add=T, col='blue')
  
  mu.post <- arr0[which.max(arr0[,3]),1]
  sig.post <- arr0[which.max(arr0[,3]),2]
  # curve(dlogitnorm(x, mu.post, sig.post), from=0, to=1, add=T, col='green')
  
  # zz %>% hist(., probability=T, ylim=c(0,4))
  # curve(dlogitnorm(x, 0, 3.16), from=0, to=1, add=T, col='blue', n=300)
  # curve(dlogitnorm(x, mu.post, sig.post), from=0, to=1, add=T, col='red', n=300)
  print(paste('est:', mu.post, sig.post))
  
  # assign(paste0("arr", ii), arr, envir = .GlobalEnv)
  arr <<- cbind(arr, arr0)
  assign("estimates", rbind(estimates, c(
    case = case_num,
    run_id = ii,
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