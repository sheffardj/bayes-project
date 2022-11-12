if(!require(rprojroot)){install.packages("rprojroot"); library(rprojroot)}
root <- find_rstudio_root_file()
setwd(paste0(root,"/Dayten"))
cat(paste("PROJECT DIR:", paste0("'", getwd(), "'"), sep='\n'))
source("seed.R")
invisible(source("libraries.R"))
print("READY TO RUN CODE IN /DAYTEN")
