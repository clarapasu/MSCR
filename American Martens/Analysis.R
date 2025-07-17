## Analysis of the pine martens dataset 

library(lubridate)
library(ggplot2)
library(dplyr)
library(tibble)
library(secr)
library(sf)
Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Sim_Func.R")
source("Functions/Fit_Func.R")


## load the American marten data
## - data collection started on the 27th of February 2017
## - traps consists of UTM coordinates, the ID of a trap is its row number
## - in data we have the capture histories consisting of three variables 
##   representing an individual's unique ID, the time of the capture
##   and the ID of the camera trap where the capture occurred.
data<-read.csv("American Martens/marten_data.csv")
traps<-read.csv("American Martens/marten_traps.csv")


## create the mask over the region; different buffers were tested 
## until the estimated parameters remained stable
trap = make.poly(x=traps$x, y=traps$y)
trap <- trap[-31,]
mask = make.mask(trap,buffer=2,spacing=0.2,type="trapbuffer")
meshmat<-as.matrix(mask)
traps<-as.matrix(traps)



## to produce Figure 1 of the paper run the following code
 # mask_plot = make.mask(trap, buffer = 2, spacing = 0.01, type = "trapbuffer")
 # # Get plot limits from the mask
 # xlim <- range(mask_plot$x)
 # ylim <- range(mask_plot$y)
 # plot(mask_plot, dots = FALSE, border = 1, ppoly = FALSE, asp = 1, 
 #      xlab = "Longitude", ylab = "Latitude", main = "Study Area",
 #      xlim = xlim, ylim = ylim)
 # plotMaskEdge(mask_plot, add = TRUE, col = "black", lwd = 1.5) 
 # points(trap, pch = 4, col = "black", cex = 1.5)
 # axis(1, at = seq(floor(xlim[1]/5)*5, ceiling(xlim[2]/5)*5, by = 5),
 #      col = "black", col.axis = "black", cex.axis = 1.2)
 # 
 # axis(2, at = seq(floor(ylim[1]/5)*5, ceiling(ylim[2]/5)*5, by = 5),
 #      col = "black", col.axis = "black", cex.axis = 1.2)
 # mtext("Easting", side = 1, line = 3, cex = 1.1)
 # mtext("Northing", side = 2, line = 3, cex = 1.1)
 # legend("topright", legend = "Traps", pch = 4, col = "black", pt.cex = 1.5,
 #        cex = 1, bty = "n", text.col = "black")
 # title(main = "Study Area with Trap Locations", col.main = "black", font.main = 2, cex.main = 1.5)
 # 

## set the length of the survey as T=12 days, the time 
## discretization with L=100 time intervals, the 
## number of observed individuals n 
T<-11
L<-100
n<-length(unique(data$id))

## make a version of the data set with the discretized time
ddf <- data.frame(t = as.numeric(), y = as.integer(), id = as.integer())
for(i in unique(data$id)){
  df <- discretize(data[data$id == i, ],T, L)
  df$id <- i
  ddf <- rbind(ddf, df)
}
ddfmat = as.matrix(ddf)
dfrows = as.numeric(table(ddf$id))

## inital parameters for MSCR (h0,sigma,beta)
theta<-c(1.5, -2 ,1.5) 


## fit the MSCR model and save the running time 
start_time <- Sys.time()
fit <- optim(theta, LikelihoodC, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T,hessian=TRUE) 
end_time <- Sys.time()
fit_time<-end_time - start_time

## extract the estimated parameters and confidence intervals
theta_est<-fit$par
param<-confint_param(fit,T,traps,mask)
N_est<-confint_pop(fit,T,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)

## initial parameters for SCR (h0,sigma)
theta<-c(1.5,-2) 

## fit the SCR model and save the running time
start_time <- Sys.time()
fit_nomem <- optim(theta, LikelihoodCnoMem, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T,hessian=TRUE) 
fit_time_nomem<-end_time <- Sys.time()
time_nomem<- end_time - start_time
## extract the estimated parameters and confidence intervals
param_nm<-confint_param(fit_nomem,T,traps,mask)
N_nm<-confint_pop(fit_nomem,T,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)


## calculate the difference in AIC
theta_est_nomem<-fit_nomem$par
(2*fit_nomem$value-2*2)-(2*fit$value - 2*3)



