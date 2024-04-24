## AC PDF plots for simulations from an OU process 
## code to obtain Figure 4 c) and d) in the paper

library(lubridate)
library(dplyr)
library(ggplot2)
library(secr)
library(sf)
Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Sim_Func.R")
source("Functions/Fit_Func.R")

## load the traps obtained after clean up in "Analysis.R"
traps<-read.csv("pine_marten_traps.csv",row.names=1)


## load a simulated dataset for which you want to plot the AC PDF
## and the MSCR and SCR models fitted to this dataset
df<-read.csv("df_01.csv",row.names=1)
ac<-read.csv("ac_01.csv",row.names=1)
fit <- get(load("model_1.Rdata"))
fit_nomem <- get(load("model_nm_1.Rdata"))

## estimated parameters for each model 
theta_est<-fit$par
theta_est_nm<-fit_nomem$par

## create the mask 
space=0.05
trap_mat<-as.matrix(traps)
trap_poly = make.poly(x=traps$x, y=traps$y)
mask = make.mask(trap_poly,buffer=2,spacing=space,type="trapbuffer")
meshmat<-as.matrix(mask)


## pick a simulated individual and calculate its AC PDF for MSCR
data<-df[df$id==3,]
trap_ind<-data.frame( x = trap_poly$x, y = trap_poly$y, seen = rep(0,31),seen_indic=rep(0,31))
for ( i in 1:length(table(data$y))){
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen<-table(data$y)[i]
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen_indic<-1
}

data<-discretize(data,12,10)
ac_posterior_ind<-apply(meshmat,1,Likelihood_ind,theta=theta_est,trap=trap_poly,data=data,memory=1)
integral<-as.double(exp(Likelihood_integrate(theta_est,trap_poly,data,mask,1)))
ac_density_ind<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind)/integral)
colnames(ac_density_ind)<-c("x","y","value")

## plot the AC PDF
fill_seen_zero <- "white"
plot_ind<- ggplot(ac_density_ind, aes(x = x, y = y, fill = value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(312,318))+
  ylim(c(4960,4965))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF") +
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill = ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",) +
  geom_point(data = data.frame(xaxis =  314.0644, yaxis =4962.585), aes(x = xaxis, y = yaxis,color="activity centre") ,size=4,inherit.aes = FALSE)+
  scale_size_discrete(range = c(5, 10))

## repeat with SCR
ac_posterior_ind_nm<-apply(meshmat,1,Likelihood_ind,theta=theta_est_nm,trap=trap_poly,data=data,memory=0)
integral_nm<-as.double(exp(Likelihood_integrate(theta_est_nm,trap_poly,data,mask,0)))
ac_density_ind_nm<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind_nm)/integral_nm)
colnames(ac_density_ind_nm)<-c("x","y","value")

fill_seen_zero <- "white"
plot_ind<- ggplot(ac_density_ind_nm, aes(x = x, y = y, fill = value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(312,318))+
  ylim(c(4960,4965))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF") +
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill = ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",) +
  geom_point(data = data.frame(xaxis =  314.0644, yaxis =4962.585), aes(x = xaxis, y = yaxis,color="activity centre") ,size=4,inherit.aes = FALSE)+
  scale_size_discrete(range = c(5, 10))


