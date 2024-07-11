## AC PDF plots of a simulated individual 
## code to obtain Figure 4 a) and b) in the paper

library(lubridate)
library(dplyr)
library(ggplot2)
library(secr)
library(sf)
Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Sim_Func.R")
source("Functions/Fit_Func.R")

## load the traps obtained after clean up in "Analysis.R"
traps<-read.csv("marten_traps.csv",row.names=1)

## load a simulated dataset for which you want to plot the AC PDF
## and the MSCR and SCR models fitted to this dataset
df<-read.csv("df_03.csv",row.names=1)
ac<-read.csv("ac_03.csv",row.names=1)
fit <- get(load("model_3.Rdata")) 
fit_nomem <- get(load("model_nm_3.Rdata"))

## obtain the estimated parameters for each model 
theta_est<-fit$par
theta_est_nm<-fit_nomem$par

## make the mask
space=0.05
trap_mat<-as.matrix(traps)
trap_poly = make.poly(x=traps$x, y=traps$y)
mask = make.mask(trap_poly,buffer=2,spacing=space,type="trapbuffer")
meshmat<-as.matrix(mask)

## pick a simulated individual and calculate its AC PDF for MSCR
data<-df[df$id==7,]
ac_ind<-ac[7,]
trap_ind<-data.frame( x = trap_poly$x, y = trap_poly$y, seen = rep(0,31),seen_indic=rep(0,31))
for ( i in 1:length(table(data$y))){
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen<-table(data$y)[i]
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen_indic<-1
}
data<-discretize(data,12,10)
ac_posterior_ind7<-apply(meshmat,1,Likelihood_ind,theta=theta_est,trap=trap_poly,data=data,memory=1)
integral7<-as.double(exp(Likelihood_integrate(theta_est,trap_poly,data,mask,1)))
ac_density_ind7<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind7)/integral7)
colnames(ac_density_ind7)<-c("x","y","value")


## plot the AC PDF 

plot_ind<- ggplot(ac_density_ind7, aes(x = x, y = y, fill = value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(314,320))+
  ylim(c(4965,4970))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF",limits=c(0,100)) +
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill = ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",) +
  geom_point(data = data.frame(xaxis =  ac_ind[1], yaxis =ac_ind[2]), aes(x = xaxis, y = yaxis,color="activity centre") ,size=4,inherit.aes = FALSE)+
  scale_size_discrete(range = c(5, 10))




## repeat for SCR without memory
ac_posterior_ind_nm7<-apply(meshmat,1,Likelihood_ind,theta=theta_est_nm,trap=trap_poly,data=data,memory=0)
integral_nm7<-as.double(exp(Likelihood_integrate(theta_est_nm,trap_poly,data,mask,0)))
ac_density_ind_nm7<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind_nm7)/integral_nm7)
colnames(ac_density_ind_nm7)<-c("x","y","value")
ac_density_ind_nm7$scaled_value <- (ac_density_ind_nm7$value - min(ac_density_ind_nm7$value)) / (max(ac_density_ind_nm7$value) - min(ac_density_ind_nm7$value))*100


plot_ind_nm<- ggplot(ac_density_ind_nm7, aes(x = x, y = y, fill = scaled_value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(314,320))+
  ylim(c(4965,4970))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF",limits=c(0,100)) + 
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill = ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",) +
  geom_point(data = data.frame(xaxis =  ac_ind[1], yaxis =ac_ind[2]), aes(x = xaxis, y = yaxis,color="activity centre") ,size=4,inherit.aes = FALSE)+
  scale_size_discrete(range = c(5, 10))


