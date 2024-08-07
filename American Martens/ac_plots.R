## Code to make the AC PDF plots (Figure 5 in the paper)


library(lubridate)
library(dplyr)
library(ggplot2)
library(secr)
library(sf)
Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/SECR_memory/Sim_Func.R")
source("Functions/Fit_Func.R")


# load the traps and data set
traps<-read.csv("marten_traps.csv",row.names=1)
df<-read.csv("marten_data.csv",row.names=1)

# load the models and results from the file "Analysis.R"
fit <- get(load("memory_model.Rdata")) 
fit_nomem <- get(load("nomemory_model.Rdata"))
theta_est<-fit$par
theta_est_nm<-fit_nomem$par

## create the mask, the variable "space" determines how coarse the mask is,
## can be adapted if code is too slow  
space=0.05
trap_mat<-as.matrix(traps)
trap_poly = make.poly(x=traps$x, y=traps$y)
mask = make.mask(trap_poly,buffer=2,spacing=space,type="trapbuffer")
meshmat<-as.matrix(mask)


## select an individual of the dataset and the traps where it was seen
data<-df[df$id==5,]
xmin<-314
xmax<-320
ymin<-4958
ymax<-4963
trap_ind<-data.frame( x = trap_poly$x, y = trap_poly$y, seen = rep(0,31),seen_indic=rep(0,31))
for ( i in 1:length(table(data$y))){
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen<-table(data$y)[i]
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen_indic<-1
}
data<-discretize(data,12,100)

## calculate its AC PDF for MSCR
ac_posterior_ind<-apply(meshmat,1,Likelihood_ind,theta=theta_est,trap=trap_poly,data=data,memory=1)
integral<-as.double(exp(Likelihood_integrate(theta_est,trap_poly,data,mask,1)))
ac_density_ind<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind)/integral)
colnames(ac_density_ind)<-c("x","y","value")

## obtain the corresponding AC PDF plot
fill_seen_zero <- "white"
plot_ind<- ggplot(ac_density_ind, aes(x = x, y = y, fill = value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(315,319))+
  ylim(c(4959,4962))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF",limits=c(0,200)) +
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill = ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",)+
  scale_size_discrete(range = c(5, 10))


## repeat with SCR, without memory
ac_posterior_ind_nm<-apply(meshmat,1,Likelihood_ind,theta=theta_est_nm,trap=trap_poly,data=data,memory=0)
integral<-as.double(exp(Likelihood_integrate(theta_est_nm,trap_poly,data,mask,0)))
ac_density_ind_nm<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind_nm)/integral)
colnames(ac_density_ind_nm)<-c("x","y","value")

plot_ind_nm <- ggplot(ac_density_ind_nm, aes(x = x, y = y, fill = value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(315,319))+
  ylim(c(4959,4962))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF",limits=c(0,200)) +  
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill = ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",)+
  scale_size_discrete(range = c(5, 10))



## repeat with another individual 
xmin<-314
xmax<-318
ymin<-4964
ymax<-4967

data<-df[df$id==7,]
trap_ind<-data.frame( x = trap_poly$x, y = trap_poly$y, seen = rep(0,31),seen_indic=rep(0,31))
for ( i in 1:length(table(data$y))){
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen<-table(data$y)[i]
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen_indic<-1
}

seen_traps_ind<-matrix(trap_mat[unique(data$y),],ncol=2)
times_seen_ind<-as.vector(table(data$y))

data<-discretize(data,12,100)
ac_posterior_ind<-apply(meshmat,1,Likelihood_ind,theta=theta_est,trap=trap_poly,data=data,memory=1)
integral<-as.double(exp(Likelihood_integrate(theta_est,trap_poly,data,mask,1)))
ac_density_ind<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind)/integral)
colnames(ac_density_ind)<-c("x","y","value")

plot_ind2<- ggplot(ac_density_ind, aes(x = x, y = y, fill = value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(xmin,xmax))+
  ylim(c(ymin,ymax))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF",limits=c(0,60)) +  
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill = ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",)+
  scale_size_discrete(range = c(5, 10))



ac_posterior_ind_nm<-apply(meshmat,1,Likelihood_ind,theta=theta_est_nm,trap=trap_poly,data=data,memory=0)
integral<-as.double(exp(Likelihood_integrate(theta_est_nm,trap_poly,data,mask,0)))
ac_density_ind_nm<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind_nm)/integral)
colnames(ac_density_ind_nm)<-c("x","y","value")
plot_ind2_nm <- ggplot(ac_density_ind_nm, aes(x = x, y = y, fill = value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(xmin,xmax))+
  ylim(c(ymin,ymax))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF",limits=c(0,60)) +  
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill = ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",)+
  scale_size_discrete(range = c(5, 10))



## and with a third individual
xmin<-314
xmax<-318
ymin<-4965
ymax<-4968

data<-df[df$id==8,]
trap_ind<-data.frame( x = trap_poly$x, y = trap_poly$y, seen = rep(0,31),seen_indic=rep(0,31))
for ( i in 1:length(table(data$y))){
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen<-table(data$y)[i]
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen_indic<-1
}

seen_traps_ind<-matrix(trap_mat[unique(data$y),],ncol=2)
times_seen_ind<-as.vector(table(data$y))

data<-discretize(data,12,100)
ac_posterior_ind<-apply(meshmat,1,Likelihood_ind,theta=theta_est,trap=trap_poly,data=data,memory=1)
integral<-as.double(exp(Likelihood_integrate(theta_est,trap_poly,data,mask,1)))
ac_density_ind<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind)/integral)
colnames(ac_density_ind)<-c("x","y","value")

plot_ind3<- ggplot(ac_density_ind, aes(x = x, y = y, fill = value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(xmin,xmax))+
  ylim(c(ymin,ymax))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF",limits=c(0,120)) +  
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill =ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",)+
  scale_size_discrete(range = c(5, 10))



ac_posterior_ind_nm<-apply(meshmat,1,Likelihood_ind,theta=theta_est_nm,trap=trap_poly,data=data,memory=0)
integral<-as.double(exp(Likelihood_integrate(theta_est_nm,trap_poly,data,mask,0)))
ac_density_ind_nm<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind_nm)/integral)
colnames(ac_density_ind_nm)<-c("x","y","value")

plot_ind3_nm <- ggplot(ac_density_ind_nm, aes(x = x, y = y, fill = value) ) +   
  geom_raster() +
  coord_fixed() +
  xlim(c(xmin,xmax))+
  ylim(c(ymin,ymax))+
  scale_fill_gradient(low = "white", high = "black", name = "AC PDF",limits=c(0,120)) +  
  scale_color_grey(name = "") + 
  ggtitle("") +
  geom_point(data = trap_ind, aes(x = x, y = y,size=as.factor(seen)), shape=23,fill = ifelse(trap_ind$seen == 0, fill_seen_zero, "grey"), color="black", inherit.aes = FALSE)+
  labs(fill = "Hazard value", size= "Number of captures",)+
  scale_size_discrete(range = c(5, 10))


