## Simulations from the MSCR model


library(proxy)
library(fdrtool)
library(secr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(MASS)
library(numDeriv)
Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Sim_Func.R")
source("Functions/Fit_Func.R")
set.seed(2)

## load the traps
traps<-read.csv("marten_traps.csv",row.names=1)

sim_mask<-make.mask(traps,buffer=2,spacing=0.05,type="trapbuffer")
sim_mesh<-as.matrix(sim_mask)
trap<-as.matrix(traps)

## set the length of survey T,
## the total number of individuals N, 
## the number of bins in the time discretization used to generate data r,
## the memory indicator m,
## and the parameters theta=(h0,sigma,beta)
T<-12 
N<-20 
r<-1000 
m<-1 
h0<- 0.5 
sigma<- -1.5 
beta<- -1	
theta<-c(h0,sigma,beta)
K <- 30


## simulate a dataset
Random_rows_1<-sample(nrow(sim_mesh),N)
ac<-sim_mesh[Random_rows_1,]
df_sim<-sim_data(N,T,r,theta,trap,ac,m)
df_sim<-re_id(df_sim,ac)
df_sim$trap_x<-trap[df_sim$y,1]
df_sim$trap_y<-trap[df_sim$y,2]
df_sim<-arrange(df_sim,id,Time)


# discretize the dataset and put it in the format we need it to optimize
r2<-10
n<-length(unique(df_sim$id))
ddf_sim <- data.frame(t = as.numeric(), y = as.integer(), id = as.integer())
for(i in unique(df_sim$id)){
  data <- discretize(df_sim[df_sim$id == i, ],T, r2)
  data$id <- i
  ddf_sim <- rbind(ddf_sim, data)
}
ddfmat_sim = as.matrix(ddf_sim)
dfrows_sim = as.numeric(table(ddf_sim$id))


## intial parameters for the optimization
theta_init<-c(-1, 0.01, 0)

## fit of the MSCR model
fit <- optim(theta_init, LikelihoodC, trap = trap,df = ddfmat_sim, dfrows = dfrows_sim, mesh = meshmat, endt = T, hessian=TRUE)
confint_param(fit,T,trap,meshmat)
confint_pop(fit,T,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)

## fit of the SCR model 
fit_nomem <- optim(theta_init[1:2], LikelihoodCnoMem, trap = trap,df = ddfmat_sim, dfrows = dfrows_sim, mesh = meshmat, endt = T, hessian=TRUE)
confint_param(fit_nomem,T,trap,meshmat)
confint_pop(fit_nomem,T,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)



## create a for loop to run 500 simulations

results_bis<-data.frame(matrix(ncol = 24, nrow = 500))
execution_times <- numeric()
execution_times_nm <- numeric()

for (j in 1:500){
  k=j
  Random_rows_1<-sample(nrow(sim_mesh),N)
  ac<-sim_mesh[Random_rows_1,]
  
  df_sim<-sim_data(N,T,r,theta,trap,ac,m)
  df_sim<-re_id(df_sim,ac)
  df_sim$trap_x<-trap[df_sim$y,1]
  df_sim$trap_y<-trap[df_sim$y,2]
  df_sim<-arrange(df_sim,id,Time)
  r2<-10
  n<-length(unique(df_sim$id))
  ddf_sim <- data.frame(t = as.numeric(), y = as.integer(), id = as.integer())
  for(i in unique(df_sim$id)){
    data <- discretize(df_sim[df_sim$id == i, ],T, r2)
    data$id <- i
    ddf_sim <- rbind(ddf_sim, data)
  }
  ddfmat_sim = as.matrix(ddf_sim)
  dfrows_sim = as.numeric(table(ddf_sim$id))
  
  start_time <- Sys.time()
  fit <- optim(theta_init, LikelihoodC, trap = trap,df = ddfmat_sim, dfrows = dfrows_sim, mesh = meshmat, endt = T, hessian=TRUE)
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  execution_times <- c(execution_times, elapsed_time)
  
  param<- confint_param(fit,T,trap,meshmat)
  N_est<-confint_pop(fit,T,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)
  
  start_time <- Sys.time()
  fit_nomem <- optim(theta[1:2], LikelihoodCnoMem, trap = trap,df = ddfmat_sim, dfrows = dfrows_sim, mesh = meshmat, endt = T, hessian=TRUE)
  end_time <- Sys.time()
  elapsed_time <-as.numeric(difftime(end_time, start_time, units = "secs"))
  execution_times_nm <- c(execution_times_nm, elapsed_time)
  
  param_nomem<-confint_param(fit_nomem,T,trap,meshmat)
  N_est_nomem<-confint_pop(fit_nomem,T,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)
  
  #data_name <- sprintf("df_%02d.csv", k) 
  #write.csv(df_sim,file=data_name)
  #fit_name <- paste("model_", k, ".RData", sep = "")
  #save(fit,file=fit_name)
  #fit_nomem_name <- paste("model_nm_", k, ".RData", sep = "")
 # save(fit_nomem,file=fit_nomem_name)
  trap_count<-df_sim %>% 
    group_by(id) %>%
    summarise(y_count = n_distinct(y), na.rm = TRUE)

  results_bis[j,]<-c(param$value,param$lower,param$upper,N_est[1],N_est[3:4],param_nomem$value,param_nomem$lower,param_nomem$upper,N_est_nomem[1],N_est_nomem[3:4],n,mean(table(df_sim$id)),mean(trap_count$y_count))
}
