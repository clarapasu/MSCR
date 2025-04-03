## Simulation Study 


require(secr)
library(ggplot2)
library(numDeriv)
library(dplyr)
library(expm)
library(fields)
library(ggplot2)
library(proxy)
library(fdrtool)
library(ggpubr)
library(MASS)
library(tidyr)


Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Sim_Func.R")
source("Functions/Fit_Func.R")


set.seed(1)

## load the traps
traps<-read.csv("marten_traps.csv",col.names = c("x","y"))
trap = make.poly(x=traps$x, y=traps$y)
trap <- trap[-31,]
cams<- read.traps(data = traps,detector="count")


## set the length of survey T,
## the total number of individuals N, 
## the number of bins in the time discretization used to generate data r,
## the memory indicator m (m=1 means we have memory),
## the parameters h0 and sigma,
## the different values of the parameter beta that we will test,
## and the initial parameter values

T<-12
N<-20
r<-1000 
m<-1 
h0<- 0.85
sigma<- -1.16
beta<- c(-1.5,-1,-0.5,0,0.5,1,1.5)
K <- 30
theta_init<- c(1.00,  0.01 , 0.00)


## simulate the activity centers we pick a very small spacing
samplemesh = make.mask(trap,buffer=2,spacing=0.05,type="trapbuffer" )
samplemask<-as.matrix(samplemesh)
mask = make.mask(trap,buffer=2,spacing=0.5,type="trapbuffer")
meshmat<-as.matrix(mask)

## obtain the area of the region
dim(mask)
a = attr(mask,"a")
D<-dim(mask)[1]
A = a*D



## create a for loop to run 100 simulations per value of beta
n_sim <- 100

## initialize a dataframe for the results
results<-data.frame(matrix(ncol = 25, nrow = n_sim*length(beta)))
colnames(results)<-c( "h0", "sigma", "beta", "SE.h0", "SE.sigma", "SE.beta",
                     "N", "SE.N", "h0.nomem", "sigma.nomem", "SE.h0.nomem", "SE.sigma.nomem", "N.nomem", "SE.N.nomem",
                     "lambda0", "sigma", "SE.lambda0", "SE.sigma", "N.secr", "SE.secr", "n", "true.beta", 
                      "time.secr", "time.memory", "time.no.memory")

j<-1

# loop over the values of beta
for (k in 1:length(beta)){
  
  theta<-c(h0,sigma,beta[k])
  
  # run 100 simulations for each value of beta
  for (i in 1:n_sim){
    
    # sample an activity center for each simulated individual
    sampling<-sample(nrow(samplemask),N,replace=TRUE)
    ac<-as.matrix(samplemask[sampling,])
    
    # simulate a dataset
    df_sim<-sim_data(N,T,r,theta,as.matrix(cams),ac,m,samplemask)
    df_sim<-re_id(df_sim)
    df_sim$trap_x<-trap[df_sim$y,1]
    df_sim$trap_y<-trap[df_sim$y,2]
    df_sim<-arrange(df_sim,id,Time)
    
    
    # define the number of bins r2 used for the time discretization
    # and obtain the discretized version of the dataset
    
    r2<-20
    n<-length(unique(df_sim$id))
    
    n_rows <- sum(table(df_sim$id))  
    ddf_sim <- data.frame(
      t = numeric(n_rows),
      y = integer(n_rows),
      id = integer(n_rows)
    )
    
    row_idx <- 1
    
    for(i in unique(df_sim$id)){
      data <- discretize(df_sim[df_sim$id == i, ],T, r2)
      data$id <- i
      n_new <- nrow(data)
      ddf_sim[row_idx:(row_idx + n_new - 1), ] <- data  
      row_idx <- row_idx + n_new
    }
    ddfmat_sim = as.matrix(ddf_sim)
    dfrows_sim = as.numeric(table(ddf_sim$id))
    
    
    # fit MSCR 
    start_time <- Sys.time()
    fit <- optim(theta_init, LikelihoodC, trap = as.matrix(cams),df = ddfmat_sim, dfrows = dfrows_sim, mesh = meshmat, endt = T, hessian=TRUE)
    end_time <- Sys.time()
    time_memory<- end_time - start_time
    N_est<-confint_pop(fit,T,trap,mask,n)


    results[j,1:3] <- fit$par
    results[j,4:6] <- sqrt(diag(solve(fit$hessian)))
    results[j,7:8] <- N_est[1:2]
    
    
    # filter the data to have at most one observation per hour and repeat the discretization
    df_sim <- df_sim %>%
      arrange(id, Time) %>% 
      filter(!(id == lag(id) & (Time - lag(Time) < 0.04166667)))
    
    n_rows <- sum(table(df_sim$id))
    ddf_sim <- data.frame(
      t = numeric(n_rows),
      y = integer(n_rows),
      id = integer(n_rows)
    )
    
    row_idx <- 1
    
    for(i in unique(df_sim$id)){
      data <- discretize(df_sim[df_sim$id == i, ],T, r2)
      data$id <- i
      n_new <- nrow(data)
      ddf_sim[row_idx:(row_idx + n_new - 1), ] <- data
      row_idx <- row_idx + n_new
      ddfmat_sim = as.matrix(ddf_sim)
      dfrows_sim = as.numeric(table(ddf_sim$id))
    }
    
    # fit CT SCR
    start_time <- Sys.time()
    fit_nomem <- optim(theta[1:2], LikelihoodCnoMem, trap = as.matrix(cams),df = ddfmat_sim, dfrows =dfrows_sim, mesh = meshmat, endt = T, hessian=TRUE)
    end_time <- Sys.time()
    time_nomemory<- end_time - start_time
    
    N_est_nomem<-confint_pop(fit_nomem,T,trap,mask,n)
    
    
    results[j,9:10] <- fit_nomem$par
    results[j,11:12] <- sqrt(diag(solve(fit_nomem$hessian)))
    results[j,13:14] <-N_est_nomem[1:2]
    

    capt <- data.frame(session = rep("test",each=dim(df_sim)[1]),
                       ID = df_sim$id,
                       occasion = ceiling(df_sim$Time),
                       trapID = as.character(df_sim$y),
                       stringsAsFactors = FALSE)
    
    
    sim_ch<-make.capthist(capt,cams)
    
    
    start_time <- Sys.time()
    afit = secr.fit(sim_ch,mask=mask,detectfn = 14,trace = FALSE)
    end_time <- Sys.time()
    time_scr<- end_time - start_time
    
    
    N_secr <- exp(coef(afit)[1,1]) * A
    SE_N_secr<- A* sqrt(diag(solve(afit$fit$hessian)))[1]* exp(coef(afit)[1,1])
    results[j,15:16] <- afit$fit$estimate[2:3]
    results[j,17:18] <- sqrt(diag(solve(afit$fit$hessian)))[2:3]
    results[j,19] <- N_secr
    results[j,20] <- SE_N_secr
    
    
    
    results[j,21]<-n
    results[j,22]<-beta[k]
    results[j,23:25]<-c(time_scr,time_memory,time_nomemory)
    
    
  
    
    j<-j+1
  }
  
}




# process the format of the results before plotting
plot_data <- results %>%
  dplyr::select(true.beta, N.secr, N, N.nomem) %>%
  pivot_longer(
    cols = c(N.secr, N, N.nomem),
    names_to = "model",
    values_to = "N"
  ) %>%
  mutate(
    model = case_when(
      model == "N.secr" ~ "DT SCR",
      model == "N" ~ "MSCR",
      model == "N.nomem" ~ "CT SCR"
    ),
    model = factor(model),
    true.beta = as.factor(true.beta)
  )



# Plot the results of the three models for the different values of beta (Figure 3 in the paper)
ggplot(plot_data, aes(x = factor(true.beta), y = N, fill = model)) +
  geom_boxplot() +
  geom_hline(yintercept = 20, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("white", "gray80", "gray50")) +
  labs(title = "Population Estimates by Model",
       x = expression(log(beta)), y = "Population Estimate (N)") +
  theme_minimal() +
  theme(legend.position = "bottom")









