## functions to fit the MSCR and SCR model to datasets in the same format as can be obtained through simulations in "Sim_Func.R"
## data frames of captures consist of
## id: the individual's unique ID 
## Time: the time of the capture (assuming that the survey started at time 0)
## y: the index of the trap where the capture occurred 

## some variables come up in many functions:
## T: total length of survey 
## trap: vector of the traps' coordinates
## s: coordinates of an activity centre
## theta: model parameters, usually will contain (h0,sigma,beta) but for SCR it will only be (h0,sigma)
## r: number of bins in the time discretization 
## memory: indicator of whether we want the MSCR (memory=1) or SCR (memory=0) model



## discretization of the capture history of an individual from the dataset into r segments
## segments are cut again at the event times 
discretize<-function(df,T,r){
  t<-round(seq(0, T, T/r),digits=3)
  y<-rep(0,r+1)
  for (j in 1:length(df$Time)){
    time<-round(df$Time[j],digits=3)
    i<-findInterval(time, t)
    if (round(time,digits=2)==0.3){ 
      y[i+1]<-df$y[j]
    }
    else if (round(time,digits=2)==0.6){ 
      y[i+1]<-df$y[j]
    }
    else if (round(time,digits=2)==0.7){  
      y[i+1]<-df$y[j]
    }
    else if(round(time,digits=3)==round(t[i],digits=3)){
      y[i]<-df$y[j]
    }
    else{
      t<-append(t,time,after=i)
      y<-append(y,df$y[j],after=i)
    }
  }
  data<-data.frame(t,y)
  return(data)
}

## likelihood of an individual, conditional on its activity centre s
## its capture history "data" is in the format obtained after 
## discretization using the previous function
Likelihood_ind<-function(theta,trap,data,s,memory){
  m<-dim(data)[1]
  L<-1
  t<-data$t
  y<-data$y
  capt_time<-t[1]
  start_hazards<-sum(apply(trap,1,halfnormal,theta=theta,s=unlist(s)))
  first_capt<-which(data$y!=0)[1]
  capt_time<-t[first_capt]
  z<-trap[y[first_capt],]
  L<-exp(-capt_time*start_hazards)*halfnormal(z,theta,s)
  j<-first_capt
  for (i in (first_capt+1):m){
    L<-L*Survival(trap,theta,t,i,j,z,s,memory)
    if (y[i]!=0){
      k<-trap[y[i],]
      L<-L*hazard(k,theta,(t[i]-capt_time),z,s,memory)
      z<-k
      capt_time<-t[i]
      j<-i
    }
  }
  return(L)
}


## integration over space of the likelihood of an individual
## mesh is the grid over space that we use for the integration
Likelihood_integrate<-function(theta,trap,data,mesh,memory){
  a = attr(mesh,"a") # grid cell area
  D<-dim(mesh)[1] # number of grid cells in the mesh
  A = a*D # area of the mesh
  mesh<-matrix(c(mesh[,1],mesh[,2]),ncol=2,nrow=D)
  L<-0
  for (i in 1:D){
    L<-L+Likelihood_ind(theta,trap,data,mesh[i,],memory)
  }
  return(log(L*a/A))
}


## probability of being observed at least once in the survey
## conditional on the activity centre s
seen<-function(T,trap,s,theta){
  h<-sum(apply(trap,1,halfnormal,theta=theta,s=s))
  U<-1-exp(-T*h)
  return(U)
}

## integration of the activity centre out of the previous function
Seen_int<-function(T,trap,mesh,theta){
  a = attr(mesh,"a")
  D<-dim(mesh)[1]
  A = a*D
  mesh<-matrix(c(mesh[,1],mesh[,2]),ncol=2,nrow=D)
  S<-0
  for (i in 1:D){
    S<-S+seen(T,trap,as.double(mesh[i,]),theta)
  }
  return(log(S*a/A))
}


## likelihood of the full dataset df, where the activity centres have been integrated out
Likelihood<-function(theta,trap,df,mesh,T,r,memory){
  n<-length(unique(df$id))
  L<-0
  for (i in 1:n){ 
    data<-discretize(df[df$id==i,],T,r)
    L<-L+Likelihood_integrate(theta,trap,data,mesh,memory)
  }
  L<-L-n*Seen_int(T,trap,mesh,theta)
  return(-L)
}


## after optimizing the likelihood we obtain a model fit, this function 
## gives the confidence intervals of the model parameters 
confint_param<-function(fit,T,trap,mesh){
  theta_est<-fit$par
  Hessian<-fit$hessian
  cov<-solve(Hessian)
  se<-sqrt(diag(cov)) 
  upper<-fit$par+1.96*se
  lower<-fit$par-1.96*se
  interval<-data.frame(value=fit$par, upper=upper, lower=lower)
  return(interval)
}

# calculates assymmetric CIs because variance is calculated on log scale and backtransformed
add.cl <- function (df, alpha, loginterval, lowerbound = 0) 
{
  z <- abs(qnorm(1 - alpha/2))
  if (loginterval) {
    delta <- df$estimate - lowerbound
    df$lcl <- delta/exp(z * sqrt(log(1 + (df$SE.estimate/delta)^2))) + 
      lowerbound
    df$ucl <- delta * exp(z * sqrt(log(1 + (df$SE.estimate/delta)^2))) + 
      lowerbound
  }
  else {
    df$lcl <- pmax(lowerbound, df$estimate - z * df$SE.estimate)
    df$ucl <- df$estimate + z * df$SE.estimate
  }
  df
}


## estimate of the population size N and its confidence interval
## n is the number of observed individuals in the dataset 
# confint_pop<-function(theta_est,cov,T,trap,mesh,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05){
#   Seen_est<-exp(Seen_int(T,trap,mesh,theta_est)) 
#   Pop_size<-function(theta,T,trap,mesh,n){
#     Seen_est_t<-exp(Seen_int(T,trap,mesh,theta)) 
#     N_est<-n/Seen_est_t
#     return(N_est)
#   }
#   N_est<-Pop_size(theta_est,T,trap,mesh,n)
#   d <- numDeriv::grad(Pop_size,theta_est, T=T,trap=trap,mesh=mesh,n=n)
#   s2 <- switch (tolower(distribution),
#                 poisson  = n/Seen_est^2,
#                 binomial = n*((1-Seen_est)/(Seen_est^2)))
#   Var <- d%*% cov %*%d + s2
#   temp <- data.frame(row.names = c('N'), estimate = N_est,  SE.estimate = sqrt(Var))
#   temp <- add.cl(temp, alpha=alpha, loginterval)
#   return(temp)
# }

## estimate of the population size N and its confidence interval
## n is the number of observed individuals in the dataset 
confint_pop<-function(fit,T,trap,mesh,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05){
  theta_est<-fit$par
  Hessian<-fit$hessian
  cov<-solve(Hessian)
  Seen_est<-exp(Seen_int(T,trap,mesh,theta_est)) 
  Pop_size<-function(theta,T,trap,mesh,n){
    Seen_est_t<-exp(Seen_int(T,trap,mesh,theta)) 
    N_est<-n/Seen_est_t
    return(N_est)
  }
  N_est<-Pop_size(theta_est,T,trap,mesh,n)
  d <- numDeriv::grad(Pop_size,theta_est, T=T,trap=trap,mesh=mesh,n=n)
  s2 <- switch (tolower(distribution),
                poisson  = n/Seen_est^2,
                binomial = n*((1-Seen_est)/(Seen_est^2)))
  Var <- d%*% cov %*%d + s2
  temp <- data.frame(row.names = c('N'), estimate = N_est,  SE.estimate = sqrt(Var))
  temp <- add.cl(temp, alpha=alpha, loginterval)
  return(temp)
}

