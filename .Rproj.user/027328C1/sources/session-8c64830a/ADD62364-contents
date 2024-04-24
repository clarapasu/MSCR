

# Let's assume that the dataset comes in continuous time like we simulated, i.e. for each individual
# we have a list of capture times and traps where the animal was seen. 
# The discretize functions creates a discetization of the time from 0 to T in r segments.
# Segments in which the individual is unseen gets the value 0, while segments where the individual is
# observed get cut in two at the event time. The segment until the event time gets the value of the trap where 
# it is seen, and so on. 
# This function does this for 1 individual. 



seen<-function(T,trap,s,theta){
  h<-sum(apply(trap,1,halfnormal,theta=theta,s=s))
  U<-1-exp(-T*h)
  return(U)
}


discretize<-function(df,T,r){
  t<-round(seq(0, T, T/r),digits=3)
  y<-rep(0,r+1)
  for (j in 1:length(df$Time)){
    time<-round(df$Time[j],digits=3)
    i<-findInterval(time, t)
    if (round(time,digits=2)==0.3){  # This is because of the 0.3 issue in computers, still weird.... 
      y[i+1]<-df$y[j]
    }
    else if (round(time,digits=2)==0.6){  # Also weird
      y[i+1]<-df$y[j]
    }
    else if (round(time,digits=2)==0.7){  # I don't understand this
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


Likelihood_ind<-function(theta,trap,data,s,memory){
  m<-dim(data)[1]
  L<-1
  t<-data$t
  y<-data$y
  seen<-0
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
  #seen_prob<-seen(t[m],trap,s,theta)
  return(L)
}



Likelihood_test<-function(theta,trap,df,ac,T,r,memory){
  n<-length(unique(df$id))
  L<-0 
  for (i in 1:n){ 
    data<-discretize(df[df$id==i,],T,r)
    seen_prob<-seen(T,trap,ac[i,],theta)
    L<-L+log(Likelihood_ind(theta,trap,data,ac[i,],memory))-log(seen_prob)
    # This is not working anymore, it did when I divided by seen_prob inside Likelihood_ind
    # so I am a bit confused
  }
  return(-L)
}

Likelihood_integrate<-function(theta,trap,data,mesh,memory){
  a = attr(mesh,"a") # grid cell area
  D<-dim(mesh)[1]
  A = a*D
  mesh<-matrix(c(mesh[,1],mesh[,2]),ncol=2,nrow=D)
  L<-0
  for (i in 1:D){
    L<-L+Likelihood_ind(theta,trap,data,mesh[i,],memory)
  }
  return(log(L*a/A))
}


Seen_int<-function(T,trap,mesh,theta){
  a = attr(mesh,"a") # grid cell area
  D<-dim(mesh)[1]
  A = a*D
  mesh<-matrix(c(mesh[,1],mesh[,2]),ncol=2,nrow=D)
  S<-0
  for (i in 1:D){
    S<-S+seen(T,trap,as.double(mesh[i,]),theta)
  }
  return(log(S*a/A))
}

# Same as Seen_int but outputs one prob per mesh point, equivalent to secr's pdot
pdot.mem <- function(T,trap,mesh,theta){
  a = attr(mesh,"a") # grid cell area
  D<-dim(mesh)[1]
  A = a*D
  mesh<-matrix(c(mesh[,1],mesh[,2]),ncol=2,nrow=D)
  S<-c()
  for (i in 1:D){
    S<-c(S,seen(T,trap,as.double(mesh[i,]),theta))
  }
  return(S)
}


Likelihood<-function(theta,trap,df,mesh,T,r,memory){
  n<-length(unique(df$id))
  L<-0
  for (i in 1:n){ 
    data<-discretize(df[df$id==i,],T,r)
    L<-L+Likelihood_integrate(theta,trap,data,mesh,memory)
  }
  L<-L-n*Seen_int(T,trap,mesh,theta)
  #L<-L-Seen_int(T,trap,mesh,theta)
  return(-L)
}




Full_likelihood<-function(par,n,trap,df,mesh,T,r,memory){
  N<-exp(par[1])+n
  theta<-par[2:4]
  p<-exp(Seen_int(T,trap,mesh,theta))
  L<-log(factorial(N))-log(factorial(N-n))+n*log(p)+(N-n)*log(1-p) # use gamma function instead of factorial 
  Lik<-0
  for (i in 1:n){ 
    data<-discretize(df[df$id==i,],T,r)
    Lik<-Lik+Likelihood_integrate(theta,trap,data,mesh,memory)
  }
  L<-L+Lik-n*log(p)
  return(-L)
}




Likelihood_onefunc<-function(theta,trap,df,mesh,T,memory){
  
  n<-length(unique(df$id))
  
  A = (max(mesh[,1]) - min(mesh[,1])) * (max(mesh[,2]) - min(mesh[,2]))
  a = A / nrow(mesh) # grid cell area
  D <- nrow(mesh)
  
  S<-0
  start_hazards <- c()
  for (i in 1:D){
    start_hazards[i] <- 0
    for(j in 1:nrow(trap)){
      start_hazards[i] <- start_hazards[i] + halfnormal(trap[j,],theta=theta,s=mesh[i,])
    }
    U <- 1 - exp(-T* start_hazards[i])
    S <- S + U
  }
  LS <- log(S*a/A)
  
  LL <- 0 
  # sums log-likelihood over individuals
  for (h in 1:n){
    data <- df[df$id == h, 1:2]
    m <- dim(data)[1]
    
    L <- 0 
    # cat("L is ", LL, "\n")
    # sums up likelihood over possible AC locations (integrating out latent ACs)
    for (j in 1:D){
      t<-data$t
      y<-data$y
      s <- mesh[j, ]
      first_capt<-which(data$y!=0)[1]
      capt_time<-t[first_capt]
      z<-trap[y[first_capt],]
      Lj<-exp(-capt_time*start_hazards[j])*halfnormal(z,theta,s)
      seen_ind <- first_capt
      
      for (i in (first_capt+1):m){
        total_hazard <- 0
        for(jj in 1:nrow(trap)){
          #total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, (t[i] - t[i-1])/2, z, s, memory)
          total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, t[i-1]-t[seen_ind]+((t[i]-t[i-1])/2), z, s, memory)
        }
        survival <- exp(-(t[i]-t[i-1]) * total_hazard)
        Lj <- Lj * survival
        if (y[i]!=0){
          k<-trap[y[i],]
          Lj <- Lj * hazard(k,theta,(t[i]-capt_time),z,s,memory)
          z<-k
          capt_time<-t[i]
          seen_ind <- i
        }
        
        
      } # close loop over time points
      L <- L + Lj
    } # close loop over mesh points
    LL <- LL + log(L * a/A)
  } # close loop over individuals
  
  LL <- LL - n * LS
  return(-LL)
}

Likelihood_onefunc_log<-function(theta,trap,df,mesh,T,memory){
  
  n<-length(unique(df$id))
  
  A = (max(mesh[,1]) - min(mesh[,1])) * (max(mesh[,2]) - min(mesh[,2]))
  a = A / nrow(mesh) # grid cell area
  D <- nrow(mesh)
  
  S<-0
  start_hazards <- c()
  for (i in 1:D){
    start_hazards[i] <- 0
    for(j in 1:nrow(trap)){
      start_hazards[i] <- start_hazards[i] + halfnormal(trap[j,],theta=theta,s=mesh[i,])
    }
    U <- 1 - exp(-T* start_hazards[i])
    S <- S + U
  }
  LS <- log(S*a/A)
  
  LL <- 0 
  # sums log-likelihood over individuals
  for (h in 1:n){
    data <- df[df$id == h, 1:2]
    m <- dim(data)[1]
    
    L <- 0 
    # cat("L is ", LL, "\n")
    # sums up likelihood over possible AC locations (integrating out latent ACs)
    for (j in 1:D){
      t<-data$t
      y<-data$y
      s <- mesh[j, ]
      first_capt<-which(data$y!=0)[1]
      capt_time<-t[first_capt]
      z<-trap[y[first_capt],]
      logLj<- (-capt_time*start_hazards[j]) + halfnormal(z,theta,s,logscale=TRUE)
      seen_ind <- first_capt
      
      for (i in (first_capt+1):m){
        total_hazard <- 0
        for(jj in 1:nrow(trap)){
          #total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, (t[i] - t[i-1])/2, z, s, memory)
          total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, t[i-1]-t[seen_ind]+((t[i]-t[i-1])/2), z, s, memory)
        }
        survival <- exp(-(t[i]-t[i-1]) * total_hazard)
        logLj <- logLj + log(survival)
        if (y[i]!=0){
          k<-trap[y[i],]
          logLj <- logLj + hazard(k,theta,(t[i]-capt_time),z,s,memory,logscale=TRUE)
          z<-k
          capt_time<-t[i]
          seen_ind <- i
        }
        
        
      } # close loop over time points
      L <- L + exp(logLj)
    } # close loop over mesh points
    LL <- LL + log(L * a/A)
  } # close loop over individuals
  
  LL <- LL - n * LS
  return(-LL)
}

confint_param<-function(fit,T,trap,mesh,n){
  theta_est<-fit$par
  Hessian<-fit$hessian
  cov<-solve(Hessian)
  se<-sqrt(diag(cov)) 
  upper<-fit$par+1.96*se
  lower<-fit$par-1.96*se
  interval<-data.frame(value=fit$par, upper=upper, lower=lower)
  return(interval)
  
}

# KEEPING VERSION
# confint_pop<-function(fit,T,trap,mesh,n){
#   theta_est<-fit$par
#   Hessian<-fit$hessian
#   cov<-solve(Hessian)
#   Seen_est<-exp(Seen_int(T,trap,mesh,theta_est)) # the estimated probability of being observed
#   # Possible error below: Pop_size fn doesn't depend on theta, so grad wrt theta always (0,0)
#   Pop_size<-function(theta,T,trap,mesh,n,Seen_est){
#     N_est<-n/Seen_est # Horvitz-Thompson estimator of the population size
#     return(N_est)
#   }
#   N_est<-Pop_size(theta_est,T,trap,mesh,n,Seen_est=Seen_est)
#   d<- grad(Pop_size,theta_est, T=T,trap=trap,mesh=mesh,n=n,Seen_est=Seen_est)
#   Var<- d%*% cov %*%d + n*(1-Seen_est)/Seen_est
#   upper<-N_est+1.96*sqrt(Var)
#   lower<-N_est-1.96*sqrt(Var)
#   interval<-data.frame(N_est,lower,upper)
#   return(interval)
# }

# fixed version 11/5/2023


confint_pop<-function(fit,T,trap,mesh,n,distribution = "poisson", loginterval = TRUE, alpha = 0.05){
  theta_est<-fit$par
  Hessian<-fit$hessian
  cov<-solve(Hessian)
  Seen_est<-exp(Seen_int(T,trap,mesh,theta_est)) # the estimated probability of being observed
  Pop_size<-function(theta,T,trap,mesh,n){
    Seen_est_t<-exp(Seen_int(T,trap,mesh,theta)) # the estimated probability of being observed
    N_est<-n/Seen_est_t # Horvitz-Thompson estimator of the population size
    return(N_est)
  }
  N_est<-Pop_size(theta_est,T,trap,mesh,n)
  d <- numDeriv::grad(Pop_size,theta_est, T=T,trap=trap,mesh=mesh,n=n)
  s2 <- switch (tolower(distribution),
                poisson  = n/Seen_est^2,
                binomial = n*((1-Seen_est)/(Seen_est^2)))
  #Var<- d%*% cov %*%d + n*(1-Seen_est)/Seen_est
  #upper<-N_est+1.96*sqrt(Var)
  #lower<-N_est-1.96*sqrt(Var)
  Var <- d%*% cov %*%d + s2
  temp <- data.frame(row.names = c('N'), estimate = N_est,  SE.estimate = sqrt(Var))
  temp <- add.cl(temp, alpha=alpha, loginterval)
  return(temp)
}

# calculates assymmetric CIs because variance is calculated on log scale and backtransformed
# this function is taken from secr:::add.cl
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

# # manual grad gives same answers
# grad2 <- function (i, est, eps, T,trap, mesh,n) {
#   temp     <- est[i]
#   if (temp != 0.0) delta <- eps * abs(temp)
#   else             delta <- eps
#   est[i]  <- temp - delta
#   fminus  <- Pop_size(est,T,trap,mesh,n)
#   est[i]  <- temp + delta
#   fplus   <- Pop_size(est,T,trap,mesh,n)
#   (fplus - fminus) / (2.0 * delta)
# }
# sapply(1:length(theta), grad2, est=theta, eps=0.001, T=T,trap=trap,mesh=mesh,n=n)
# grad(Pop_size,theta, T=T,trap=trap,mesh=mesh,n=n)

# pac_plot <- function(theta,trap,df,mesh,T,memory){
#   
#   n<-length(unique(df$id))
#   
#   # set up list to contain Pr(AC_i|omega_t) for each mesh point i at time step t (omega_t is capt hist up to t)
#   pac <- list()
#   
#   A = (max(mesh[,1]) - min(mesh[,1])) * (max(mesh[,2]) - min(mesh[,2]))
#   a = A / nrow(mesh) # grid cell area
#   D <- nrow(mesh)
#   
#   S<-0
#   start_hazards <- c()
#   for (i in 1:D){
#     start_hazards[i] <- 0
#     for(j in 1:nrow(trap)){
#       start_hazards[i] <- start_hazards[i] + halfnormal(trap[j,],theta=theta,s=mesh[i,])
#     }
#     U <- 1 - exp(-T* start_hazards[i])
#     S <- S + U
#   }
#   LS <- log(S*a/A)
#   
#   LL <- 0 
#   df_ids <- unique(df$id)
#   # sums log-likelihood over individuals
#   for (h in 1:n){
#     data <- df[df$id == df_ids[h], 1:2]
#     m <- dim(data)[1]
#     
#     # set up matrix to store hazard of detection from mesh point at each time step
#     pac[[h]] <- matrix(0, nrow = D, ncol = m)
#     
#     L <- 0 
#     # cat("L is ", LL, "\n")
#     # sums up likelihood over possible AC locations (integrating out latent ACs)
#     for (j in 1:D){
#       t<-data$t
#       y<-data$y
#       s <- mesh[j, ]
#       first_capt<-which(data$y!=0)[1]
#       capt_time<-t[first_capt]
#       z<-trap[y[first_capt],]
#       Lj<-exp(-capt_time*start_hazards[j])*halfnormal(z,theta,s)
#       seen_ind <- first_capt
#       
#       for(i in 1:first_capt){
#         pac[[h]][, i] <- start_hazards
#       }
#       
#       for (i in (first_capt+1):m){
#         total_hazard <- 0
#         for(jj in 1:nrow(trap)){
#           #total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, (t[i] - t[i-1])/2, z, s, memory)
#           total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, t[i-1]-t[seen_ind]+((t[i]-t[i-1])/2), z, s, memory)
#         }
#         survival <- exp(-(t[i]-t[i-1]) * total_hazard)
#         Lj <- Lj * survival
#         if (y[i]!=0){
#           k<-trap[y[i],]
#           Lj <- Lj * hazard(k,theta,(t[i]-capt_time),z,s,memory)
#           z<-k
#           capt_time<-t[i]
#           seen_ind <- i
#         }
#         
#         pac[[h]][j, i] <- survival
#         
#       }
#       L <- L + Lj
#     }
#     LL <- LL + log(L * a/A)
#   }
#   
#   LL <- LL - n * LS
#   return(list(LL = -LL, pac = pac))
# }



pac_plot<-function(theta,trap,df,mesh,T,memory){
  
  n<-length(unique(df$id))
  
  A = (max(mesh[,1]) - min(mesh[,1])) * (max(mesh[,2]) - min(mesh[,2]))
  a = A / nrow(mesh) # grid cell area
  D <- nrow(mesh)
  
  pcapthist_given_ac <- matrix(0, nrow = n, ncol = D)
  
  S<-0
  start_hazards <- c()
  for (i in 1:D){
    start_hazards[i] <- 0
    for(j in 1:nrow(trap)){
      start_hazards[i] <- start_hazards[i] + halfnormal(trap[j,],theta=theta,s=mesh[i,])
    }
    U <- 1 - exp(-T* start_hazards[i])
    S <- S + U
  }
  LS <- log(S*a/A)
  
  LL <- 0 
  # sums log-likelihood over individuals
  for (h in 1:n){
    data <- df[df$id == h, 1:2]
    m <- dim(data)[1]
    
    L <- 0 
    # cat("L is ", LL, "\n")
    # sums up likelihood over possible AC locations (integrating out latent ACs)
    for (j in 1:D){
      t<-data$t
      y<-data$y
      s <- mesh[j, ]
      first_capt<-which(data$y!=0)[1]
      capt_time<-t[first_capt]
      z<-trap[y[first_capt],]
      Lj<-exp(-capt_time*start_hazards[j])*halfnormal(z,theta,s)
      seen_ind <- first_capt
      
      for (i in (first_capt+1):m){
        total_hazard <- 0
        for(jj in 1:nrow(trap)){
          #total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, (t[i] - t[i-1])/2, z, s, memory)
          total_hazard <- total_hazard + hazard(trap[jj,], theta=theta, t[i-1]-t[seen_ind]+((t[i]-t[i-1])/2), z, s, memory)
        }
        survival <- exp(-(t[i]-t[i-1]) * total_hazard)
        Lj <- Lj * survival
        if (y[i]!=0){
          k<-trap[y[i],]
          Lj <- Lj * hazard(k,theta,(t[i]-capt_time),z,s,memory)
          z<-k
          capt_time<-t[i]
          seen_ind <- i
        }
        
        
      } # close loop over time points
      L <- L + Lj
      pcapthist_given_ac[h,j] <- Lj
    } # close loop over mesh points
    LL <- LL + log(L * a/A)
  } # close loop over individuals
  
  LL <- LL - n * LS
  return(list(LL = -LL, pcapthist_given_ac = pcapthist_given_ac))
}





# Old code for the Likelihood of an individual
# Likelihood_ind_old<-function(theta,trap,data,s,memory){
#   m<-dim(data)[1]
#   L<-1
#   t<-data$t
#   y<-data$y
#   seen<-0
#   capt_time<-t[1]
#   start_hazards<-sum(apply(trap,1,halfnormal,theta=theta,s=unlist(s)))
#   for (i in 2:m){
#     if (seen==0){
#       L<-L*exp(-(t[i]-t[i-1])*start_hazards)
#       if (y[i]!=0){
#         k<-trap[y[i],]
#         L<-L*halfnormal(k,theta,s)
#         z<-k
#         capt_time<-t[i]
#         j<-i
#         seen<-1
#       }
#     }
#     else if (seen==1){
#       L<-L*Survival(trap,theta,t,i,j,z,s,memory)
#       
#       if (y[i]!=0){
#         k<-trap[y[i],]
#         L<-L*hazard(k,theta,(t[i]-capt_time),z,s,memory)
#         z<-k
#         capt_time<-t[i]
#         j<-i
#       }
#     }
#   }
#   #seen_prob<-seen(t[m],trap,s,theta)
#   return(L)
# }

