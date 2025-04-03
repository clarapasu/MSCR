# Simulating MSCR and continuous-time SCR data


#' @param theta List of model parameters (h0, sigma, beta). Ignore beta if using continuous-time SCR model. 
#' @param k List of coordinates of a camera trap
#' @param s List of coordinates of an activity center
#' @param m Memory indicator, 1 when using MSCR and 0 for continous-time SCR
#' @param N: Total population size
#' @param ac Matrix of dimension Nx2, containing the activity center coordinates for the N individuals.
#' @param trap Matrix of all of the camera trap coordinates. 




#' halfnormal distribution
halfnormal<-function(k, theta, s, logscale = FALSE){
  kx <- k[1]
  ky <- k[2]
  h0 <- exp(theta[1])
  sigma <- exp(theta[2])
  h <- h0*exp(-( (kx - s[1])^2 + (ky - s[2])^2 ) / (2*sigma))
  if(logscale){h <- log(h)}
  return(as.numeric(h))
}

## hazard function with an OU shape, t represents the time since the last capture 
hazard<-function(k, theta, t, z, s, m, logscale = FALSE){
  kx <- k[1]
  ky <- k[2]
  h0<-exp(theta[1])   
  sigma <-exp(theta[2])
  if (m==0){
    beta<--exp(100)}
  if (m==1){
    beta <- -exp(theta[3])}
  B<-exp(beta*t)  
  mu<-B*z+(1-B)*s
  if(logscale){
    h <- log(h0) + ((-(1/2) * ( (kx - mu[1])^2 + (ky - mu[2])^2 )  / (sigma-B^2*sigma)))
  } else {
    h <- h0 * exp(-(1/2) * ( (kx - mu[1])^2 + (ky - mu[2])^2 )  / (sigma-B^2*sigma))
    }
  return(h)
}

## cumulative hazard function
total_hazard<-function(trap,theta,t,z,s,m){
  h.<-sum(apply(trap,1,hazard,theta=theta,t=t,z=unlist(z),s=unlist(s),m=m))
  return(h.)
}

## Approximation of the survival function between times t[i] and t[i-1]
## as a constant, where the last capture occurred at time t[j]
Survival<-function(trap,theta,t,i,j,z,s,m){
  Surv<-exp(-(t[i]-t[i-1])*total_hazard(trap,theta, t[i-1]-t[j]+((t[i]-t[i-1])/2),z,s,m))
  return(Surv)
}


## simulation of N capture histories over time [0,T]
## only returns the capture histories of observed individuals
sim_data<-function(N,T,r,theta,trap,ac,m,mesh){
  data<-data.frame(matrix(ncol = 3, nrow = 0))
  names<-c("id","t","y")
  colnames(data)<-names
  time<-seq(0,T,T/r)
  K<-dim(trap)[1]
  for (n in 1:N){
    s<-ac[n,]
    y<-c()
    capt_time<-time[1] 
    start_hazards<-apply(mesh,1,halfnormal,theta=theta,s=s)
    z<-mesh[sample(1:nrow(mesh), size = 1, prob = start_hazards),]
    j<-1

    for (i in 2:length(time)){
      
          ti<-time[i]
          proba<-Survival(trap,theta,time,i,j,z,s,m)
          u<-runif(1)
          if (u<=proba){
            x<-0
          }
          else{
            traps_proba<-apply(trap,1,hazard,theta=theta,t=(ti-capt_time),z=z,s=s,m=m)
            p<-traps_proba/sum(traps_proba)
            x<-sample(length(traps_proba),1,prob=p)
            z<-trap[x,]
            capt_time<-ti
            j<-i
          }
          y<-c(y,x)

      }
    
    id<-rep(n,length(time)-1)
    Time<-time[-1]
    dat<-data.frame(id,Time,y)
    data<-rbind(data,dat)
  }
  df<-data
  for (j in 1:N){
    if (identical(df[df$id==j,]$y,rep(0,length(time)-1))){
      df<-df[!(df$id==j),]
    }
  }
  df<-df[!df$y==0,]
  rownames(df) <- NULL
  return(df)
}


## re-number the IDs of the previously simulated dataset starting from 1 
re_id<-function(df){
  n<-length(unique(df$id)) 
  ids<-as.data.frame(table(df$id))
  ids$Var1<-c(1:n)
  ids <- ids[rep(ids$Var1,ids$Freq),1:(ncol(ids)-1)]
  df$id<-ids
  return(df)
}
