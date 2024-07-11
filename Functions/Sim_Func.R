## functions to simulate a dataset from the MSCR model
## in this file we encounter
## theta: model parameters, where theta=(h0,sigma,beta)
## k: coordinates of a trap 
## s: coordinates of an activity centre
## m: memory indicator, 1 if we want the model with memory and 0 without
## z: coordinates of a trap, in practice the one where the individual was last seen at
## N: total population size
## ac: vector of N activity centre coordinates
## trap: vector of all of the trap coordinates



## halfnormal distribution
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
sim_data<-function(N,T,r,theta,trap,ac,m){
  data<-data.frame(matrix(ncol = 3, nrow = 0))
  names<-c("id","t","y")
  colnames(data)<-names
  time<-seq(0,T,T/r)
  K<-dim(trap)[1]
  for (n in 1:N){
    s<-ac[n,]
    y<-c()
    capt_time<-time[1]     
    obs<-0
    for (i in 2:length(time)){
      
      if (obs==0){
        start_hazards<-apply(trap,1,halfnormal,theta=theta,s=s)
        start_surv<-exp(-(time[i]-time[i-1])*sum(start_hazards))
        u<-runif(1)
        if (u<=start_surv){
          x<-0
        }
        else{
          traps_proba<-start_hazards
          p<-traps_proba/sum(traps_proba)
          x<-sample(length(traps_proba),1,prob=p)
          z<-trap[x,]
          capt_time<-time[i]
          obs<-1
          j<-i
        }
        y<-c(y,x)
      }
      else{
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
re_id<-function(df,ac){
  #seen_animals<-unique(df$id)
  #seen_ac<-ac[seen_animals,]
  n<-length(unique(df$id)) 
  ids<-as.data.frame(table(df$id))
  ids$Var1<-c(1:n)
  ids <- ids[rep(ids$Var1,ids$Freq),1:(ncol(ids)-1)]
  df$id<-ids
  return(df)
}

## simulation of a trajectory from an OU process
## theta represents different parameters than previously in this case 
sim<-function(theta, s,t,T){
  beta<- -exp(theta[2])
  sigma<- exp(theta[1])^2
  I<-diag(2)
  start<-mvrnorm(n=1, mu=s,Sigma=sigma*I)
  trk<-data.frame(start[1],start[2])
  names(trk)<-c("x","y")
  for (i in 1:T){
    z<-trk[i,]
    mu<-as.double(exp(beta*t)*z+(1-exp(beta*t))*s)
    Cov<- sigma*I-(sigma)*exp(2*beta*t)*I
    mvrnorm(n=1, mu=mu,Sigma=Cov)
    bivariate_data <-mvrnorm(n=1, mu=mu,Sigma=Cov)
    trk[nrow(trk) + 1,] = c(bivariate_data[1],bivariate_data[2])
  }
  return(trk)
}


## function to simulate N capture histories using the previous function to simulate trajectories,
## a location closer than 50m from a trap is a capture
simulate_dataOU<-function(T,theta,sim_mesh,N){
  capthist<-matrix(nrow = 0, ncol = 3)
  colnames(capthist)<-c("id","Time","y")
  Random_rows_1<-sample(nrow(sim_mesh),N)
  ac<-sim_mesh[Random_rows_1,]
  print(ac) # remove this line 
  for (i in 1:N){
    s<-ac[i,]
    trk<-sim(theta,s,0.001,T)
    for (j in 1:T){
      dist<-proxy::dist(as.matrix(traps), trk)[,j]
      obs<-which(dist<0.05) #change the value here if we want to change the distance to a trap that is a capture
      if (length(obs)!=0 ){
        capthist<-rbind(capthist,c(i,j,obs))
      }
    }
  }
  df<-as.data.frame(capthist)
  df<-df%>% mutate(id=dense_rank(id))
  df$Time<-df$Time
    return(df)
}

