

# no matrix
halfnormal<-function(k, theta, s, logscale = FALSE){
  # k <- c(kx, ky)
  # theta<-c(h0,sigma)
  kx <- k[1]
  ky <- k[2]
  h0 <- exp(theta[1])
  sigma <- exp(theta[2])
  h <- h0*exp(-( (kx - s[1])^2 + (ky - s[2])^2 ) / (2*sigma))
  if(logscale){h <- log(h)}
  return(as.numeric(h))
}

# # older version
# halfnormal<-function(k,theta,s){
#   h0<-exp(theta[1])
#   sigma<-exp(theta[2])
#   h<-h0*exp(-t(k-s)%*%(k-s)/(2*sigma))
#   return(h)
# }

# Hazard function for a trap at location k at time t.
# Theta are the parameters, s the individual's activity centre
# and z the location where it was last observed.


# removing matrix ops makes faster
hazard<-function(k, theta, t, z, s, m, logscale = FALSE){
  #theta<-c(h0,sigma,beta)
  #k<-c(kx,ky)
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


# # older version
# hazard<-function(k,theta,t,z,s){
#   h0<-exp(theta[1])
#   sigma<-exp(theta[2])
#   beta<--exp(theta[3])
#   I<-diag(2)
#   B<-exp(beta*t)*I
#   Lambda<-sigma*I
#   mu<-B%*%z+(I-B)%*%s
#   Sigma<-Lambda-B%*%Lambda%*%B
#   h<-h0*exp(-(1/2)*t(k-mu)%*%solve(Sigma)%*%(k-mu))
#   return(h)
# }

# Total hazard function for an individual with activity centre s
# last seen at location z. Theta are the parameters and 
# trap contains the locations of all the traps in the study. 
total_hazard<-function(trap,theta,t,z,s,m){
  h.<-sum(apply(trap,1,hazard,theta=theta,t=t,z=unlist(z),s=unlist(s),m=m))
  return(h.)
}

# Approximation of the survival function between times t1 and t2 
# for an individual with activity center s last seen at location z.
# The parameters are theta and trap contains the locations of 
# all the traps in the study. 
Survival<-function(trap,theta,t,i,j,z,s,m){
  Surv<-exp(-(t[i]-t[i-1])*total_hazard(trap,theta, t[i-1]-t[j]+((t[i]-t[i-1])/2),z,s,m))
  return(Surv)
}


# Simulates capture times and traps for each of the N individuals (although some will remain unobserved). 
# ac contains the N activity centres, trap the locations of all the traps in the study and theta the parameters.
# The study is run from time 0 to T. To simulate the capture times we discretize this time frame in r segments of
# equal length. 
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
  # Now we remove from the dataset the unobserved individuals
  df<-data
  #i<-1
  for (j in 1:N){
    if (identical(df[df$id==j,]$y,rep(0,length(time)-1))){
      df<-df[!(df$id==j),]
    }
    #else{
    #  df[df$id==j,]$id<-i
    #  i<-i+1
    #}
  }
  df<-df[!df$y==0,]
  rownames(df) <- NULL
  return(df)
}

re_id<-function(df,ac){
  #record the seen animals activity centers
  seen_animals<-unique(df$id)
  seen_ac<-ac[seen_animals,]
  
  #re-number the seen animals id's starting from 1 
  n<-length(unique(df$id)) # number of observed animals
  ids<-as.data.frame(table(df$id))
  ids$Var1<-c(1:n)
  ids <- ids[rep(ids$Var1,ids$Freq),1:(ncol(ids)-1)]
  df$id<-ids
  return(df)
}


# For simulating from OU process
sim<-function(theta, s,t,T){
  beta<--exp(theta[2])
  sigma<-exp(theta[1])^2
  I<-diag(2)
  start<-mvrnorm(n=1, mu=s,Sigma=sigma*I)
  trk<-data.frame(start[1],start[2])
  names(trk)<-c("x","y")
  for (i in 1:T){
    z<-trk[i,]
    mu<-as.double(exp(beta*t)*z+(1-exp(beta*t))*s)
    Cov<- sigma*I-(1/sigma)*exp(2*beta*t)*I
    mvrnorm(n=1, mu=mu,Sigma=Cov)
    bivariate_data <-mvrnorm(n=1, mu=mu,Sigma=Cov)
    trk[nrow(trk) + 1,] = c(bivariate_data[1],bivariate_data[2])
  }
  return(trk)
}

simulate_capthist<-function(theta, N, T, K){
  traps = make.grid(sqrt(K), sqrt(K), spacex = 1, detector = "multi")
  s_xy<-runif(2*N, min = -1, max = sqrt(K)+1) # simulate the location of the N activity centers
  ac<-matrix(s_xy,ncol=2,nrow=N) # place them in a matrix
  new_row_names <- 1:nrow(traps)
  rownames(traps) <- new_row_names
  df<-data.frame(matrix(nrow = 0, ncol = 3))
  
  for (i in 1:N){
    s<-ac[i,]
    trk<-sim(theta,s,1,T) # what happens when we change t here?? 
    for (j in 1:T){
      dist<-proxy::dist(as.matrix(traps), trk)[,j]
      seen<-which(dist<0.1)
      if (length(seen)!=0){
        df<-rbind(df,c(i,j,seen))
      }
    }
  }
  colnames(df)<-c("id","Time","y")
  df<-df%>% mutate(id=dense_rank(id))
  capt <- data.frame(session = rep("test",each=dim(df)[1]),
                     ID = df$id,
                     occasion = ceiling(df$Time),
                     trapID = as.character(df$y),
                     stringsAsFactors = FALSE)
  capthist<-make.capthist(capt,traps)
  return(capthist)
}

capthist_to_df<-function(capthist){
  df<-data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1: dim(capthist)[1]){
    var<-which(capthist[i,,]!=0,arr.ind=TRUE)
    for (j in 1: dim(var)[1]){
      df<-rbind(df,c(i,var[j,]))
    }
  }
  colnames(df)<-c("id","Time","y")
  df<-df[order(df$id, df$Time), ]   
  return(df)
}
