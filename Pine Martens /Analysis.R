## Analysis of the pine martens dataset 

library(lubridate)
library(dplyr)
library(ggplot2)
library(tibble)
library(secr)
library(sf)
Rcpp::sourceCpp("/Users/clara/Documents/RProjects/SECR_memory/LikelihoodC.cpp")
source("/Users/clara/Documents/RProjects/SECR_memory/Sim_Func.R")
source("/Users/clara/Documents/RProjects/SECR_memory/Fit_Func.R")


## load the original data 
data<-read.csv("/Users/clara/Documents/RProjects/Pine_Martens/Data analysis/pine_martens.csv")
traps<-read.csv("/Users/clara/Documents/RProjects/Pine_Martens/Data analysis/CameraLocations.csv")

## prepare the capture histories and traps data:
## - keep only the data from the 2017 field season of the study area number 3. 
## - remove the unmarked individuals
## - set the starting date as the 27th of February 2017 and create a variable "Time" that
##   gives the time of a capture as the days since the starting date. 
## - for the traps we only keep the coordinates. The ID of a trap is its row number.
## - for the capture histories keep only three variables representing an individual's 
##   unique ID, the time of the capture and the ID of the camera trap where the capture
##   occurred. 
data <- data %>%
  filter(!is.na(FieldSeason),
         FieldSeason == 2017,
         IndividualID != 0,
         StudyAreaID == 3) %>%
  mutate(
    date_obs = mdy(gsub("/", "-", DefaultStart)),
    start_date = ymd(20170227),
    day_obs = as.numeric(date_obs - start_date),
    time = hms(StartTime),
    hours = hour(time),
    minutes = minute(time),
    Time = round(day_obs + (hours * 60 + minutes) / 1440, digits = 3),
    id = match(IndividualID, unique(IndividualID)),
    y = LocationID
  ) %>%
  select(id, Time, y) %>%
  arrange(id, Time)

traps <- traps %>%
  filter(!is.na(FieldSeason)) %>%
  filter(FieldSeason == 2017) %>%
  filter(StudyAreaID == 3) %>%
  select(LocationID, UTM_E, UTM_N) %>%
  column_to_rownames("LocationID") %>%
  select(UTM_E, UTM_N) %>%
  rename(x = UTM_E, y = UTM_N)%>%
  mutate(
    x = x / 1000,
    y = y / 1000,
    name = rownames(.)
  ) %>%
  select(x, y)

unique_traps<-unique(rownames(traps))
data$y<-match(data$y,unique_traps)
rownames(traps)<-c(1:30)
traps<-subset(traps, select = c(x,y))


## create the mask over the region; different buffers were tested 
## until the estimated parameters remained stable
trap = make.poly(x=traps$x, y=traps$y)
mask = make.mask(trap,buffer=2,spacing=0.2,type="trapbuffer")
meshmat<-as.matrix(mask)
traps<-as.matrix(traps)

## produce Figure 1 of the paper
#plot(mask3,dots=FALSE,border=0,xlim=c(310,325),ylim=ylim,ppoly=FALSE,asp=1)
#plotMaskEdge(mask3,add=TRUE)
#points(trap3,pch=18)
#title("Study Area")
#legend_label <- "Camera Traps"
#legend("topright", legend = legend_label, pch = 18, col = "black", cex = 0.8, inset = c(0.05,0.12))


## set the length of the survey as T=12 days, the time 
## discretization with L=100 time intervals, the 
## number of observed individuals n 
T<-12
L<-100
n<-length(unique(data$id))

## make a version of the data set with the discretized time
ddf <- data.frame(t = as.numeric(), y = as.integer(), id = as.integer())
for(i in unique(data$id)){
  df <- discretize(data[data$id == i, ],T, L)
  df$id <- i
  ddf <- rbind(ddf, df)
}
ddfmat = as.matrix(ddf)
dfrows = as.numeric(table(ddf$id))

## inital parameters for MSCR (h0,sigma,beta)
theta<-c(1.5,-2,1.3) 


## fit the MSCR model and save the running time 
start_time <- Sys.time()
fit <- optim(theta, LikelihoodC, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T,hessian=TRUE) 
end_time <- Sys.time()
fit_time<-end_time - start_time

## extract the estimated parameters and confidence intervals
theta_est<-fit$par
param<-confint_param(fit,T,traps,mask,n)
N_est<-confint_pop(fit,T,traps,mask,n,distribution = "poisson", loginterval = TRUE, alpha = 0.05)



## initial parameters for SCR (h0,sigma)
theta<-c(1.5,-2) 

## fit the SCR model and save the running time
start_time <- Sys.time()
fit_nomem <- optim(theta, LikelihoodCnoMem, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T,hessian=TRUE) 
fit_time_nomem<-end_time <- Sys.time()
end_time - start_time

## extract the estimated parameters and confidence intervals
param_nm<-confint_param(fit_nomem,T,traps,mask,n)
N_nm<-confint_pop(fit_nomem,T,traps,mask,n,distribution = "poisson", loginterval = TRUE, alpha = 0.05)


## calculate the difference in AIC
theta_est_nomem<-fit_nomem$par
(2*fit_nomem$value-2*2)-(2*fit$value - 2*3)



