## Code to make the AC PDF plots

library(lubridate)
library(dplyr)
library(ggplot2)
library(secr)
library(sf)
library(patchwork)
library(scales)
library(grid)    

Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Sim_Func.R")
source("Functions/Fit_Func.R")

# load the traps and data set
df<-read.csv("American Martens/marten_data.csv")
traps<-read.csv("American Martens/marten_traps.csv")

# load the models and results from the file "Analysis.R"
fit <- get(load("American Martens/memory_model.Rdata")) 
fit_nomem <- get(load("American Martens/nomemory_model.Rdata"))
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
# in the paper we present the AC PDF plots for individual 2, 3 and 7
data<-df[df$id==2,]
trap_ind<-data.frame( x = trap_poly$x, y = trap_poly$y, seen = rep(0,31),seen_indic=rep(0,31))
for ( i in 1:length(table(data$y))){
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen<-table(data$y)[i]
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen_indic<-1
}
data<-discretize(data,11,100)

## calculate the AC PDF from the MSCR model
ac_posterior_ind<-apply(meshmat,1,Likelihood_ind,theta=theta_est,trap=trap_poly,data=data,memory=1)
integral<-as.double(exp(Likelihood_integrate(theta_est,trap_poly,data,mask,1)))
ac_density_ind<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind)/integral)
colnames(ac_density_ind)<-c("x","y","value")


## repeat with SCR, without memory
ac_posterior_ind_nm<-apply(meshmat,1,Likelihood_ind,theta=theta_est_nm,trap=trap_poly,data=data,memory=0)
integral<-as.double(exp(Likelihood_integrate(theta_est_nm,trap_poly,data,mask,0)))
ac_density_ind_nm<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind_nm)/integral)
colnames(ac_density_ind_nm)<-c("x","y","value")



# start plotting 

# define a shared scale for the legend
common_fill <- scale_fill_gradient(
  low     = "white",
  high    = "black",
  name    = "AC PDF",
  limits  = c(0, 1000),
  oob     = squish,
  trans   = "sqrt",
  breaks  = c(0, 250, 500, 750, 1000),
  labels  = c("0", "250", "500", "750", "1000"),
  guide   = guide_colourbar(
    barwidth  = unit(8, "cm"),   
    barheight = unit(0.5, "cm"),  
    ticks     = TRUE,
    frame.colour = "black"
  )
)

# produce the plot with and without memory

p1 <- ggplot(ac_density_ind, aes(x, y, fill = value)) +
  geom_raster() +
  coord_fixed() +
  common_fill +
  xlim(324, 329) + ylim(5007, 5012) +
  geom_point(
    data        = trap_ind,
    aes(x, y, size = as.factor(seen)),
    shape       = 23,
    inherit.aes = FALSE,
    color       = "black",
    fill        = ifelse(trap_ind$seen == 0, "white", "grey")
  ) +
  scale_size_discrete(
    range = c(5, 10),
    guide = guide_legend(
      override.aes = list(
        shape = 23,
        fill  = c("white", rep("grey", length(levels(as.factor(trap_ind$seen))) - 1))
      )
    )
  ) +
  labs( title = "Surface obtained from the MSCR model",
        x     = "Easting",
        y     = "Northing",
        size = "Number of captures") +
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "grey90", colour = NA)
  )

p2 <- ggplot(ac_density_ind_nm, aes(x, y, fill = value)) +
  geom_raster() +
  coord_fixed() +
  common_fill +
  xlim(324, 329) + ylim(5007, 5012) +
  geom_point(
    data        = trap_ind,
    aes(x, y, size = as.factor(seen)),
    shape       = 23,
    inherit.aes = FALSE,
    color       = "black",
    fill        = ifelse(trap_ind$seen == 0, "white", "grey")
  ) +
  scale_size_discrete(
    range = c(5, 10),
    guide = guide_legend(
      override.aes = list(
        shape = 23,
        fill  = c("white", rep("grey", length(levels(as.factor(trap_ind$seen))) - 1))
      )
    )
  ) +
  labs(title = "Surface obtained from the CT SCR model",
       x     = "Easting",
       y     = "Northing",
       size = "Number of captures") +
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "grey90", colour = NA)
  )


# combine them with a shared legend

spacer <- plot_spacer()


(p1 + spacer + p2) +
  plot_layout(
    ncol   = 3,
    widths = c(1, 0.1, 1),  
    guides = "collect"
  ) &
  theme(legend.position = "bottom")





# Repeat for the second individual:

data<-df[df$id==3,]
trap_ind<-data.frame( x = trap_poly$x, y = trap_poly$y, seen = rep(0,31),seen_indic=rep(0,31))
for ( i in 1:length(table(data$y))){
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen<-table(data$y)[i]
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen_indic<-1
}
data<-discretize(data,11,100)


ac_posterior_ind<-apply(meshmat,1,Likelihood_ind,theta=theta_est,trap=trap_poly,data=data,memory=1)
integral<-as.double(exp(Likelihood_integrate(theta_est,trap_poly,data,mask,1)))
ac_density_ind<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind)/integral)
colnames(ac_density_ind)<-c("x","y","value")


ac_posterior_ind_nm<-apply(meshmat,1,Likelihood_ind,theta=theta_est_nm,trap=trap_poly,data=data,memory=0)
integral<-as.double(exp(Likelihood_integrate(theta_est_nm,trap_poly,data,mask,0)))
ac_density_ind_nm<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind_nm)/integral)
colnames(ac_density_ind_nm)<-c("x","y","value")



common_fill <- scale_fill_gradient(
  low       = "white",
  high      = "black",
  name      = "AC PDF",
  limits    = c(0, 200),
  oob       = squish,
  trans     = "sqrt",
  breaks    = c(0, 50, 100, 150, 200),
  labels    = c("0", "50", "100", "150", "200"),
  na.value  = "grey90",
  guide     = guide_colourbar(
    barwidth     = unit(8, "cm"),
    barheight    = unit(0.5, "cm"),
    ticks        = TRUE,
    frame.colour = "black",
    order        = 2
  )
)


n_seen_levels <- length(levels(as.factor(trap_ind$seen)))
fill_ovr       <- c("white", rep("grey", n_seen_levels - 1))
size_guide     <- guide_legend(
  override.aes = list(
    shape = 23,
    fill  = fill_ovr
  ),
  order = 1
)


p1 <- ggplot(ac_density_ind, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  coord_fixed() +
  common_fill +
  xlim(324.5, 329.5) +
  ylim(5010.5, 5015.5) +
  geom_point(
    data        = trap_ind,
    aes(x = x, y = y, size = as.factor(seen)),
    shape       = 23,
    color       = "black",
    fill        = ifelse(trap_ind$seen == 0, "white", "grey"),
    inherit.aes = FALSE
  ) +
  scale_size_discrete(range = c(5, 10), guide = size_guide) +
  labs(
    title = "Surface obtained from the MSCR model",
    x     = "Easting",
    y     = "Northing",
    size  = "Number of captures"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "grey90", colour = NA)
  )

p2 <- ggplot(ac_density_ind_nm, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  coord_fixed() +
  common_fill +
  xlim(324.5, 329.5) +
  ylim(5010.5, 5015.5) +
  geom_point(
    data        = trap_ind,
    aes(x = x, y = y, size = as.factor(seen)),
    shape       = 23,
    color       = "black",
    fill        = ifelse(trap_ind$seen == 0, "white", "grey"),
    inherit.aes = FALSE
  ) +
  scale_size_discrete(range = c(5, 10), guide = size_guide) +
  labs(
    title = "Surface obtained from the CT SCR model",
    x     = "Easting",
    y     = "Northing",
    size  = "Number of captures"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "grey90", colour = NA)
  )


spacer <- plot_spacer()

(p1 + spacer + p2) +
  plot_layout(
    ncol   = 3,
    widths = c(1, 0.1, 1),   
    guides = "collect"
  ) &
  theme(legend.position = "bottom")



#Now for a third individual

data<-df[df$id==7,]
trap_ind<-data.frame( x = trap_poly$x, y = trap_poly$y, seen = rep(0,31),seen_indic=rep(0,31))
for ( i in 1:length(table(data$y))){
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen<-table(data$y)[i]
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen_indic<-1
}
data<-discretize(data,11,100)


ac_posterior_ind<-apply(meshmat,1,Likelihood_ind,theta=theta_est,trap=trap_poly,data=data,memory=1)
integral<-as.double(exp(Likelihood_integrate(theta_est,trap_poly,data,mask,1)))
ac_density_ind<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind)/integral)
colnames(ac_density_ind)<-c("x","y","value")



ac_posterior_ind_nm<-apply(meshmat,1,Likelihood_ind,theta=theta_est_nm,trap=trap_poly,data=data,memory=0)
integral<-as.double(exp(Likelihood_integrate(theta_est_nm,trap_poly,data,mask,0)))
ac_density_ind_nm<-data.frame(meshmat[,1],meshmat[,2],unlist(ac_posterior_ind_nm)/integral)
colnames(ac_density_ind_nm)<-c("x","y","value")



common_fill <- scale_fill_gradient(
  low       = "white",
  high      = "black",
  name      = "AC PDF",
  limits   = c(0, 200),
  oob      = squish,
  trans    = "sqrt",
  breaks   = c(0, 50, 100, 150, 200),
  labels   = c("0", "50", "100", "150", "200"),
  na.value  = "grey90",
  guide     = guide_colourbar(
    barwidth    = unit(8, "cm"),
    barheight   = unit(0.5, "cm"),
    ticks       = TRUE,
    frame.colour= "black",
    order=2
  )
)


n_seen_levels <- length(levels(as.factor(trap_ind$seen)))
fill_ovr       <- c("white", rep("grey", n_seen_levels - 1))


size_guide <- guide_legend(
  override.aes = list(
    shape = 23,
    fill  = c("white", rep("grey", n_seen_levels - 1))
  ),
  order = 1                 
)


p1 <- ggplot(ac_density_ind, aes(x, y, fill = value)) +
  geom_raster() +
  coord_fixed() +
  common_fill +
  xlim(325, 330) +
  ylim(5004, 5009) +
  geom_point(
    data        = trap_ind,
    aes(x, y, size = as.factor(seen)),
    shape       = 23,
    color       = "black",
    fill        = ifelse(trap_ind$seen == 0, "white", "grey"),
    inherit.aes = FALSE
  ) +
  scale_size_discrete(range = c(5,10), guide = size_guide) +
  labs(
    title = "Surface obtained from the MSCR model",
    x     = "Easting",
    y     = "Northing",
    size  = "Number of captures"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "grey90", colour = NA)
  )

p2 <- ggplot(ac_density_ind_nm, aes(x, y, fill = value)) +
  geom_raster() +
  coord_fixed() +
  common_fill +
  xlim(325, 330) +
  ylim(5004, 5009) +
  geom_point(
    data        = trap_ind,
    aes(x, y, size = as.factor(seen)),
    shape       = 23,
    color       = "black",
    fill        = ifelse(trap_ind$seen == 0, "white", "grey"),
    inherit.aes = FALSE
  ) +
  scale_size_discrete(range = c(5,10), guide = size_guide) +
  labs(
    title = "Surface obtained from the CT SCR model",
    x     = "Easting",
    y     = "Northing",
    size  = "Number of captures"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "grey90", colour = NA)
  )


spacer <- plot_spacer()

(p1 + spacer + p2) +
  plot_layout(
    ncol   = 3,
    widths = c(1, 0.1, 1),  
    guides = "collect"
  ) &
  theme(legend.position = "bottom")
