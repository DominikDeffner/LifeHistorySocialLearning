# Code for Fig. 3, birth rates and effective vital rates, you need simulation output with individual-level data, times of
#environmental changes and population sizes


#Calculate suff from output matrices

for (i in 1:nrow(seq)) {
  for (z in 1:10) {
    
    result[[i]][[z]][["TimeSinceChange"]] <- c()        #Time since last switch in environment
    result[[i]][[z]][["Babies"]] <- c()                 #Number of juveniels
    result[[i]][[z]][["Adults"]] <- c()                 #Number of adults   
    result[[i]][[z]][["RelNoJuveniles"]] <- c()         #Proportion of juveniles J_i
    result[[i]][[z]][["AdultAdapt"]] <- c()             #Adaptation level among adults q_A,i
    result[[i]][[z]][["JuvenileAdapt"]] <- c()          #Adaptation level among juveniles q_J,i
    result[[i]][[z]][["JuvenileXi"]] <- c()             #Mean xi among juveniles
    
    for (x in 1:7000) {
      
      if (all(result[[i]][[z]][["EnvChange"]][1:x] == 0)) {
        result[[i]][[z]][["TimeSinceChange"]][x] <- NA
      } else {result[[i]][[z]][["TimeSinceChange"]][x] <- x- max(which(result[[i]][[z]][["EnvChange"]][1:x]==1))}
      
      result[[i]][[z]][["Babies"]][[x]]        <- length(which(result[[i]][[z]][["Age"]][x,]==1))
      
      result[[i]][[z]][["Adults"]][[x]]        <- result[[i]][[z]][["Populationsize"]][[x]] - result[[i]][[z]][["Babies"]][[x]]
      
      result[[i]][[z]][["RelNoJuveniles"]][x]  <- result[[i]][[z]][["Babies"]][[x]]/result[[i]][[z]][["Populationsize"]][[x]]
      
      result[[i]][[z]][["AdultAdapt"]][x]    <- mean(result[[i]][[z]][["Adapt"]][x,][which(result[[i]][[z]][["Age"]][x,]!=1)], na.rm = TRUE)
      
      result[[i]][[z]][["JuvenileAdapt"]][x] <- mean(result[[i]][[z]][["Adapt"]][x,][which(result[[i]][[z]][["Age"]][x,]==1)], na.rm = TRUE)
      
      result[[i]][[z]][["JuvenileXi"]][x]    <- mean(result[[i]][[z]][["Xi"]][x,][which(result[[i]][[z]][["Age"]][x,]==1)], na.rm = TRUE)
      
    }
    result[[i]][[z]][["MeanAdapt"]] <- apply(result[[i]][[z]][["Adapt"]] , 1, mean, na.rm=TRUE) 
  }
}


#Use only last 5000 time steps

for (i in 1:nrow(seq)) {
  for (z in 1:10) {
    result[[i]][[z]][["Adapt"]]           <- result[[i]][[z]][["Adapt"]][2001:7000,]
    result[[i]][[z]][["Age"]]             <- result[[i]][[z]][["Age"]][2001:7000,]
    result[[i]][[z]][["EnvChange"]]       <- result[[i]][[z]][["EnvChange"]][2001:7000]
    result[[i]][[z]][["Populationsize"]]  <- result[[i]][[z]][["Populationsize"]][2001:7000]
    result[[i]][[z]][["TimeSinceChange"]] <- result[[i]][[z]][["TimeSinceChange"]][2001:7000]
    result[[i]][[z]][["Babies"]]          <- result[[i]][[z]][["Babies"]][2001:7000]
    result[[i]][[z]][["MeanAdapt"]]       <- result[[i]][[z]][["MeanAdapt"]][2001:7000]
    result[[i]][[z]][["AdultAdapt"]]      <- result[[i]][[z]][["AdultAdapt"]][2001:7000]
    result[[i]][[z]][["JuvenileAdapt"]]   <- result[[i]][[z]][["JuvenileAdapt"]][2001:7000]
    result[[i]][[z]][["JuvenileXi"]]      <- result[[i]][[z]][["JuvenileXi"]][2001:7000]
  }
}



#Calculate adaptation levels, juveniles and population size conditional on time since last switch
for (i in 1: nrow(seq)) {
  for (z in 1:10) {
    
    result[[i]][[z]][["AdaptPerTime"]] <- c()
    result[[i]][[z]][["BabiesPerTime"]] <- c()
    result[[i]][[z]][["NPerTime"]] <- c()
    
    for (x in unique(result[[i]][[z]][["TimeSinceChange"]])) {
      result[[i]][[z]][["AdaptPerTime"]][x]  <- mean(result[[i]][[z]][["MeanAdapt"]][which(result[[i]][[z]][["TimeSinceChange"]]==x)])
      result[[i]][[z]][["BabiesPerTime"]][x] <- mean(result[[i]][[z]][["Babies"]][which(result[[i]][[z]][["TimeSinceChange"]]==x)])
      result[[i]][[z]][["NPerTime"]][x]      <- mean(result[[i]][[z]][["Populationsize"]][which(result[[i]][[z]][["TimeSinceChange"]]==x)])
    }
  }
}


#calculate effective vital rates
for (i in 1:nrow(seq)) {
  s<- function(L_hat){(L_hat-1)/L_hat}
  b<- function(s,delta, N_hat){(1-s)/(s*exp(-delta*N_hat))}
  gamma <- function(s,b, delta){delta*(log(1/(s*(1+b)))/log((1-s)/(b*s)))}
  
  #Define vital rate parameters
  
  s <- s(seq$L_hat[i])                   
  b <- b(s=s,delta=1/1000,N_hat=500)  
  
  if (seq$regulation[i] == "Survival"){
    gamma<-gamma(s=s, b=b, delta=1/1000)  
    delta<-0                 
  }
  
  if (seq$regulation[i] == "Fertility"){
    delta<-1/1000                 
    gamma<-0                 
  }                
  
  for (z in 1:10) {
    result[[i]][[z]][["EffSurvival"]] <- c()
    result[[i]][[z]][["EffFertility"]] <- c()
    
    for (x in unique(result[[i]][[z]][["TimeSinceChange"]])) {
      
      J <-  mean(result[[i]][[z]][["RelNoJuveniles"]][which(result[[i]][[z]][["TimeSinceChange"]]==x)])    # Proportion Juveniles
      
      qA <- mean(result[[i]][[z]][["AdultAdapt"]][which(result[[i]][[z]][["TimeSinceChange"]]==x)])        # Proportion of adaptive behavior among adults
      
      qJ <- mean(result[[i]][[z]][["JuvenileAdapt"]][which(result[[i]][[z]][["TimeSinceChange"]]==x)])     # Proportion of adaptive behavior among juveniles
      
      Dens_surv <- s*exp(-gamma*result[[i]][[z]][["NPerTime"]][x])                                         # Density-dependent survival
      
      Dens_fert <- b*exp(-delta*result[[i]][[z]][["NPerTime"]][x])                                         # Density-dependent fertility
      
      Xi_J <-  mean(result[[i]][[z]][["JuvenileXi"]][which(result[[i]][[z]][["TimeSinceChange"]]==x)])     # Mean Xi among juveniles
      
      result[[i]][[z]][["EffFertility"]][x] <- J * 0 + (1-J)*(qA*1.1 + (1-qA))*Dens_fert  #Proportion adapted
      
      result[[i]][[z]][["EffSurvival"]][x] <- (1-J)*(qA*1.1 + (1-qA))*Dens_surv +
        J * ((1-Xi_J) * (qJ*1.1 + (1-qJ)) * Dens_surv +
               Xi_J  * (qJ*1.1 + (1-qJ)) * 0.95 * Dens_surv)
    }
  }
}




# Reduce results to necessary variables

for (i in 1: nrow(seq)) {
  for (z in 1:10) {
    result[[i]][[z]] <- result[[i]][[z]][c("EffFertility","EffSurvival","AdaptPerTime","BabiesPerTime","NPerTime")]
  }  
}


#color stuff
library(scales)
require(RColorBrewer)#load package
x <- seq(from=0, to=1, by=0.2) # fake data
col.pal <- brewer.pal(length(x), "Dark2") #create a pallette which you loop over for corresponding values


#Limit to first 1/u timesteps

for (i in 1: nrow(seq)) {
  for (x in 1:10) {
    result[[i]][[x]][["AdaptPerTime"]]  <- result[[i]][[x]][["AdaptPerTime"]][1:(1/seq$u[i])]
    result[[i]][[x]][["BabiesPerTime"]] <- result[[i]][[x]][["BabiesPerTime"]][1:(1/seq$u[i])]
    result[[i]][[x]][["NPerTime"]]      <- result[[i]][[x]][["NPerTime"]][1:(1/seq$u[i])]
    result[[i]][[x]][["EffFertility"]]  <- result[[i]][[x]][["EffFertility"]][1:(1/seq$u[i])]
    result[[i]][[x]][["EffSurvival"]]   <- result[[i]][[x]][["EffSurvival"]][1:(1/seq$u[i])]
  }
}

# Calculate mean and upper/lower limits over 10 independent simulations

for (i in 1: nrow(seq)) {
  result[[i]][["MeanN"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["MeanN"]][z] <- mean(c(result[[i]][[1]][["NPerTime"]][z],
                                        result[[i]][[2]][["NPerTime"]][z],
                                        result[[i]][[3]][["NPerTime"]][z],
                                        result[[i]][[4]][["NPerTime"]][z],
                                        result[[i]][[5]][["NPerTime"]][z],
                                        result[[i]][[6]][["NPerTime"]][z],
                                        result[[i]][[7]][["NPerTime"]][z],
                                        result[[i]][[8]][["NPerTime"]][z],
                                        result[[i]][[9]][["NPerTime"]][z],
                                        result[[i]][[10]][["NPerTime"]][z]))}

for (i in 1: nrow(seq)) {
  result[[i]][["UpperN"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["UpperN"]][z] <- max(c(result[[i]][[1]][["NPerTime"]][z],
                                        result[[i]][[2]][["NPerTime"]][z],
                                        result[[i]][[3]][["NPerTime"]][z],
                                        result[[i]][[4]][["NPerTime"]][z],
                                        result[[i]][[5]][["NPerTime"]][z],
                                        result[[i]][[6]][["NPerTime"]][z],
                                        result[[i]][[7]][["NPerTime"]][z],
                                        result[[i]][[8]][["NPerTime"]][z],
                                        result[[i]][[9]][["NPerTime"]][z],
                                        result[[i]][[10]][["NPerTime"]][z]))}

for (i in 1: nrow(seq)) {
  result[[i]][["LowerN"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["LowerN"]][z] <- min(c(result[[i]][[1]][["NPerTime"]][z],
                                         result[[i]][[2]][["NPerTime"]][z],
                                         result[[i]][[3]][["NPerTime"]][z],
                                         result[[i]][[4]][["NPerTime"]][z],
                                         result[[i]][[5]][["NPerTime"]][z],
                                         result[[i]][[6]][["NPerTime"]][z],
                                         result[[i]][[7]][["NPerTime"]][z],
                                         result[[i]][[8]][["NPerTime"]][z],
                                         result[[i]][[9]][["NPerTime"]][z],
                                         result[[i]][[10]][["NPerTime"]][z]))}


for (i in 1: nrow(seq)) {
  result[[i]][["MeanFert"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["MeanFert"]][z] <- mean(c(result[[i]][[1]][["EffFertility"]][z],
                                                    result[[i]][[2]][["EffFertility"]][z],
                                                    result[[i]][[3]][["EffFertility"]][z],
                                                    result[[i]][[4]][["EffFertility"]][z],
                                                    result[[i]][[5]][["EffFertility"]][z],
                                                    result[[i]][[6]][["EffFertility"]][z],
                                                    result[[i]][[7]][["EffFertility"]][z],
                                                    result[[i]][[8]][["EffFertility"]][z],
                                                    result[[i]][[9]][["EffFertility"]][z],
                                                    result[[i]][[10]][["EffFertility"]][z]))}

for (i in 1: nrow(seq)) {
  result[[i]][["UpperFert"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["UpperFert"]][z] <- max(c(result[[i]][[1]][["EffFertility"]][z],
                                           result[[i]][[2]][["EffFertility"]][z],
                                           result[[i]][[3]][["EffFertility"]][z],
                                           result[[i]][[4]][["EffFertility"]][z],
                                           result[[i]][[5]][["EffFertility"]][z],
                                           result[[i]][[6]][["EffFertility"]][z],
                                           result[[i]][[7]][["EffFertility"]][z],
                                           result[[i]][[8]][["EffFertility"]][z],
                                           result[[i]][[9]][["EffFertility"]][z],
                                           result[[i]][[10]][["EffFertility"]][z]))}

for (i in 1: nrow(seq)) {
  result[[i]][["LowerFert"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["LowerFert"]][z] <- min(c(result[[i]][[1]][["EffFertility"]][z],
                                           result[[i]][[2]][["EffFertility"]][z],
                                           result[[i]][[3]][["EffFertility"]][z],
                                           result[[i]][[4]][["EffFertility"]][z],
                                           result[[i]][[5]][["EffFertility"]][z],
                                           result[[i]][[6]][["EffFertility"]][z],
                                           result[[i]][[7]][["EffFertility"]][z],
                                           result[[i]][[8]][["EffFertility"]][z],
                                           result[[i]][[9]][["EffFertility"]][z],
                                           result[[i]][[10]][["EffFertility"]][z]))}


for (i in 1: nrow(seq)) {
  result[[i]][["MeanSurv"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["MeanSurv"]][z] <- mean(c(result[[i]][[1]][["EffSurvival"]][z],
                                           result[[i]][[2]][["EffSurvival"]][z],
                                           result[[i]][[3]][["EffSurvival"]][z],
                                           result[[i]][[4]][["EffSurvival"]][z],
                                           result[[i]][[5]][["EffSurvival"]][z],
                                           result[[i]][[6]][["EffSurvival"]][z],
                                           result[[i]][[7]][["EffSurvival"]][z],
                                           result[[i]][[8]][["EffSurvival"]][z],
                                           result[[i]][[9]][["EffSurvival"]][z],
                                           result[[i]][[10]][["EffSurvival"]][z]))}

for (i in 1: nrow(seq)) {
  result[[i]][["UpperSurv"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["UpperSurv"]][z] <- max(c(result[[i]][[1]][["EffSurvival"]][z],
                                           result[[i]][[2]][["EffSurvival"]][z],
                                           result[[i]][[3]][["EffSurvival"]][z],
                                           result[[i]][[4]][["EffSurvival"]][z],
                                           result[[i]][[5]][["EffSurvival"]][z],
                                           result[[i]][[6]][["EffSurvival"]][z],
                                           result[[i]][[7]][["EffSurvival"]][z],
                                           result[[i]][[8]][["EffSurvival"]][z],
                                           result[[i]][[9]][["EffSurvival"]][z],
                                           result[[i]][[10]][["EffSurvival"]][z]))}

for (i in 1: nrow(seq)) {
  result[[i]][["LowerSurv"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["LowerSurv"]][z] <- min(c(result[[i]][[1]][["EffSurvival"]][z],
                                           result[[i]][[2]][["EffSurvival"]][z],
                                           result[[i]][[3]][["EffSurvival"]][z],
                                           result[[i]][[4]][["EffSurvival"]][z],
                                           result[[i]][[5]][["EffSurvival"]][z],
                                           result[[i]][[6]][["EffSurvival"]][z],
                                           result[[i]][[7]][["EffSurvival"]][z],
                                           result[[i]][[8]][["EffSurvival"]][z],
                                           result[[i]][[9]][["EffSurvival"]][z],
                                           result[[i]][[10]][["EffSurvival"]][z]))}



for (i in 1: nrow(seq)) {
  result[[i]][["MeanAdaptPerTime"]] <- c()
for (z in 1: (1/seq$u[i]))
result[[i]][["MeanAdaptPerTime"]][z] <- mean(c(result[[i]][[1]][["AdaptPerTime"]][z],
                                               result[[i]][[2]][["AdaptPerTime"]][z],
                                               result[[i]][[3]][["AdaptPerTime"]][z],
                                               result[[i]][[4]][["AdaptPerTime"]][z],
                                               result[[i]][[5]][["AdaptPerTime"]][z],
                                               result[[i]][[6]][["AdaptPerTime"]][z],
                                               result[[i]][[7]][["AdaptPerTime"]][z],
                                               result[[i]][[8]][["AdaptPerTime"]][z],
                                               result[[i]][[9]][["AdaptPerTime"]][z],
                                               result[[i]][[10]][["AdaptPerTime"]][z]))}
  
for (i in 1: nrow(seq)) {
  result[[i]][["UpperAdaptPerTime"]] <- c()
  for (z in 1: (1/seq$u[i]))
  result[[i]][["UpperAdaptPerTime"]][z] <- max(c(result[[i]][[1]][["AdaptPerTime"]][z],
                                                 result[[i]][[2]][["AdaptPerTime"]][z],
                                                 result[[i]][[3]][["AdaptPerTime"]][z],
                                                 result[[i]][[4]][["AdaptPerTime"]][z],
                                                 result[[i]][[5]][["AdaptPerTime"]][z],
                                                 result[[i]][[6]][["AdaptPerTime"]][z],
                                                 result[[i]][[7]][["AdaptPerTime"]][z],
                                                 result[[i]][[8]][["AdaptPerTime"]][z],
                                                 result[[i]][[9]][["AdaptPerTime"]][z],
                                                 result[[i]][[10]][["AdaptPerTime"]][z]))}

for (i in 1: nrow(seq)) {
  result[[i]][["LowerAdaptPerTime"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["LowerAdaptPerTime"]][z] <- min(c(result[[i]][[1]][["AdaptPerTime"]][z],
                                                   result[[i]][[2]][["AdaptPerTime"]][z],
                                                   result[[i]][[3]][["AdaptPerTime"]][z],
                                                   result[[i]][[4]][["AdaptPerTime"]][z],
                                                   result[[i]][[5]][["AdaptPerTime"]][z],
                                                   result[[i]][[6]][["AdaptPerTime"]][z],
                                                   result[[i]][[7]][["AdaptPerTime"]][z],
                                                   result[[i]][[8]][["AdaptPerTime"]][z],
                                                   result[[i]][[9]][["AdaptPerTime"]][z],
                                                   result[[i]][[10]][["AdaptPerTime"]][z]))}


#Average over 5 sims
for (i in 1: nrow(seq)) {
  result[[i]][["MeanBabiesPerTime"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["MeanBabiesPerTime"]][z] <- mean(c(result[[i]][[1]][["BabiesPerTime"]][z],
                                                   result[[i]][[2]][["BabiesPerTime"]][z],
                                                   result[[i]][[3]][["BabiesPerTime"]][z],
                                                   result[[i]][[4]][["BabiesPerTime"]][z],
                                                   result[[i]][[5]][["BabiesPerTime"]][z],
                                                   result[[i]][[6]][["BabiesPerTime"]][z],
                                                   result[[i]][[7]][["BabiesPerTime"]][z],
                                                   result[[i]][[8]][["BabiesPerTime"]][z],
                                                   result[[i]][[9]][["BabiesPerTime"]][z],
                                                   result[[i]][[10]][["BabiesPerTime"]][z]))}

for (i in 1: nrow(seq)) {
  result[[i]][["UpperBabiesPerTime"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["UpperBabiesPerTime"]][z] <- max(c(result[[i]][[1]][["BabiesPerTime"]][z],
                                                    result[[i]][[2]][["BabiesPerTime"]][z],
                                                    result[[i]][[3]][["BabiesPerTime"]][z],
                                                    result[[i]][[4]][["BabiesPerTime"]][z],
                                                    result[[i]][[5]][["BabiesPerTime"]][z],
                                                    result[[i]][[6]][["BabiesPerTime"]][z],
                                                    result[[i]][[7]][["BabiesPerTime"]][z],
                                                    result[[i]][[8]][["BabiesPerTime"]][z],
                                                    result[[i]][[9]][["BabiesPerTime"]][z],
                                                    result[[i]][[10]][["BabiesPerTime"]][z]))}

for (i in 1: nrow(seq)) {
  result[[i]][["LowerBabiesPerTime"]] <- c()
  for (z in 1: (1/seq$u[i]))
    result[[i]][["LowerBabiesPerTime"]][z] <- min(c(result[[i]][[1]][["BabiesPerTime"]][z],
                                                    result[[i]][[2]][["BabiesPerTime"]][z],
                                                    result[[i]][[3]][["BabiesPerTime"]][z],
                                                    result[[i]][[4]][["BabiesPerTime"]][z],
                                                    result[[i]][[5]][["BabiesPerTime"]][z],
                                                    result[[i]][[6]][["BabiesPerTime"]][z],
                                                    result[[i]][[7]][["BabiesPerTime"]][z],
                                                    result[[i]][[8]][["BabiesPerTime"]][z],
                                                    result[[i]][[9]][["BabiesPerTime"]][z],
                                                    result[[i]][[10]][["BabiesPerTime"]][z]))}
  

# Plotting function

Babies_EffRates_plot_fct <- function(w,u, N_hat, L_hat, Mxi, c){
  par(mfrow = c(2,2),
      oma = c(4,2,1,3),
      mar = c(1, 3.5, 1.5, 1.5))
  
  L_hat <- L_hat
  
  plot(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat & seq$c==c & seq$regulation== "Survival")]][["MeanAdaptPerTime"]], type = "l",xlab = "", ylab = "", lwd=3)
  polygon(c(1:(1/u),(1/u):1), c(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Survival")]][["UpperAdaptPerTime"]],rev(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$L_hat==L_hat & seq$regulation== "Survival")]][["LowerAdaptPerTime"]])), col=alpha("black",alpha = 0.3), border = NA)
  par(new=TRUE)
  plot(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Survival")]][["MeanBabiesPerTime"]], col=col.pal[3], type = "l",yaxt="n",xlab = "", ylab = "", lwd=3, ylim = c(0, max(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$L_hat==L_hat & seq$regulation== "Survival")]][["MeanBabiesPerTime"]], ra.rm=TRUE)))
  polygon(c(1:(1/u),(1/u):1), c(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat & seq$c==c &seq$regulation== "Survival")]][["UpperBabiesPerTime"]],rev(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$L_hat==L_hat & seq$regulation== "Survival")]][["LowerBabiesPerTime"]])), col=alpha(col.pal[3],alpha = 0.3), border = NA)
  axis(side = 4, col = col.pal[3], col.axis=col.pal[3] )
  title("Mortality Regulation")
  mtext(side = 2, "Proportion of adapted individuals", line = 2.7, cex = 1)

  plot(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi &seq$N_hat==N_hat &seq$L_hat==L_hat &seq$c==c & seq$regulation== "Fertility")]][["MeanAdaptPerTime"]], type = "l",xlab = "", ylab = "", lwd=3)
  polygon(c(1:(1/u),(1/u):1), c(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi &seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Fertility")]][["UpperAdaptPerTime"]],rev(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$L_hat==L_hat & seq$regulation== "Fertility")]][["LowerAdaptPerTime"]])), col=alpha("black",alpha = 0.3), border = NA)
  par(new=TRUE)
  plot(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi &seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Fertility")]][["MeanBabiesPerTime"]], col=col.pal[3], type = "l",yaxt="n",xlab = "", ylab = "",ylim=c(0,88), lwd=3)
  polygon(c(1:(1/u),(1/u):1), c(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi &seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Fertility")]][["UpperBabiesPerTime"]],rev(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$L_hat==L_hat & seq$regulation== "Fertility")]][["LowerBabiesPerTime"]])), col=alpha(col.pal[3],alpha = 0.3), border = NA)
  axis(side = 4, col = col.pal[3], col.axis=col.pal[3] )
  title("Fertility Regulation")
  mtext(side = 4, "Number of births", col = col.pal[3], line = 2.7, cex = 1)
  


  plot(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi &seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Survival")]][["MeanSurv"]],col=col.pal[1], type = "l",xlab = "", ylab = "", lwd=3,axes = FALSE)
  polygon(c(1:(1/u),(1/u):1), c(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi &seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Survival")]][["UpperSurv"]],
                                rev(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi &seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Survival")]][["LowerSurv"]])), yaxt='n',ann=FALSE,col=alpha(col.pal[1],alpha = 0.3), border = NA)
  axis(side = 2, col = col.pal[1], col.axis=col.pal[1] )
  par(new=TRUE)
  plot(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi &seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Survival")]][["MeanFert"]], col=col.pal[3], type = "l",yaxt="n",xlab = "", ylab = "",ylim=c(0.21,0.232), lwd=3)
  polygon(c(1:(1/u),(1/u):1), c(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Survival")]][["UpperFert"]],
                                rev(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat &seq$c==c & seq$regulation== "Survival")]][["LowerFert"]])),yaxt='n', ann=FALSE, col=alpha(col.pal[3],alpha = 0.3), border = NA)
  axis(side = 4, col = col.pal[3], col.axis=col.pal[3] )
  mtext(side = 2, "Effective Survival rate", line = 2.7, cex = 1, col = col.pal[1])

  
  plot(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat & seq$c==c & seq$regulation== "Fertility")]][["MeanSurv"]],yaxt='n', ann=FALSE, col=col.pal[1], type = "l",xlab = "", ylab = "", lwd=3,axes = FALSE)
  polygon(c(1:(1/u),(1/u):1), c(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat & seq$c==c & seq$regulation== "Fertility")]][["UpperSurv"]],
                                rev(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat &  seq$c==c &seq$regulation== "Fertility")]][["LowerSurv"]])),yaxt='n', ann=FALSE,col=alpha(col.pal[1],alpha = 0.3), border = NA)
  axis(side = 2, col = col.pal[1], col.axis=col.pal[1] )
  par(new=TRUE)
  plot(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat & seq$regulation== "Fertility")]][["MeanFert"]],yaxt='n', ann=FALSE,col=col.pal[3], type = "l",yaxt="n",xlab = "", ylab = "", ylim=c(0.05,0.095), lwd=3)
  polygon(c(1:(1/u),(1/u):1), c(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat & seq$c==c &seq$regulation== "Fertility")]][["UpperFert"]],
                                rev(result[[which(seq$u==u & seq$w==w & seq$Mxi==Mxi & seq$N_hat==N_hat & seq$L_hat==L_hat & seq$c==c & seq$regulation== "Fertility")]][["LowerFert"]])), yaxt='n',ann=FALSE,col=alpha(col.pal[3],alpha = 0.3), border = NA)
  axis(side = 4, col = col.pal[3], col.axis=col.pal[3] )
  mtext(side = 4, "Effective Fertility rate", col = col.pal[3], line = 2.7, cex = 1)
  mtext(side = 1, "Time since switch in the environment", line = 2.3, cex = 1, outer = TRUE)
}  

#Choose parameter values and plot
Babies_EffRates_plot_fct(0.9, 0.01, 500, 7.5, 0.01, 0.95)