
#Code for Figs. 2, S4 and S5
#Takes results from file "Simulation_Code_Cluster"

#color stuff
library(scales)
require(RColorBrewer)#load package


x <- seq(from=0, to=1, by=0.2) # fake data
col.pal <- brewer.pal(length(x), "Dark2") #create a pallette which you loop over for corresponding values

#Exclude first 2000 timesteps

for (x in 1:nrow(seq)){
  for (y in 1:10){
    result[[x]][[y]][["xi"]]<- result[[x]][[y]][["xi"]][-(1:2000)]
  }
}

#Plotting function that computes averages over all values of other parameters and plots results from 10 simulations and their average

plot_fct <- function(c){
  
  par(mfrow = c(2,2),
      oma = c(0,2.5,0,0) + 0.1,
      mar = c(4, 2, 1, 1))
  
  # 
  #
  # N_hat
  #
  #
  
  MeansSmallMort<-matrix(NA, nrow = 10, ncol = length(which(seq$N_hat==200 & seq$regulation=="Survival"& seq$c==c)))
  x_seq <- which(seq$N_hat==200 & seq$regulation=="Survival" & seq$c==c)
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansSmallMort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallSmallMort<- apply(MeansSmallMort, 2, mean)
  
  MeansSmallFert<-matrix(NA, nrow = 10, ncol = length(which(seq$N_hat==200 & seq$regulation=="Fertility"& seq$c==c)))
  x_seq <- which(seq$N_hat==200 & seq$regulation=="Fertility"& seq$c==c)
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansSmallFert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallSmallFert<- apply(MeansSmallFert, 2, mean)
  
  MeansIntermMort<-matrix(NA, nrow = 10, ncol = length(which(seq$N_hat==350 & seq$regulation=="Survival"& seq$c==c)))
  x_seq <- which(seq$N_hat==350 & seq$regulation=="Survival"& seq$c==c)
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansIntermMort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallIntermMort<- apply(MeansIntermMort, 2, mean)
  
  MeansIntermFert<-matrix(NA, nrow = 10, ncol = length(which(seq$N_hat==350 & seq$regulation=="Fertility"& seq$c==c)))
  x_seq <- which(seq$N_hat==350 & seq$regulation=="Fertility"& seq$c==c)
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansIntermFert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallIntermFert<- apply(MeansIntermFert, 2, mean)
  
  
  MeansLargeMort<-matrix(NA, nrow = 10, ncol = length(which(seq$N_hat==500 & seq$regulation=="Survival"& seq$c==c)))
  x_seq <- which(seq$N_hat==500 & seq$regulation=="Survival"& seq$c==c)
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansLargeMort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallLargeMort<- apply(MeansLargeMort, 2, mean)
  
  MeansLargeFert<-matrix(NA, nrow = 10, ncol = length(which(seq$N_hat==500 & seq$regulation=="Fertility"& seq$c==c)))
  x_seq <- which(seq$N_hat==500 & seq$regulation=="Fertility"& seq$c==c)
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansLargeFert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallLargeFert<- apply(MeansLargeFert, 2, mean)

  plot(c(mean(MeansSmallFert[y,]),mean(MeansIntermFert[y,]), mean(MeansLargeFert[y,])), type = "n", ylim = c(0,0.6), ylab = "", xaxt = "n", xlab=expression(hat(N)) )
  Overallmean <- c()
  for (y in 1:10) {
    
    lines(c(mean(MeansSmallFert[y,]),mean(MeansIntermFert[y,]), mean(MeansLargeFert[y,])), col=alpha(col.pal[3], alpha = 0.4))
    lines(c(mean(MeansSmallMort[y,]),mean(MeansIntermMort[y,]), mean(MeansLargeMort[y,])), col=alpha(col.pal[1], alpha = 0.4))
    
  }
  lines(c(mean(OverallSmallFert),mean(OverallIntermFert), mean(OverallLargeFert)), col=alpha(col.pal[3], alpha = 1), lwd=2)
  lines(c(mean(OverallSmallMort),mean(OverallIntermMort), mean(OverallLargeMort)), col=alpha(col.pal[1], alpha = 1), lwd=2)
  
  axis(1, at=1:3, labels=c(200,350,500))
  legend("topleft", "A", bty = "n")

  #
  #
  # L_hat
  #
  #

  x_seq <- which(seq$L_hat==3 & seq$regulation=="Survival"& seq$c==c)
  MeansFastMort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansFastMort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallFastMort<- apply(MeansFastMort, 2, mean)
  
  x_seq <- which(seq$L_hat==3 & seq$regulation=="Fertility"& seq$c==c)
  MeansFastFert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansFastFert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallFastFert<- apply(MeansFastFert, 2, mean)

  x_seq <- which(seq$L_hat==5 & seq$regulation=="Survival"& seq$c==c)
  MeansIntermMort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansIntermMort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallIntermMort<- apply(MeansIntermMort, 2, mean)
  
  x_seq <- which(seq$L_hat==5 & seq$regulation=="Fertility"& seq$c==c)
  MeansIntermFert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansIntermFert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallIntermFert<- apply(MeansIntermFert, 2, mean)

  x_seq <- which(seq$L_hat==7.5 & seq$regulation=="Survival"& seq$c==c)
  MeansSlowMort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansSlowMort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  
  OverallSlowMort<- apply(MeansSlowMort, 2, mean)
  
  x_seq <- which(seq$L_hat==7.5 & seq$regulation=="Fertility"& seq$c==c)
  MeansSlowFert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansSlowFert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallSlowFert<- apply(MeansSlowFert, 2, mean)
  
  plot(c(mean(MeansFastFert[y,]),mean(MeansIntermFert[y,]), mean(MeansSlowFert[y,])), type = "n", ylim = c(0,0.6), ylab = "", xaxt = "n", xlab= expression(hat(L)) )
  Overallmean <- c()
  for (y in 1:10) {
    
    lines(c(mean(MeansFastFert[y,]),mean(MeansIntermFert[y,]), mean(MeansSlowFert[y,])), col=alpha(col.pal[3], alpha = 0.4))
    lines(c(mean(MeansFastMort[y,]),mean(MeansIntermMort[y,]), mean(MeansSlowMort[y,])), col=alpha(col.pal[1], alpha = 0.4))
    
  }
  lines(c(mean(OverallFastFert),mean(OverallIntermFert), mean(OverallSlowFert)), col=alpha(col.pal[3], alpha = 1), lwd=2)
  lines(c(mean(OverallFastMort),mean(OverallIntermMort), mean(OverallSlowMort)), col=alpha(col.pal[1], alpha = 1), lwd=2)
  
  axis(1, at=1:3, labels=c(3,5,7.5))
  legend("topright", c("Mortality Regulation", "Fertility Regulation"), col=c(col.pal[1],col.pal[3]), lty=1, lwd=2, bty = "n", cex = 1)
  legend("topleft", "B", bty = "n")

  #
  #
  # Omega
  #
  #
  
  x_seq <- which(seq$u==0.001 & seq$regulation=="Survival"& seq$c==c)
  MeansSlowMort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansSlowMort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallSlowMort<- apply(MeansSlowMort, 2, mean)
  
  x_seq <- which(seq$u==0.001 & seq$regulation=="Fertility"& seq$c==c)
  MeansSlowFert <-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansSlowFert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallSlowFert<- apply(MeansSlowFert, 2, mean)

  x_seq <- which(seq$u==0.01 & seq$regulation=="Survival"& seq$c==c)
  MeansIntermMort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansIntermMort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallIntermMort<- apply(MeansIntermMort, 2, mean)
  
  x_seq <- which(seq$u==0.01 & seq$regulation=="Fertility"& seq$c==c)
  MeansIntermFert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansIntermFert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallIntermFert<- apply(MeansIntermFert, 2, mean)

  x_seq <- which(seq$u==0.1 & seq$regulation=="Survival"& seq$c==c)
  MeansFastMort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansFastMort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallFastMort<- apply(MeansFastMort, 2, mean)
  
  x_seq <- which(seq$u==0.1 & seq$regulation=="Fertility"& seq$c==c)
  MeansFastFert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      MeansFastFert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  OverallFastFert<- apply(MeansFastFert, 2, mean)
  
  plot(c(mean(MeansSlowFert[y,]),mean(MeansIntermFert[y,]), mean(MeansFastFert[y,])), type = "n", ylim = c(0,0.6), ylab = "", xaxt = "n", xlab=expression(Omega) )
  Overallmean <- c()
  for (y in 1:10) {
    
    lines(c(mean(MeansSlowFert[y,]),mean(MeansIntermFert[y,]), mean(MeansFastFert[y,])), col=alpha(col.pal[3], alpha = 0.4))
    lines(c(mean(MeansSlowMort[y,]),mean(MeansIntermMort[y,]), mean(MeansFastMort[y,])), col=alpha(col.pal[1], alpha = 0.4))
    
  }
  lines(c(mean(OverallSlowFert),mean(OverallIntermFert), mean(OverallFastFert)), col=alpha(col.pal[3], alpha = 1), lwd=2)
  lines(c(mean(OverallSlowMort),mean(OverallIntermMort), mean(OverallFastMort)), col=alpha(col.pal[1], alpha = 1), lwd=2)
  
  axis(1, at=1:3, labels=c(1000,100,10))
  legend("topleft", "C", bty = "n")

  #
  #
  # w
  #
  #
  
  x_seq <- which(seq$w==0.01 & seq$regulation=="Survival"& seq$c==c)
  Means01Mort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means01Mort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall01Mort<- apply(Means01Mort, 2, mean)
  
  x_seq <- which(seq$w==0.01 & seq$regulation=="Fertility"& seq$c==c)
  Means01Fert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means01Fert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall01Fert<- apply(Means01Fert, 2, mean)
  
  x_seq <- which(seq$w==0.1 & seq$regulation=="Survival"& seq$c==c)
  Means1Mort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means1Mort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall1Mort<- apply(Means1Mort, 2, mean)
  
  x_seq <- which(seq$w==0.1 & seq$regulation=="Fertility"& seq$c==c)
  Means1Fert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means1Fert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall1Fert<- apply(Means1Fert, 2, mean)

  x_seq <- which(seq$w==0.5 & seq$regulation=="Survival"& seq$c==c)
  Means5Mort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means5Mort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall5Mort<- apply(Means5Mort, 2, mean)
  
  x_seq <- which(seq$w==0.5 & seq$regulation=="Fertility"& seq$c==c)
  Means5Fert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means5Fert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall5Fert<- apply(Means5Fert, 2, mean)

  x_seq <- which(seq$w==0.9 & seq$regulation=="Survival"& seq$c==c)
  Means9Mort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means9Mort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall9Mort<- apply(Means9Mort, 2, mean)
  
  x_seq <- which(seq$w==0.9 & seq$regulation=="Fertility"& seq$c==c)
  Means9Fert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means9Fert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall9Fert<- apply(Means9Fert, 2, mean)
  
  x_seq <- which(seq$w==0.99 & seq$regulation=="Survival"& seq$c==c)
  Means99Mort<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means99Mort[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall99Mort<- apply(Means99Mort, 2, mean)
  
  x_seq <- which(seq$w==0.99 & seq$regulation=="Fertility"& seq$c==c)
  Means99Fert<-matrix(NA, nrow = 10, ncol = length(x_seq))
  for (x in x_seq) {
    for (y in 1:10) {
      
      Means99Fert[y,which(x_seq==x)]<- mean(result[[x]][[y]][["xi"]])
      
    }
  }
  Overall99Fert<- apply(Means99Fert, 2, mean)
  
  plot(c(mean(Means01Fert[y,]),mean(Means1Fert[y,]),mean(Means5Fert[y,]), mean(Means9Fert[y,]),mean(Means99Fert[y,])),type = "n", ylim = c(0,0.6), ylab = "", xaxt = "n", xlab="w")
  Overallmean <- c()
  for (y in 1:10) {
    
    lines(c(mean(Means01Fert[y,]),mean(Means1Fert[y,]),mean(Means5Fert[y,]), mean(Means9Fert[y,]),mean(Means99Fert[y,])), col=alpha(col.pal[3], alpha = 0.4))
    lines(c(mean(Means01Mort[y,]),mean(Means1Mort[y,]),mean(Means5Mort[y,]), mean(Means9Mort[y,]),mean(Means99Mort[y,])), col=alpha(col.pal[1], alpha = 0.4))
    
  }
  lines(c(mean(Overall01Fert),mean(Overall1Fert),mean(Overall5Fert), mean(Overall9Fert),mean(Overall99Fert)), col=alpha(col.pal[3], alpha = 1), lwd=2)
  lines(c(mean(Overall01Mort),mean(Overall1Mort),mean(Overall5Mort), mean(Overall9Mort),mean(Overall99Mort)), col=alpha(col.pal[1], alpha = 1), lwd=2)
  
  axis(1, at=1:5, labels=c(0.01,0.1,0.5,0.9, 0.99))
  legend("topleft", "D", bty = "n")
  mtext(side = 2, expression(paste("Propensity for Individual Learning  ", xi)),line=1, outer = TRUE, cex = 1)
}

plot_fct(0.95)




