
#Code for Fig. S1:  Population dynamics for adults (blue) and juveniles (red) according to analytical models

#color stuff
library(scales)
require(RColorBrewer)#load package

x <- seq(from=0, to=1, by=0.2) # fake data
col.pal <- brewer.pal(length(x), "Dark2") #create a pallette which you loop over for corresponding values


L_hat <- 7.5
delta <- 1/1000

par(mfrow=c(3,2),
    oma=c(2.5,2.5,0,0),
    mar=c(1,1,1,1)+1)
#My new approach
#Choose, N_hat, L_hat and delta
delta <- 1/1000
for (N_hat in c(200,350,500)) {
  

#Define functions for s and b

s<- function(L_hat){(L_hat-1)/L_hat}

b<- function(s,delta, N_hat){(1-s)/(s*exp(-delta*N_hat))}

#Determine gamma from delta

gamma <- function(s,b, delta){delta*(log(1/(s*(1+b)))/log((1-s)/(b*s)))}

#Define parameter values

s <- s(L_hat)                #Baseline adult survival

b <- b(s=s,delta=delta,N_hat=N_hat)                   #Baseline fertility

gamma<- gamma(s=s,b=b, delta=delta)     #Strength of Pop regulation through adult survival (0 = no regulation)  


N0 <- 1
N1<- 1


PopSize<-c()
PopSize[1]<- N0+N1

AdultSize<- c()
AdultSize[1]<- N1

JuvenileSize<- c()
JuvenileSize[1]<- N0

    #Fertility Regulation
for (i in 1:500) {
JuvenileSize[i+1] <- AdultSize[i]*b*exp(-delta*PopSize[i])

AdultSize[i+1]    <- (AdultSize[i]+JuvenileSize[i])*s
  
PopSize[i+1]<-JuvenileSize[i+1]+AdultSize[i+1]
}

plot(PopSize, type = "n", ylab="N",ylim = c(0, max(PopSize)), xlab = "Timestep")
lines(PopSize, col=col.pal[6], lwd=3)
lines(JuvenileSize, col= col.pal[1], lwd=3)
lines(AdultSize, col= col.pal[4], lwd=3)
if (N_hat== 200) title("Fertility Regulation")
if (N_hat== 200) legend("topleft", c("Juveniles", "Adults", "All"), col=c(col.pal[1], col.pal[4], col.pal[6]),lwd=3, lty=1, cex=1, bty="n")


N0 <- 1
N1<- 1

PopSize<-c()
PopSize[1]<- N0+N1

AdultSize<- c()
AdultSize[1]<- N1

JuvenileSize<- c()
JuvenileSize[1]<- N0

#Mortality Regulation
for (i in 1:500) {
  JuvenileSize[i+1] <- AdultSize[i]*b
  
  AdultSize[i+1]    <- (AdultSize[i]+JuvenileSize[i])*s*exp(-gamma*PopSize[i])
  
  PopSize[i+1]<-JuvenileSize[i+1]+AdultSize[i+1]
}

plot(PopSize, type = "n", ylab="N",ylim = c(0, max(PopSize)), xlab = "Timestep")
lines(PopSize, col=col.pal[6], lwd=3)
lines(JuvenileSize, col= col.pal[1], lwd=3)
lines(AdultSize, col= col.pal[4], lwd=3)
if (N_hat== 200) title("Mortality Regulation")
mtext("Timestep", side = 1, outer = TRUE, line=1.2)
mtext("Population Size N", side = 2, outer = TRUE, line=1.2)

}
