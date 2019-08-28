

# Simulation code for single run and plotting code to reproduce Fig. S3 (for Tmax=1000)

#Choose, N_hat, L_hat and delta
N_hat <- 500
L_hat <- 7.5
delta <- 1/1000

#Define functions for s and b

s<- function(L_hat){(L_hat-1)/L_hat}
b<- function(s,delta, N_hat){(1-s)/(s*exp(-delta*N_hat))}

#Determine gamma from delta

gamma <- function(s,b, delta){delta*(log(1/(s*(1+b)))/log((1-s)/(b*s)))}

#Define vital rate parameters

s <- s(L_hat)                           #Baseline survival
b <- b(s=s,delta=delta,N_hat=N_hat)     #Baseline fertility


#Select mode of population regulation

Regulation <- "Fertility" # Set to "Fertility" or "Mortality"


if (Regulation=="Fertility"){
    gamma<- 0
} else {
    gamma <- gamma(s=s,b=b, delta=delta)     
    delta=0
}
  
#Define other parameter values

  Tmax=1000                     #Simulation length
  K=2000                        #Maximum size of Pop
  beta<-1.1                     #Fertility of adapted individuals B=b*beta
  sigma<-  1.1                  #Survival of adapted adults S = s * sigma

  w<-0.9                        #P(Newborn innovation produces adaptive behavior)
  u<-0.01                       #P(Environmental change)
  c<-0.9                       #Recruitment cost for IL
  
  Mxi<-0.01                     #Mutation rate for learning strategy

  
      #Set up Population 
  
      N<- 500 #Starting Pop size
  
      Pop <- data.frame(id=1:N, xi=NA, Strategy=NA, Adapt=NA, Signal=NA, Age=NA,Survival=NA, Fertility=NA, Parent=NA)
      
      #Create ID vector for juveniles
      ID_vector <- N+1:5e6
      
      #Define starting values
      
      Pop$Age <- sample(0:20, size = length(Pop$Age), replace = TRUE)
      Pop$Adapt <- 0
      Pop$xi <- runif(length(Pop$xi), min = 0,max = 1)

      #Create output vectors/matrices
      ID_matrix <- matrix(nrow = Tmax, ncol = K)
      Adapt_matrix<- matrix(nrow = Tmax, ncol = K)
      xi_matrix<- matrix(nrow = Tmax, ncol = K)
      Age_matrix<- matrix(nrow = Tmax, ncol = K)
      EnvironmentalChange <- rep(0, Tmax)
      TimesSinceChange <- NA
      PopSize<-NA
      
      #Start simulation loop
      for (i in 1:Tmax){
        
        # 1.Reproduction
        # Determine fertility
        
        Pop$Fertility[which(Pop$Adapt==1)]<-b*beta*exp(-delta*N)
        Pop$Fertility[which(Pop$Adapt==0)]<-b*exp(-delta*N)
        
        #Select parents for reproduction
        
        Pop$Parent <- rbinom(length(Pop$Parent), 1, Pop$Fertility)
        
        #Produce offspring
        Newborns<-subset(Pop, Pop$Parent==1)
        
        #Sample new IDs for newborns from ID vector
        Newborns$id<-sample(ID_vector, size = length(Newborns$id), replace = FALSE)
        
        #Remove IDs from ID vector
        ID_vector <- ID_vector[! ID_vector %in% Newborns$id]
        
        #Assign other values
        Newborns$Age<-0
        Newborns$Adapt<-0
        Newborns$Strategy<-NA
        Newborns$Survival<-NA
        Newborns$Fertility<-NA
        Newborns$Parent<-NA
        
        #Mutation on xi; phenotypes change continually
          OldXi <- Newborns$xi
          Newborns$xi <- OldXi + rnorm(length(OldXi), mean = 0, sd=Mxi)
          OutOfRange <- which(Newborns$xi < 0 | Newborns$xi > 1)
          if (length(OutOfRange)==1){
            Old <- Newborns$xi[OutOfRange]
            Newborns$xi[OutOfRange] <- Old + rnorm(1, mean = 0, sd=Mxi)
            while (Newborns$xi[OutOfRange] < 0 | Newborns$xi[OutOfRange] >1) {
              Newborns$xi[OutOfRange] <- Old + rnorm(1, mean = 0, sd=Mxi)
            }
            
          } else {
             for (x_mut in OutOfRange){
                Old <- Newborns$xi[x_mut]
                Newborns$xi[x_mut] <- Old + rnorm(1, mean = 0, sd=Mxi)
                while (Newborns$xi[x_mut] < 0 | Newborns$xi[x_mut] >1) {
                Newborns$xi[x_mut] <- Old + rnorm(1, mean = 0, sd=Mxi)
            }
           }
          }
          
        #Add offspring to Pop
        Pop<-rbind(Pop,Newborns)
        
        #Updating Pop size
        N<-length(Pop$id)
        
        
        #2. Learning
        
        #Assign Learning strategy based on value of xi
        for (x in Pop$id[which(Pop$Age== 0)]) {
          Pop$Strategy[which(Pop$id==x)] <- sample(c("Individual", "Social"), size = 1, prob = c(Pop$xi[which(Pop$id==x)],1-Pop$xi[which(Pop$id==x)]))
        } 
        
        #Adults don't learn
        Pop$Strategy[which(Pop$Age>0)]<-"None"
        
        #2.1.Social learning
        
        #Select potential models; Juveniles cannot copy other juveniles as no one has any adaptive information yet
        ModelPool<-subset(Pop,Pop$Age>0)
        
        #Loop over all juveniles
        for(x_soc in Pop$id[which(Pop$Strategy=="Social")]){
        Pop$Adapt[which(Pop$id==x_soc)]<-sample(ModelPool$Adapt,size = 1)
        }

        #1.1.2 Innovation
        Pop$Adapt[which(Pop$Strategy=="Individual")]<- rbinom(length(Pop$Adapt[which(Pop$Strategy=="Individual")]), 1, w)

        # 3. Survival 
        
        #Adults and socially learning juveniles
        Pop$Survival[which(Pop$Adapt==1 & Pop$Strategy!="Individual")]<- rbinom(length(Pop$id[which(Pop$Adapt==1 & Pop$Strategy!="Individual")]),1,(sigma*s)*exp(-gamma*N))
        Pop$Survival[which(Pop$Adapt==0 & Pop$Strategy!="Individual")]<- rbinom(length(Pop$id[which(Pop$Adapt==0 & Pop$Strategy!="Individual")]),1,s*exp(-gamma*N))
        
        #Innovating juveniles
        Pop$Survival[which(Pop$Adapt==1 & Pop$Strategy=="Individual")]<- rbinom(length(Pop$id[which(Pop$Adapt==1 & Pop$Strategy=="Individual")]),1,(c*sigma*s)*exp(-gamma*N))
        Pop$Survival[which(Pop$Adapt==0 & Pop$Strategy=="Individual")]<- rbinom(length(Pop$id[which(Pop$Adapt==0 & Pop$Strategy=="Individual")]),1,(c*s)*exp(-gamma*N))
        
        
        #Retain survivors in the Population
        Pop <-subset(Pop,Pop$Survival==1)
        
        #Ageing
        #All survivors age
        Pop$Age<-Pop$Age+1
        
        
        #4.Environmental change
        ProbEnvChange<-runif(1)
        
        if(ProbEnvChange<u){
          Pop$Adapt<-0
          EnvironmentalChange[i]<-1
        }
        
        #Record state variables of Population
        
        for (x in Pop$id){
          ID_matrix[i,which(Pop$id==x)]<- x
          Adapt_matrix[i,which(Pop$id==x)]<- Pop$Adapt[which(Pop$id==x)]
          Age_matrix[i,which(Pop$id==x)]<- Pop$Age[which(Pop$id==x)]
          xi_matrix[i,which(Pop$id==x)]<- Pop$xi[which(Pop$id==x)]
        }
        
        PopSize[i]<- N
        print(i)
      }# end simulation loop
      
      
      #Calculate mean adaptation levels for selected age classes
      MeanAdapt<-matrix(nrow = 1000, ncol = 10)
      for (x in 1:1000) {
        for (y in 1:10) {
          MeanAdapt[x,y]<- mean(Adapt_matrix[x,][which(Age_matrix[x,]==y)])
        }
      }
      
      #Plotting script for Fig. S3
      par(mfrow=c(1,1),
          oma=c(3,3.5,0,3.5),
          mar=c(0,0,1,0))
      plot(MeanAdapt[1:1000,1], type = "n", ylim = c(min(PopSize[1:900]),max(PopSize[1:900])), xlab = "", ylab = "")
      lines(PopSize[1:900],lwd = 1)
      
      par(new=TRUE)
      library(scales)
      plot(MeanAdapt[1:1000,1], col=alpha("orange", alpha = 0.3), ylim = c(0,1), type = "n", yaxt="n", ylab = "", xlab = "")
      
      lines(MeanAdapt[1:900,1], col=alpha("orange", alpha = 0.3))
      lines(MeanAdapt[1:900,3], col=alpha("red", alpha = 0.3))
      lines(MeanAdapt[1:900,5], col=alpha("blue", alpha = 0.3))
      lines(MeanAdapt[1:900,7], col=alpha("pink", alpha = 0.3))
      lines(MeanAdapt[1:900,10], col=alpha("green", alpha = 0.3))
      
      abline(v=which(EnvironmentalChange==1)[1:5], lty=2, lwd=2)
      axis(side = 4)
      
      legend("topright", title = "Age classes", c("1","3", "5","7", "10"), col=c("orange", "red", "blue", "pink", "green"), lwd=3,cex=0.8, bty="n",lty=1)
      legend("right", "Population size", col="black", cex=0.8, bty="n",lty=1, lwd=3)
      mtext("Proportion of adapted individuals", side = 4, line = 2.3)
      mtext("Timestep", side = 1, line = 2)
      mtext("Population size", side = 2, line = 2.3)