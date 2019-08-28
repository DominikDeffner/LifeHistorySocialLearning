
#Simulation code for all parameter combinations reported in the paper
# First we create grid with all 810 parameter combinations, then we define simulation function and finally we pass parameter values to
#function via mclapply, a parallelized version of lapply. You should have access to some sort of computer cluster, otherwise
#it will run for very long. 

#For main results (e.g. Fig. 2), we only save population means of xi, if you want individual-level values, population size and changes
#in the environment (as needed for Fig. 3) you can change to outcommented output files

#Create matrix with all combinations of parameter values

seq<-expand.grid(Nsim=10, Tmax=7000, K=2000,u=c(0.001, 0.01, 0.1), N_hat = c(200,350,500), beta=1.1,
                 L_hat=c(3,5,7.5),sigma=1.1, 
                 w=c(0.01, 0.1,0.5,0.9, 0.99),regulation = c("Survival", "Fertility"),c=c(0.9,0.95,0.99),  Mxi= 0.01)


#Define functions for vital rate and  mortality regulation parameters

s<- function(L_hat){(L_hat-1)/L_hat}

b<- function(s,delta, N_hat){(1-s)/(s*exp(-delta*N_hat))}

gamma <- function(s,b, delta){delta*(log(1/(s*(1+b)))/log((1-s)/(b*s)))}


#Define simulation function

sim.funct <- function(Nsim, Tmax, K, u, N_hat, beta, L_hat, sigma, regulation, w, c,Mxi){
  
  
  #Define vital rate parameters
  
  s <- s(L_hat)                   
  b <- b(s=s,delta=1/1000,N_hat=N_hat)  
  
  if (regulation == "Survival"){
    gamma<-gamma(s=s, b=b, delta=1/1000)  
    delta<-0                 
  }
  
  if (regulation == "Fertility"){
    delta<-1/1000                 
    gamma<-0                 
  }
  
  
  #Create empty list for output matrices of all N=Nsim simulation runs
  Combined_list <- list()
  
  #Loop over all N=Nsim separate simulations
  for(xsim in 1:Nsim){
    
    
    #Set up Population and output matrices
    
    N<- 500 #Starting Pop size
    
    Pop <- data.frame(id=1:N, xi=NA, Strategy=NA, Adapt=NA, Signal=NA, Age=NA,Survival=NA, Fertility=NA, Parent=NA)
    
    #Create ID vector for newborns
    ID_vector <- N+1:1e7
    
    #Define starting values
    
    Pop$Age <- sample(0:20, size = length(Pop$Age), replace = TRUE)
    Pop$Adapt <- 0
    Pop$xi <- runif(length(Pop$xi), min = 0,max = 1)
    
    #Create output vector
    xi_vector<- c()
    
    #For Fig. 3, we need individual-level data
    
    #ID_matrix <- matrix(nrow = Tmax, ncol = K)
    #Adapt_matrix<- matrix(nrow = Tmax, ncol = K)
    #Age_matrix<- matrix(nrow = Tmax, ncol = K)
    #xi_matrix<- matrix(nrow = Tmax, ncol = K)
    #Popsize <- c()
    #EnvironmentalChange <- rep(0, Tmax)


    
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
        
        #For Fig. 3
        #EnvironmentalChange[i]<-1
        
      }
      
      #Record state variables of Population
      xi_vector[i]   <- mean(Pop$xi, na.rm = TRUE)
      
      
      
      #For Fig. 3, we need individual-level data
      
      #for (x in Pop$id){
     #   ID_matrix[i,which(Pop$id==x)]<- x
     #   Adapt_matrix[i,which(Pop$id==x)]<- Pop$Adapt[which(Pop$id==x)]
     #   Age_matrix[i,which(Pop$id==x)]<- Pop$Age[which(Pop$id==x)]
     #   xi_matrix[i,which(Pop$id==x)]<- Pop$xi[which(Pop$id==x)]
     #  Popsize[i] <- length(Pop$id)
      #}
      
      
      
    }# end simulation loop
    
    #Combine matrices into one output list
    Output_list<-list(Xi= xi_vector) 
    
    #For Fig. 3, we need individual-level data
    #Output_list<-list(Adapt= Adapt_matrix, Age=Age_matrix, Xi = xi_matrix, EnvChange=EnvironmentalChange, Populationsize=PopSize) 
    
    
    #store combined output list as sublist of master output list
    Combined_list[[xsim]]<- Output_list
    
  }# end different simulations loop
  
  # Return master list as function output
  
  return(Combined_list)  
  
} #end simulation function


# pass to mclapply

library(parallel)

 result <- mclapply(
  1:nrow(seq) ,
  function(i) sim.funct(seq$Nsim[i], seq$Tmax[i], seq$K[i], seq$u[i], seq$N_hat[i], 
                        seq$beta[i], seq$L_hat[i], seq$sigma[i],seq$regulation[i],seq$w[i],
                        seq$c[i],seq$Mxi[i]),
  mc.cores=1)


