#Code to recreate Figure S2:  Isoclines for different equilibrium population sizes

gradient_maker <- function(start=NA, stop=NA, cols=c("darkorange", "white", "darkcyan"), vis=FALSE, n=1000){
  if(is.na(start) | is.na(stop)) stop("need to specify start and stop points on a numerical scale")
  colfunc <- colorRampPalette(cols)
  color.list <- colfunc(n)
  color.locations <- seq(start, stop, length=n)
  names(color.locations) <- color.list
  if(vis==TRUE) plot(color.locations, rep(1, n), col=color.list, pch="|", ylim=c(0.9, 1.1), cex=5)
  return(color.locations)
}

data_gradient <- function(data, colors=c("darkorange", "white", "darkcyan"), my.start=NA, my.stop=NA){
  if(is.na(my.start)) my.start <- min(data, na.rm=TRUE)
  if(is.na(my.stop)) my.stop <- max(data, na.rm=TRUE)
  my.gradient <- gradient_maker(start=my.start, stop=my.stop, cols=colors)
  if(any(data > max(my.gradient), na.rm=T) | any(data < min(my.gradient), na.rm=T)) warning("data is not within gradient range")
  data.colors <- rep(NA, length(data))
  for(i in 1:length(data)){
    if(!is.na(data[i])) data.colors[i] <- names(my.gradient)[which.min(abs(data[i]-my.gradient))]
  }
  data.colors
}


N_hat_fertility <- function(s,b,delta){log((1-s)/(b*s))/(-delta)} 

N_hat_mortality <- function(s,b,gamma){log(1/(s*(1+b)))/(-gamma)} 

b_fertility <- function(s) {(1-s)/(exp(-delta*N)*s)}
b_mortality <- function(s) {(1/(exp(-gamma*N)*s)-1)}

range01 <- function(x, xmin = min(x), xmax = max(x)){
  y <- (x-xmin)/(xmax-xmin)
  y[y < 0] <- 0
  y[y > 1] <- 1
  return(y)
}

par(mfrow = c(2,2),
    oma = c(2,2.5,0,0),
    mar = c(2, 2, 0.3, 0.3))

#Fertility regulation
delta=1/1500

b_seq<-seq(from=0.0001, to=1, by=0.001)
s_seq<-seq(from=0.5, to=0.9995, by=0.0005)

d<- expand.grid(b_seq, s_seq)
colnames(d) <- c("b", "s")
d$z <- N_hat_fertility(d$s, d$b, delta = delta)
d$z01 <- range01(d$z, 0,1000)
#d$color <- grey(1-d$z01)
d$color<- data_gradient(d$z01, colors = c("white", "lightskyblue","skyblue", "deepskyblue","darkblue"))


plot(d$s, d$b, col=d$color, pch=15,xlab = "",ylab = "", bty="n")
legend("bottomleft", "A", bty = "n", cex=1.3)


delta=1/550

b_seq<-seq(from=0.0001, to=1, by=0.001)
s_seq<-seq(from=0.5, to=0.9995, by=0.0005)

d<- expand.grid(b_seq, s_seq)
colnames(d) <- c("b", "s")
d$z <- N_hat_fertility(d$s, d$b, delta = delta)
d$z01 <- range01(d$z, 0,1000)
#d$color <- grey(1-d$z01)
d$color<- data_gradient(d$z01, colors = c("white", "lightskyblue","skyblue", "deepskyblue","darkblue"))


plot(d$s, d$b, col=d$color, pch=15,xlab = "",ylab = "", bty="n")

legend("bottomleft", "B", bty = "n", cex=1.3)


#Mortality regulation
gamma=1/3000

b_seq<-seq(from=0.0001, to=1, by=0.001)
s_seq<-seq(from=0.5, to=0.9995, by=0.0005)

d<- expand.grid(b_seq, s_seq)
colnames(d) <- c("b", "s")
d$z <- N_hat_mortality(d$s, d$b, gamma = gamma)
d$z01 <- range01(d$z, 0,1000)
#d$color <- grey(1-d$z01)
d$color<- data_gradient(d$z01, colors = c("white", "lightskyblue","skyblue", "deepskyblue","darkblue"))

plot(d$s, d$b, col=d$color, pch=15,xlab = "",ylab = "", bty="n")

legend("bottomleft", "C", bty = "n", cex=1.3)

gamma=1/1000

b_seq<-seq(from=0.0001, to=1, by=0.005)
s_seq<-seq(from=0.5, to=0.9995, by=0.0005)

d<- expand.grid(b_seq, s_seq)
colnames(d) <- c("b", "s")
d$z <- N_hat_mortality(d$s, d$b, gamma = gamma)
d$z01 <- range01(d$z, 0,1000)
#d$color <- grey(1-d$z01)
d$color<- data_gradient(d$z01, colors = c("white", "lightskyblue","skyblue", "deepskyblue","darkblue"))


plot(d$s, d$b, col=d$color, pch=15,xlab = "",ylab = "", bty="n")

legend("bottomleft", "D", bty = "n", cex=1.3)

mtext(side = 1, expression(paste("Survival rate ",italic("s"))),line = 1, cex = 1.1, outer = TRUE)
mtext(side = 2, expression(paste("Fertility rate ", italic("b"))), line = 1, cex = 1.1, outer = TRUE)





