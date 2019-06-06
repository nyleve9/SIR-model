
##This code fits a biological model to data taken from the first 40 days of an epidemic in Kinamali,
##using a combination of direct search and descent methods and sums of squares. 

#----------SIR_ts()--------------------

# This is the SIR model function, which calculates a time series of predictions for S, I, and R:

# alpha = probability that an infected individual infects another (S to I)
# gamma = probability of recovery (I to R)
# P = initial population size
# I0 = number of initially infected individuals
# maxt = time to run the simulation

SIR_ts <- function(alpha, gamma, I0, maxt, P) {

    S <- rep(NA,maxt+1)
    I <- rep(NA,maxt+1)
    R <- rep(NA,maxt+1)
    S[1] <- P - I0
    I[1] <- I0
    R[1] <- 0
        
    for ( t in 1:maxt ) {
        S[t+1] <- S[t] - S[t] * ( 1 - (1-alpha)^I[t] )
        I[t+1] <- I[t] + S[t] * ( 1 - (1-alpha)^I[t] ) - gamma * I[t]
        R[t+1] <- R[t] + gamma * I[t]
    }
    return(list(S=S,I=I,R=R))
}
    
#-------ssq.calc()----------------

# This function will calculate the total sum of squares as SSQ(I) + SSQ(R), for given parameters:

# p = parameters alpha, gamma, and I0
# Iobs = observed # of infected individuals
# Robs = observed # of removed individuals
# maxt = time of collected data

ssq.calc <- function( p, Iobs, Robs, maxt) {
    expected <- SIR_ts(alpha=p[1], gamma=p[2], I0=p[3], maxt, P)
    Sexp <- expected$S
    Iexp <- expected$I
    Idiff <- Iobs-Iexp
    ssqI <- sum(Idiff^2)
    Rexp <- expected$R
    Rdiff <- Robs-Rexp
    ssqR <- sum(Rdiff^2)
    ssq.total <- ssqR + ssqI
    return(ssq.total)
}


#----------------Main program----------------------

# This program will use the SIR time series function, sum of squares calculating function, and
# R's native optim function to find the values of parameters alpha, gamma, and I0 that correspond
# to the best fit to the data.  It then calculates the predicted values for susceptible, infectious,
# and removed individuals using those parameters.

# Import the data and assign the values for time, S, I, and R to vectors:
setwd("C:/Users/Cheng/Desktop/QuantBIO")
SIRdata <- read.csv("student02_ass5.csv",header=TRUE)   
time <- SIRdata$t
Sobs <- SIRdata$S
Iobs <- SIRdata$I
Robs <- SIRdata$R

# Visualize the observed data in plots:
par(mfrow=c(1,3))
plot(time,Sobs,xlab="Days",ylab="observed S")
lines(time,Sobs)
plot(time,Iobs,xlab="Days",ylab="observed I")
lines(time,Iobs)
plot(time,Robs,xlab="Days",ylab="observed R")
lines(time,Robs)

# Plot observed data on the same plot:
par(mfrow=c(1,1))
plot(time,Sobs,col="green",type="l",xlab="Days",ylab="Population",ylim=c(0,6200),main="Observed data from the epidemic")
lines(time,Iobs,col="red")
lines(time,Robs,col="black")
legend("left", c("Susceptible","Infectious","Removed"),pch=1,cex=0.6,col=c("green","red","black"))


# Initialize parameters found using optim:
alpha <- 0.00007399
gamma <- 0.1887
I0 <- 0.9546
maxt <- 39   
P <- 6153

# Find the minimum SSQ:

# Initialize parameters for optim use (Use ssq profiles to choose starting parameter values):
alpha <- 0.00003
gamma <- 0.8
I0 <- 0.3
maxt <- 39   
P <- 6153

# Create a vector to store parameters to find:
par <- c(alpha,gamma,I0) 

# Run optim() to find the parameters that correspond with the best fit:
fit <- optim( par, ssq.calc, Iobs=Iobs, Robs=Robs, maxt=maxt )
fit  

# Calculate fitted values:
fitted <- SIR_ts(alpha=fit$par[1],gamma=fit$par[2],I0=fit$par[3],maxt,P)

# Use a log scale to narrow down the order of magnitude from which to find the minimum SSQ values:
alpha_range <- exp(seq(log(1e-10),log(1e-1),length.out=60))
gamma_range <- exp(seq(log(1e-10),log(1e-1),length.out=60))
I0_range <- exp(seq(log(1e-10),log(200),length.out=60))

# Use this to center in on a region of I0 ssq's, then use a log scale to find non-integer restricted values:
# I0_range <- seq(0,6200,length.out=60)

# Create a matrix to store the SSQ's for the different combinations of parameters,
# and fill matrix with SSQ's for all combinations of parameters in a given range:
ssqmat <- matrix(NA,(length(alpha_range)*length(gamma_range)*length(I0_range)),4)
rowindex <- 1
for ( a in alpha_range ) {
    for ( g in gamma_range ) {
        for (i in I0_range ) {
            ssqmat[rowindex,] <- c(a,g,i,NA)
            rowindex <- rowindex + 1
        }
    }
}

for (x in 1:length(ssqmat[,1])) {
    ssqmat[x,4] <- ssq.calc(p=ssqmat[x,1:3],Iobs,Robs,maxt)
}

# SSQ profile for alpha values:
plot(ssqmat[,1],ssqmat[,4])
plot(ssqmat[,1],ssqmat[,4],xlim=c(0,0.001),ylim=c(8000000,100000000),xlab="alpha",ylab="ssq",main="SSQ profile for alpha") # 0<alpha<0.001
plot(ssqmat[,1],ssqmat[,4],xlim=c(0,0.0001),ylim=c(8000000,60000000),xlab="alpha",ylab="ssq",main="SSQ profile for alpha") # 0.00002<alpha<0.00004

# SSQ profile for gamma values:
plot(ssqmat[,2],ssqmat[,4])
plot(ssqmat[,2],ssqmat[,4],xlim=c(0,0.5),ylim=c(8000000,60000000),xlab="gamma",ylab="ssq",main="SSQ profile for gamma")   # 0.5<gamma<1

# SSQ profile for I0 values:
plot(ssqmat[,3],ssqmat[,4],xlim=c(0,6200)) 
plot(ssqmat[,3],ssqmat[,4],xlim=c(0,1),ylim=c(1000000,10000000),xlab="I0",ylab="ssq",main="SSQ profile for I0")   # I0=100

#########################################

# Plot fitted model over the data, and show observation error using green lines:

par(mfrow=c(1,3))
plot(time,Sobs,xlab="time",ylab="S",main="Observation error fit for S")
lines(time,fitted$S,col="red")
segments(time,fitted$S,time,Sobs,col="green")

plot(time,Iobs,xlab="time",ylab="I",main="Observation error fit for I")
lines(time,fitted$I,col="red")
segments(time,fitted$I,time,Iobs,col="green")

plot(time,Robs,xlab="time",ylab="R",main="Observation error fit for R")
lines(time,fitted$R,col="red")
segments(time,fitted$R,time,Robs,col="green")


# Plot models only on the same plot:
par(mfrow=c(1,1))
days <- 0:maxt    
plot(days,fitted$S,col="green",type="l",xlab="Days",ylab="Population",ylim=c(0,6200),main="A Model of the Epidemic")
lines(days,fitted$I,col="red")
lines(days,fitted$R,col="black")
legend("left", c("Susceptible","Infectious","Removed"),pch=1,cex=0.7,col=c("green","red","black"))

###################################################

# Run model for an extended time to predict what will happen in the future:
par(mfrow=c(1,1))
t.extend <- 120
extend <- SIR_ts(alpha=fit$par[1], gamma=fit$par[2], I0=fit$par[3], t.extend, 6153)  
plot(extend$S,xlab="Days",ylab="Model Population",col="green",type="l",ylim=c(0,6200),main="Estimated Course of the Epidemic")
lines(extend$I,col="red",type="l")
lines(extend$R,col="black",type="l")
legend("right", c("Susceptible","Infectious","Removed"),pch=1,cex=0.7,col=c("green","red","black"))


# Estimate how long the epidemic will run by finding out which day the number of infected individuals 
# drops below one:
which( extend$I < 1 )    #96 days 

# Estimate how many people will utimately die by calling the last element in the "Removed" vector:
extend$R[t.extend]       #5417 people die
extend$S[t.extend]       #736 people survive

# Decide if epidemic has already peaked by using which.max to find if day of peak infected individuals
# is greater than 40:
which.max( extend$I )   #epidemic peaks at day 39

# Number of peak infected individuals as predicted by the model:
extend$I[which.max(extend$I)]     #1355 individuals
extend$R[which.max(extend$I)]
