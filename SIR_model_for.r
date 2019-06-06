# Brett Melbourne
# 2 Oct 2013
# brett.melbourne@colorado.edu

# Discrete time SIR model with mass action force of infection.
# Individuals become infected by other individuals with probability
# alpha, and recover with probability gamma (i.e. in average time 1/gamma).
# The deterministic model is:
#   S[t+1] <- S[t] - S[t]*(1-(1-alpha)^I[t])
#   I[t+1] <- I[t] + S[t]*(1-(1-alpha)^I[t]) - gamma*I[t]
#   R[t+1] <- R[t] + gamma*I[t]

# This version of the program uses counter controlled repetition. There
# are several advantages of this approach:
# 1) Entire storage vectors are initialized at once - no concatonation.
#    This is much faster, an issue for model fitting.
# 2) The equations read very naturally because of the indexing.


# Initialize parameters
alpha <- 0.005  #Probability that an infected individual infects another
gamma <- 0.5    #Probability of recovery (1/mean duration of illness)
P <- 500        #Total population size
I0 <- 1         #Number of initially infected individuals
maxt <- 23      #Maximum time to run the simulation

# Initialize storage and variables
S <- rep(NA,maxt+1)   #make sure you think about why we need to +1 here
I <- rep(NA,maxt+1)   #to prevent "off-by-1" error
R <- rep(NA,maxt+1)
S[1] <- P - I0
I[1] <- I0
R[1] <- 0

# Run the model (step through time)
for ( t in 1:maxt ) {
    S[t+1] <- S[t] - S[t] * ( 1 - (1-alpha)^I[t] )
    I[t+1] <- I[t] + S[t] * ( 1 - (1-alpha)^I[t] ) - gamma * I[t]
    R[t+1] <- R[t] + gamma * I[t]
}

# Plot the variables through time
# This is an illustration of some of R's graphing functions
# See also R_quick_reference.txt - PLOTTING
days <- 0:t    #nb the object t still exists after exiting the for loop
plot(NA,NA,type="n",xlim=c(0,maxt),ylim=c(0,P),xlab="Day",ylab="N")
lines(days,S,col="green")
points(days,S,col="green")
lines(days,R,col="black")
points(days,R,col="black")
lines(days,I,col="red")
points(days,I,col="red")
legend("right", c("Susceptible","Infectious","Removed"),pch=1,
       col=c("green","red","black"))


# We could replot this with different options for line width and
# symbol fills (nb I prefer the version above; below is just for
# illustration).
#windows() #Start a new device on Windows (not necessary in Rstudio)
#quartz() #Start a new device in Mac
#this step below just sets up the empty box, label, axes
plot(NA,NA,type="n",xlim=c(0,maxt),ylim=c(0,P),xlab="Day",ylab="N")
lines(days,S,col="green",lwd=2)   #just plots the green line
points(days,S,col="green",lwd=2,pch=16)  #just plots the points
lines(days,R,col="black",lwd=2,pch=16)
points(days,R,col="black",lwd=2,pch=16)
lines(days,I,col="red",lwd=2,pch=16)
points(days,I,col="red",lwd=2,pch=16)
legend("right", c("Susceptible","Infectious","Removed"),pch=rep(16,3),
       col=c("green","red","black"))   #ggplot will set up legend colors automatically


#Alternatively we could use matplot (plot columns of a matrix)
matplot(days,cbind(S,I,R),xlim=c(0,max(days)),ylim=c(0,P),
        xlab="Day",ylab="N",pch=1,col=c("green","red","black"))
matlines(days,cbind(S,I,R),col=c("green","red","black"),lty=1)
legend("right", c("Susceptible","Infectious","Removed"),pch=1,col=3:1)

