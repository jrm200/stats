#read in the data
rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#define likelihood function to minimise
ratlikelihood <- function (theta) {
  #create subset of data since for the purposes of calculating the optimal theta since when status = 0,
  #the survival function does not depend on theta. Need to optimise this function using optim.
  t1 <- rats[which(rats$status==1),,]$time
  rx1 <- rats[which(rats$status==1),,]$rx
  
  #compute the negative log likelihood
  loglik1 <- -sum(log(dweibull(t1, shape=1/exp(theta[3]), scale=exp(theta[1] + theta[2]*rx1))))
  
  t0 <- rats[which(rats$status==0),,]$time
  rx0 <- rats[which(rats$status==0),,]$rx 
  
  a <- 1/exp(theta[3])
  b <- exp(theta[1] + theta[2]*rx0)
  
  loglik2 <- sum((t0/b)^a)
  
  loglik <- loglik1 + loglik2
  
  loglik
}

#to pick initial point, vary paramaters of theta individually and fix the others to find a rough minimum

#then to find the standard error, take the sqrt of the diag of the hessian
#5% CI = +- 1.96*(std error)

theta <- c(0,0,0)
beta0 <- theta[1]
beta1 <- theta[2]
logsig <- theta[3]
beta0var <- 0
beta1var <- 0
logsigvar <- 0

for (i in 0:1000) {
  theta[1] <- i*0.005
  beta0[i+1] <- theta[1]
  beta0var[i+1] <- ratlikelihood(theta)
}

plot(beta0, beta0var, type="l") #estimate optimal beta0 to be ~4.25

theta <- c(0,0,0)
for (i in 0:1000) {
  theta[2] <- i*0.005
  beta1[i+1] <- theta[2]
  beta1var[i+1] <- ratlikelihood(theta)
}

plot(beta1, beta1var, type="l") #estimate optimal beta1 to be ~4.25

theta <- c(0,0,0)
for (i in 0:1000) {
  theta[3] <- i*0.005
  logsig[i+1] <- theta[3]
  logsigvar[i+1] <- ratlikelihood(theta)
}

plot(logsig, logsigvar, type="l") #estimate optimal logsigma to be ~1.75

theta <- c(4.5,4.5,1.75)
ratlik <- optim(theta, ratlikelihood, hessian=TRUE, method="Nelder-Mead")
hess <- solve(ratlik$hessian)
stderr <- sqrt(diag(hess))
# this needs work dont use normal distribution
CI <- c(ratlik$par[2]-1.96*stderr[2], ratlik$par[2]+1.96*stderr[2])


x1 <- seq(from=0,to=200,by=0.5)
y1 <- dweibull(x1, shape = 1/exp(ratlik$par[3]), scale = exp(ratlik$par[1] + ratlik$par[2]))
plot(x1,y1,type="l")
x0 <- seq(from=0,to=200,by=0.5)
y0 <- dweibull(x0, shape = 1/exp(ratlik$par[3]), scale = exp(ratlik$par[1]))
plot(x0,y0,type="l")
