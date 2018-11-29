#read in the data
rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#define likelihood function to minimise
ratlikelihood <- function (theta) {
  #create subset of data since for the purposes of calculating the optimal theta since when status = 0,
  #the survival function does not depend on theta. Need to optimise this function using optim.
  t1 <- rats[which(rats$status==1),,]$time
  rx1 <- rats[which(rats$status==1),,]$rx
  t0 <- rats[which(rats$status==0),,]$time
  rx0 <- rats[which(rats$status==0),,]$rx
  
  #compute the negative log likelihood
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(theta[1] + theta[2]*rx0) #scale
  
  loglik1 <- -sum(log(dweibull(t1, shape=1/exp(theta[3]), scale=exp(theta[1] + theta[2]*rx1))))
  
  loglik2 <- sum((t0/lambda)^k)
  
  loglik <- loglik1 + loglik2
  
  loglik
}

#to pick initial point, vary paramaters of theta individually and fix the others to find a rough minimum

#initialize variables
theta <- c(0,0,0)
beta0 <- theta[1]
beta1 <- theta[2]
logsig <- theta[3]
beta0var <- 0
beta1var <- 0
logsigvar <- 0

for (i in 0:1000) {
  theta[1] <- 4.5 + i*0.0025
  beta0[i+1] <- theta[1]
  beta0var[i+1] <- ratlikelihood(theta)
}

plot(beta0, beta0var, type="l") #estimate optimal beta0 to be ~5.75

theta <- c(0,0,0)
for (i in 0:1000) {
  theta[2] <- 4 + i*0.0025
  beta1[i+1] <- theta[2]
  beta1var[i+1] <- ratlikelihood(theta)
}

plot(beta1, beta1var, type="l") #estimate optimal beta1 to be ~5.25

theta <- c(0,0,0)
for (i in 0:1000) {
  theta[3] <- 1.5 + i*0.005
  logsig[i+1] <- theta[3]
  logsigvar[i+1] <- ratlikelihood(theta)
}

plot(logsig, logsigvar, type="l") #estimate optimal logsigma to be ~2.75

theta <- c(5.75,5.25,2.75)
ratlik <- optim(theta, ratlikelihood, hessian=TRUE, method="Nelder-Mead")
weibullapprox <- ratlik$par

#then to find the standard error, take the sqrt of the diag of the inverse hessian
hess <- solve(ratlik$hessian)
stderr <- sqrt(diag(hess))

#95% CI = parameter estimate +- ~1.96*stderr
CI <- c(ratlik$par[2]+qnorm(0.025)*stderr[2], ratlik$par[2]+qnorm(0.975)*stderr[2])

a <- 1/exp(ratlik$par[3]) #shape
b <- exp(ratlik$par[1] + ratlik$par[2]) #scale (received treatment)

x1 <- seq(from=0,to=200,by=0.5)
y1 <- dweibull(x1, shape = a, scale = b)
plot(x1,y1,type="l")

b <- exp(ratlik$par[1]) #redefine scale (no treatment)
x0 <- seq(from=0,to=200,by=0.5)
y0 <- dweibull(x0, shape = a, scale = b)
plot(x0,y0,type="l")

#plots appear to show that rats which received medicine died sooner.................