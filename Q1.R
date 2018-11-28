#read in the data
rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#define likelihood function to minimise
ratlikelihood <- function (theta) {
  #create subset of data since for the purposes of calculating the optimal theta since when status = 0,
  #the survival function does not depend on theta. Need to optimise this function using optim.
  t <- rats[which(rats$status==1),,]$time
  
  #compute the negative log likelihood
  loglik <- -sum(log(dweibull(t, shape=1/exp(theta[3]), scale=exp(theta[1] + theta[2]))))
  loglik
}

#to pick initial point, vary paramaters of theta individually and fix the others to find a rough minimum

#then to find the standard error, take the sqrt of the diag of the hessian
#5% CI = +- 1.96*(std error)

theta <- c(0,0,0)
beta0 <- theta[1]
beta1 <- theta[2]
logsig <- theta[3]

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
