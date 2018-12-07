theta <- c(4.4168380,  0.1069816, -2.3558563, 0) #from Q6

nll.weibull <- function (theta, s) {
  
  s1 <- fatigue[which(fatigue$ro==1),,]$s
  s0 <- fatigue[which(fatigue$ro==0),,]$s
  
  #compute the negative log likelihood
  k <- 1/exp(theta[3]) #shape
  lambda1 <- exp(theta[1])*(s1 - exp(theta[4]))^theta[2] #scale
  lambda0 <- exp(theta[1])*(s0 - exp(theta[4]))^theta[2] #scale
  
  loglik0 <- -sum(dweibull(s0, shape=k, scale=lambda0, log=TRUE))
  
  loglik1 <- sum((s1/lambda1)^k)
  
  loglik <- loglik0 + loglik1
  
  loglik
}

opt.gama <- optim(theta, nll.weibull, s=fatigue$s, method="BFGS", hessian=TRUE)
opt.theta <- opt.gama$par
print(gamahat)

shape <- 1/exp(opt.theta[3])
scale <- exp(opt.theta[1])*(s-exp(opt.theta[4]))^theta[2]

N <- exp(theta[1])*(s-80)^theta[2]*log(2)^exp(theta[3])
plot(N,s)
