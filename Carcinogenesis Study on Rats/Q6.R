fatigue<- read.table("http://people.bath.ac.uk/kai21/ASI/fatigue.txt")

gama <- 80
theta <- c(3.5,0.3,-3)
s <- fatigue$s

nll.weibull <- function (theta, gama, s) {
  
  s1 <- fatigue[which(fatigue$ro==1),,]$s
  s0 <- fatigue[which(fatigue$ro==0),,]$s
  
  #compute the negative log likelihood
  k <- 1/exp(theta[3]) #shape
  lambda1 <- exp(theta[1])*(s1 - gama)^theta[2] #scale
  lambda0 <- exp(theta[1])*(s0 - gama)^theta[2] #scale
  
  loglik0 <- -sum(dweibull(s0, shape=k, scale=lambda0, log=TRUE))
  
  loglik1 <- sum((s1/lambda1)^k)
  
  loglik <- loglik0 + loglik1
  
  loglik
}

maxlik <- optim(theta, nll.weibull, gama=gama, s=s, method="BFGS", hessian=TRUE)
theta <- maxlik$par
se <- sqrt(diag(solve(maxlik$hessian)))

CI1 <- c(maxlik$par[1]+qnorm(0.025)*se[1], maxlik$par[1]+qnorm(0.975)*se[1])
CI2 <- c(maxlik$par[2]+qnorm(0.025)*se[2], maxlik$par[2]+qnorm(0.975)*se[2])
CI3 <- c(maxlik$par[3]+qnorm(0.025)*se[3], maxlik$par[3]+qnorm(0.975)*se[3])
CI <- matrix(c(CI1,CI2,CI3), ncol=3, nrow=2)

f.vals <- rep(0,80)
f.par <- matrix(0,nrow=3, ncol=80)
for (i in 1:80) {
  gama <- i
  f.vals[i] <- nll.weibull(theta,gama,s)
  maxlik <- optim(theta, nll.weibull, gama=gama, s=s, hessian=TRUE)
  f.par[,i] <- maxlik$par
}

gama <- seq(from=1,to=80)
plot(gama, f.par[1,])
plot(gama, f.par[2,])
plot(gama, f.par[3,])

#need to show the Ni is a weibull dist with shape epsilon and scale alpha(s-gama)^delta,
#then can compute likelihood as in Q1.