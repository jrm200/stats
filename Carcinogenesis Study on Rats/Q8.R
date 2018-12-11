fatigue<- read.table("http://people.bath.ac.uk/kai21/ASI/fatigue.txt")

gama <- 33 #arbitrary gamma
theta <- c(1,1,1) #arbitrary initial theta
s <- fatigue$s #stress

nll.weibull <- function (theta, gama, s) {
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(theta[1])*(s - gama)^theta[2] #scale
  ro <- fatigue$ro
  #compute the negative log likelihood
  nll0 <- -sum(dweibull(s, shape=k, scale=lambda, log=TRUE)^(1-ro))
  nll1 <- sum(((s/lambda)^k)^ro)
  nll <- nll0 + nll1
  nll
}

maxlik <- optim(theta, nll.weibull, gama=gama, s=s, method="BFGS", hessian=TRUE)
se <- sqrt(diag(solve(maxlik$hessian)))
theta.opt <- maxlik$par

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
plot(gama, f.vals) #QUESTION 7

gama.opt <- 33
theta.opt <- c(theta.opt, gama.opt)
s <- fatigue$s
N <- fatigue$N
s <- seq(from=34, to=150)
q <- qweibull(0.1, shape = 1/exp(theta.opt[3]), scale = exp(theta.opt[1])*(s-theta.opt[4])^theta.opt[2])
plot(s,log(N))

#need to show the Ni is a weibull dist with shape epsilon and scale alpha(s-gama)^delta,
#then can compute likelihood as in Q1.
