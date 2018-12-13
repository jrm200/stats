fatigue<- read.table("http://people.bath.ac.uk/kai21/ASI/fatigue.txt")
#need to show the Ni is a weibull dist with shape epsilon and scale alpha(s-gama)^delta,
#then can compute likelihood as in Q1.
#Question 2.1
gama <- 50 #arbitrary gamma
theta <- c(1,1,1) #arbitrary initial theta
s <- fatigue$s #stress
N <- fatigue$N

nll.weibull <- function (theta, gama, s, N) {
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(theta[1])*(s - gama)^theta[2] #scale
  ro <- fatigue$ro
  #compute the negative log likelihood
  nll0 <- -sum(dweibull(N, shape=k, scale=lambda, log=TRUE)^(1-ro))
  nll1 <- sum(((N/lambda)^k)^ro)
  nll <- nll0 + nll1
  nll
}

maxlik <- optim(theta, nll.weibull, gama=gama, s=s, N=N, method="BFGS", hessian=TRUE)
se <- sqrt(diag(solve(maxlik$hessian)))
theta.opt <- maxlik$par

CI1 <- c(maxlik$par[1]+qnorm(0.025)*se[1], maxlik$par[1]+qnorm(0.975)*se[1])
CI2 <- c(maxlik$par[2]+qnorm(0.025)*se[2], maxlik$par[2]+qnorm(0.975)*se[2])
CI3 <- c(maxlik$par[3]+qnorm(0.025)*se[3], maxlik$par[3]+qnorm(0.975)*se[3])
CI <- matrix(c(CI1,CI2,CI3), ncol=3, nrow=2)

#How variations in gama affects the theta
f.vals <- rep(0,80)
f.par <- matrix(0,nrow=3, ncol=80)
for (i in 1:80) {
  gama <- i
  f.vals[i] <- nll.weibull(theta, gama, s, N)
  maxlik <- optim(theta, nll.weibull, gama=gama, s=s, N=N, hessian=TRUE)
  f.par[,i] <- maxlik$par
}

gama <- seq(from=1,to=80)
plot(gama, f.par[1,])
plot(gama, f.par[2,])
plot(gama, f.par[3,])

#Q7 we can estimate gama by including it within theta and passing it to optim
#redefine nll.weibull with theta[4] = gama
gama <- 50 #reset gama
theta.opt <- c(theta.opt, gama)

nll.weibull2 <- function (theta, s, N) {
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(theta[1])*(s - theta[4])^theta[2] #scale
  ro <- fatigue$ro
  #compute the negative log likelihood
  nll0 <- -sum(dweibull(N, shape=k, scale=lambda, log=TRUE)^(1-ro))
  nll1 <- sum(((N/lambda)^k)^ro)
  nll <- nll0 + nll1
  nll
}

maxlik2 <- optim(theta.opt, nll.weibull2, s=s, N=N, hessian=TRUE)
theta.opt2 <- maxlik2$par
gama.opt <- theta.opt2[4] #gama looks like it should be around 66

#Question 8
#function to compute quantile of N
N.quantile <- function (theta, s, quantile) {
  ZQ <- qweibull(quantile, shape=1, scale=1)
  NQ <- exp(theta[1])*(s-theta[4])^theta[2]
  Q <- NQ*ZQ^exp(theta[3])
  Q
}

plot(s, log(N))
S <- seq(from=min(s), to=max(s), by=0.1)
N.lower.quantile <- N.quantile(theta=theta.opt2, s=s, quantile=0.1)
lines(s, log(N.lower.quantile), col="red")
N.middle.quantile <- N.quantile(theta=theta.opt2, s=s, quantile=0.5)
lines(s, log(N.middle.quantile), col="green")
