log.posterior <- function (theta, times = y, X.mat=X, Z.mat=Z) {
  #log likelihood
  prior <- dexp(x=exp(theta[4]), rate=5, log=TRUE)
  prior[prior==-Inf] <- -100
  l.lik <- -lal(theta, times, X.mat, Z.mat) + prior
  l.lik
}

n.rep <- 10000
sigma.prop <- 0.1275 #began with 0.05, ar = 0.538. Tuned to 0.1275(ish) by varying sigma until ar ~ 0.23

MH <- function (theta0, sigma.prop, n.rep) {
  theta.vals <- matrix(0, n.rep, 4) #matrix to save generated values
  theta.vals[1,] <- theta0
  lp0 <- log.posterior(theta.vals[1,])
  alpha <- rep(0, n.rep)
  
  for (i in 2:n.rep) {
    current_theta <- theta.vals[i-1,]
    proposed_theta <- current_theta + rnorm(4, sd=sigma.prop)
    lp1 <- log.posterior(proposed_theta)
    acc <- exp(min(0, lp1-lp0))
    test.value <- runif(1)
    if (test.value <= acc) {
      theta.vals[i,] <- proposed_theta
      lp0 <- lp1
      alpha[i] <- 1
    } else {
      theta.vals[i,] <- current_theta
      lp1 <- lp0
    }
  }
  accept.rate <- sum(alpha)/n.rep
  list(theta=theta.vals, accept.rate=accept.rate, alpha=alpha)
}

theta0 <- c(1,0,0,0)
mh <- MH(theta0, sigma.prop, n.rep)
theta <- mh$theta
accept.rate <- mh$accept.rate

acf(theta[,4], lag.max=1000) #looks like the 'stickiest' parameter is theta[4], acf = 0 after about 200 iterations, effective
acf(theta[,3], lag.max=200) #sample size = n.rep/200
acf(theta[,2], lag.max=200)
acf(theta[,1], lag.max=200)

theta <- theta[1001:n.rep,]

par(mfrow=c(3,2),mar=c(4,4,1,1))
plot(theta[,1],theta[,2],xlab="B0",ylab="B1",pch=20,cex=1)
plot(theta[,1],theta[,3],xlab="B0",ylab="log(sig)",pch=20,cex=1)
plot(theta[,1],theta[,4],xlab="B0",ylab="log(sig.b)",pch=20,cex=1)
plot(theta[,2],theta[,3],xlab="B1",ylab="log(sig)",pch=20,cex=1)
plot(theta[,2],theta[,4],xlab="B1",ylab="log(sig.b)",pch=20,cex=1)
plot(theta[,3],theta[,4],xlab="log(sig)",ylab="log(sig.b)",pch=20,cex=1)
