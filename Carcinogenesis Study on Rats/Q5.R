rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

set.seed(4)

#theta = [beta0, beta1, log(sigma), log(sigmab)]
X <- model.matrix(~ rx + status, rats) #define X matrix
rats$litter <- factor(rats$litter)
Z <- model.matrix(~ litter - 1, rats) #define Z matrix
b <- rnorm(50,0,0.1) #set arbitrary b
y <- rats$time #define y vector
theta <- c(4.9831505, -0.2384417, -1.3324342, 0) #using optimal paramater estimates from Q1, & small log(sig.b)


log.posterior <- function (b, theta, times = y, X.mat=X, Z.mat=Z) {
  #log likelihood
  prior <- dexp(x=exp(theta[4]), rate=5, log=TRUE)
  prior[prior==-Inf] <- -10
  l.lik <- -lfyb(b, times, theta, X.mat, Z.mat) + prior
  l.lik
}

lfyb <- function (b, y, theta, X, Z) {
  #function to compute the joint log density of y and b, which is the
  #(conditional log likelihood of y given b) + (likelihood of b)
  #first compute conditional density of y given b
  s <- X[,3] #status
  beta <- c(theta[1:2],0) #define beta vector with a 0 in the last position since status is not relevant for eta
  eta <- as.numeric(X%*%beta + Z%*%b) #define eta vector
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(eta) #scale vector
  #log conditional likelihood of y given b
  nll1 <- k*log(lambda) - log(k) - (k-1)*log(y) + (y/lambda)^k #scale = 1
  nll0 <- (y/lambda)^k #scale = 0
  lfy_b <- sum(s*nll1 + (1-s)*nll0)
  #negative log likelihood of b
  sig.b <- exp(theta[4])
  #lfb <- -sum(dnorm(b, 0, sig.b, log=TRUE))
  lfb <- dnorm(b, 0, sig.b)
  lfb <- log(lfb)
  lfb[lfb==-Inf] <- -10
  lfb <- -sum(lfb)
  #compute and output negative joint log density
  lf <- lfy_b + lfb
  lf
}

n.rep <- 100000
sigma.prop <- c(0.15, 0.055, 0.055, 0.1) #tuning possible; (0.15, 0.055, 0.055, 0.1) looks to be very good in terms of acceptance rate

MH <- function (theta, sigma.prop, n.rep) {
  theta.vals <- matrix(0, n.rep, 4) #matrix to save generated values
  theta.vals[1,] <- theta
  b <- rep(0, 50)
  b.vals <- matrix(0, n.rep, 50)
  b.vals[1,] <- b
  lp0 <- log.posterior(b, theta.vals[1,])
  accept.th <- 0
  accept.b <- 0

  for (i in 2:n.rep) {
    #update theta
    theta <- theta+rnorm(4, 0, sigma.prop)
    lp1 = log.posterior(b, theta)
    if (runif(1) < exp(lp1 - lp0)){
      accept.th <- accept.th+1
      lp0 <- lp1
    }else{
      theta <- theta.vals[i-1,]
    }
    #update random effects
    b <- b + rnorm(50)*0.045 #tuning possible; 0.045 seems to be good
    lp1 <- log.posterior(b, theta)
    if (runif(1) < exp(lp1 - lp0)){
      accept.b <- accept.b+1
      lp0 <- lp1
    }else{
      b <- b.vals[i-1,]
    }
    theta.vals[i,] <- theta
    b.vals[i,] <- b
  }
  accept.rate <- c(accept.th/n.rep, accept.b/n.rep)
  list(theta=theta.vals, accept.rate=accept.rate)
}

theta0 <- c(1,0,0.1,0.2)
mh <- MH(theta0, sigma.prop, n.rep)
theta <- mh$theta
accept.rate <- mh$accept.rate

par(mfrow=c(2,2),mar=c(4,4,1,1))
acf(theta[,4], lag.max=2000) #looks like the 'stickiest' parameter is theta[4], acf = 0 after about 200 iterations, effective
acf(theta[,3], lag.max=2000) #sample size = n.rep/200
acf(theta[,2], lag.max=2000)
acf(theta[,1], lag.max=2000)

lower <- n.rep/5 + 1
theta <- theta[lower:n.rep,]

par(mfrow=c(3,2),mar=c(4,4,1,1))
plot(theta[,1],theta[,2],xlab="B0",ylab="B1",pch=20,cex=1)
plot(theta[,1],theta[,3],xlab="B0",ylab="log(sig)",pch=20,cex=1)
plot(theta[,1],theta[,4],xlab="B0",ylab="log(sig.b)",pch=20,cex=1)
plot(theta[,2],theta[,3],xlab="B1",ylab="log(sig)",pch=20,cex=1)
plot(theta[,2],theta[,4],xlab="B1",ylab="log(sig.b)",pch=20,cex=1)
plot(theta[,3],theta[,4],xlab="log(sig)",ylab="log(sig.b)",pch=20,cex=1)

par(mfrow=c(2,2),mar=c(4,4,1,1))
for (i in 1:4){
  plot(theta[,i])
}
