rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#note now theta = (beta0, beta1, log(sigma), log(sigmab))
X <- model.matrix(~ rx + status, rats) #define X matrix
rats$litter <- factor(rats$litter)
Z <- model.matrix(~ litter - 1, rats) #define Z matrix
b <- rnorm(50,0,0.1) #set arbitrary b
y <- rats$time #define y vector
theta <- c(4.9831505, -0.2384417, -1.3324342, 0) #using weibullapprox from Q1, & small log(sigmab)

lfyb <- function (b, y, theta, X, Z) {
  #computes joint log density of y and b
  #first compute conditional density of y given b
  status <- X[,3]
  beta <- c(theta[1:2],0) #define beta vector
  eta <- as.numeric(X%*%beta + Z%*%b) #define eta vector
  loglik <- rep(0,150) #initialise loglik
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(eta) #scale vector
  
  #log conditional likelihood of y given b
  for (i in 1:length(y)) {
    if (status[i] == 0){
      loglik[i] <- (y[i]/lambda[i])^k
    }
    else {
      loglik[i] <- -log((dweibull(y[i], shape = k, scale = lambda[i])))
    }
  }
  lfy_b <- sum(loglik)
  
  #log likelihood of b
  sig.b <- exp(theta[4])
  lfb <- -sum(dnorm(b, 0, sig.b, log=TRUE))
  
  #joint density is (log likelihood of y given b) + (log likelihood of b)
  lf <- lfy_b + lfb
  #output lf
  lf
}

Grad1 <- function (theta, y, X, Z, b) {
  k <- 1/exp(theta[3]) #shape
  beta <- c(theta[1:2],0)
  b.aux <- Z%*%b
  eta <- as.numeric(X%*%beta + b.aux)
  y.vec <- (y^k)*k*exp(-k*eta) - b.aux/(3*exp(2*theta[4]))
  y.aux <- t(matrix(c(y.vec, rep(-k,150)), nrow=150, ncol=2))
  Grad <- diag(X[,c(1,3)]%*%y.aux)
  Grad <- matrix(Grad, nrow=3, ncol=50)
  Grad <- colSums(Grad)
  Grad
}

lal <- function (theta, y, X, Z) {
  #compute the negative log likelihood of theta using the laplace approximation
  #b <- rep(0,50)
  ratz <- optim(b, fn=lfyb, gr=Grad1, theta=theta, y=y, X=X, Z=Z, method="BFGS", hessian=TRUE)
  bhat <- ratz$par
  hess <- ratz$hessian
  lapprox <- -(length(bhat)/2)*log(2*pi) + (1/2)*sum(log(diag(hess))) + lfyb(b=bhat, y, theta, X, Z)
  lapprox
}

thetahat <- optim(theta, fn=lal, y=y, X=X, Z=Z, bhat=bhat, hess=hess, method="BFGS", hessian=TRUE)
se <- sqrt(diag(solve(thetahat$hessian)))
