rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#note now theta = (beta0, beta1, log(sigma), log(sigmab))
X <- model.matrix(~ rx, rats) #define X matrix

rats$litter <- factor(rats$litter)
Z <- model.matrix(~ litter - 1, rats) #define Z matrix

b <- rnorm(50,0,0.3950646) #set arbitrary b
y <- rats$time #define y vector
theta <- c(4.9831505, -0.2384417, -1.3324342, 0) #using weibullapprox from Q1, & small log(sigmab)

lfyb <- function (b, y, theta, X, Z) {
  #computes joint log density of y and b
  #first compute conditional density of y given b
  
  beta <- c(theta[1:2]) #define beta vector
  eta <- as.numeric(X%*%beta + Z%*%b) #define eta vector
  loglik <- 0 #initialise loglik
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(eta) #scale vector
  
  #log conditional likelihood of y given b
  for (i in 1:length(y)) {
    if (rats$status[i] == 0){
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

ratz <- optim(b, lfyb, y=y,X=X,theta=theta,Z=Z, method="BFGS")

lal <- function (theta, y, X, Z, ratz) {
  #compute the negative log likelihood of theta using the laplace approximation
  sig.b <- exp(theta[4])
  beta <- c(theta[1:2])
  eta <- X%*%beta
  bhat <- ratz$par
  
  hess <- matrix(0, nrow=50, ncol=50) #make hess
  for (j in 1:length(y)) {
    hess[floor((2.5+j)/3),floor((2.5+j)/3)] <- hess[floor((2.5+j)/3)] + (k^2)*(y[j]^k)*exp(-k*eta[j])
  }
  
  lapprox <- -(length(bhat)/2)*log(2*pi) + (1/2)*log(prod(diag(hess))) + lfyb(b=bhat, y, theta, X, Z)
  -lapprox
}

Grad1 <- function (y, rats, k, eta) {
  Grad <- matrix(0, nrow=50,ncol=1)
  for (j in 1:length(y)) {
    if (rats$status[j] == 1) {
      Grad[floor((2.5+j)/3)] <- Grad[floor((2.5+j)/3)] + (y[j]^k)*k*exp(-k*eta[j]) - k
    } else {
      Grad[floor((2.5+j)/3)] <- Grad[floor((2.5+j)/3)] + (y[j]^k)*k*exp(-k*eta[j])
    }
  }
}

#TOO MANY NANS
thetahat <- optim(par=theta, fn=lal, gr=Grad1, y=y,X=X,Z=Z,ratz=ratz, method="BFGS",hessian=TRUE)
