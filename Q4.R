rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#note now theta = (beta0, beta1, log(sigma), log(sigmab))
X <- model.matrix(~ rx, rats) #define X matrix

rats$litter <- factor(rats$litter)
Z <- model.matrix(~ litter - 1, rats) #define Z matrix

b <- rep(0, ncol(Z)) #set arbitrary b
y <- rats$time #define y vector
theta <- c(0,0,0,0) #arbitrary

lfyb <- function (theta, y, b, X, Z) {
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
      loglik[i] <- -(dweibull(y[i], shape = k, scale = lambda[i], log=TRUE))
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

#STARTED BUT NEEDS WORK
lal <- function (theta, y, X, Z) {
  #compute the negative log likelihood of theta using the laplace approximation
  sig.b <- exp(theta[4])
  beta <- c(theta[1:2])
  eta <- X%*%beta
  b <- rep(0,ncol(Z))
  
  lf <- lfyb(theta, y, b, X, Z)
  optim(par=b, fn=lfyb)
}
