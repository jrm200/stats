sigma <- 1

b <- rnorm(50, sd = sigma)

X <- model.matrix(~ rx, rats)

rats$litter <- factor(rats$litter)

Z <- model.matrix(~ litter - 1, rats)

y <- rats$time

lfyb <- function (theta, y, b, X, Z) {
  beta <- c(theta[1], theta[2])
  eta <- as.numeric(X%*%beta + Z%*%b)
  loglik <- 0
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(eta)
  
  for (i in 1:150){
    if (rats$status[i] == 0){
      loglik <- loglik + (y[i]/lambda[i])^k
    }
    else{
      loglik <- loglik - log(dweibull(y[i], shape = k, scale = lambda[i]))
    }
  }
  
  loglik
}
