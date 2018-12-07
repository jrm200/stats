rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#note now theta = (beta0, beta1, log(sigma), log(sigmab))
X <- model.matrix(~ rx, rats) #define X matrix

rats$litter <- factor(rats$litter)
Z <- model.matrix(~ litter - 1, rats) #define Z matrix

#b <- rnorm(50,0,0.01) #set arbitrary b
b <- rep(0,50)
y <- rats$time #define y vector
theta <- c(1,1,0,0)
#theta <- c(4.9831505, -0.2384417, -1.3324342, 0) #using weibullapprox from Q1, & small log(sigmab)

#Check result with other group?
lfyb <- function (b, y, theta, X, Z) {
  #computes joint log density of y and b
  #first compute conditional density of y given b
  
  beta <- c(theta[1:2]) #define beta vector
  eta <- as.numeric(X%*%beta + Z%*%b) #define eta vector
  loglik <- 0 #initialise loglik
  k <- 1/exp(theta[3])#shape
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

#hess always PSD!

##Newton loop
b <- rnorm(50,0,0.01) #initial point
#b <- rep(1,50)
k <- 1/exp(theta[3])
beta <- c(theta[1:2])
eta <- as.numeric(X%*%beta + Z%*%b)
ly <- lfyb(b, y, theta, X, Z)
iter.max <- 100 #set max iterations
f.vals <- rep(0, iter.max)
sig.b <- exp(theta[4])

for (i in 1:iter.max) {
  
  f.vals[i] <- ly
  Hess <- matrix(0, nrow=50, ncol=50) #make hess
  for (j in 1:length(y)) {
    Hess[floor((2.5+j)/3),floor((2.5+j)/3)] <- Hess[floor((2.5+j)/3)] + (k^2)*(y[j]^k)*exp(-k*eta[j]) + 1/(3*sig.b^2)
  }
  
  Grad <- Grad1(theta,y,X,Z,b)
  
  if (max(abs(Grad)) <= (abs(ly)+0.001)*1e-10) { 
    break
  } else {
  #step <- -(1/Hess)%*%Grad
  step <- -solve(Hess,Grad)
  step <- step
  ly1 <- lfyb((b+step), y, theta, X, Z)
  m <- 0
  while (ly1 > ly && m < 50) { #backtracking
    #print("step being halved")
    step <- step/2
    m <- m + 1
    ly1 <- lfyb(b+step, y, theta, X, Z)
  }
  ly <- ly1
  b <- b+step
  print(m)
}
}
bhat <- b

chat <- optim(b, fn=lfyb, gr=Grad1, theta=theta, y=y, X=X, Z=Z, method="BFGS", hessian=TRUE)





lal <- function (theta, y, X, Z, bhat) {
  #compute the negative log likelihood of theta using the laplace approximation
  
  ## Function doesn't work because optim is inefficient in 50 dimensions, use Newton loop or pass gradient to optim
  ## Also can calculate hessian algebraically to save time
  sig.b <- exp(theta[4])
  beta <- c(theta[1:2])
  eta <- X%*%beta
  lapprox <- (2*pi)^(length(bhat)/2)/sqrt(prod(diag(hess)))*lfyb(b=bhat, y, theta, X, Z)
  lapprox
}

#TOO MANY NANS
#This optim should be fine as is only in 4 dimensions. Should work once lal is resolved
thetahat <- optim(theta, lal, y=y,X=X,Z=Z,method="BFGS",hessian=TRUE)

