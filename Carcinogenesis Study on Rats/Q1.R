#read in the data
rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#define negative log likelihood function to minimise
nll.weibull <- function (theta) {
  #initialise variables from dataset
  s <- rats$status
  rx <- rats$rx
  t <- rats$time
  #define shape and scale of Weibull distribution
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(theta[1] + theta[2]*rx) #scale
  #compute negative log likelihood
  loglik <- -sum(log(((k/lambda)*(t/lambda)^(k-1)*exp(-(t/lambda)^k))^s*(exp(-(t/lambda)^k))^(1-s)))
  loglik
}

#to pick initial theta, vary paramaters of theta individually and fix the others to find an approximate minimum

#initialize variables
theta <- c(0,0,0)
beta0 <- theta[1]
beta1 <- theta[2]
logsig <- theta[3]
beta0var <- 0
beta1var <- 0
logsigvar <- 0

for (i in 0:1000) {
  theta[1] <- 4.5 + i*0.0025
  beta0[i+1] <- theta[1]
  beta0var[i+1] <- nll.weibull(theta)
}

plot(beta0, beta0var, type="l") #estimate optimal beta0 to be ~5.75

theta <- c(0,0,0)
for (i in 0:1000) {
  theta[2] <- 4 + i*0.0025
  beta1[i+1] <- theta[2]
  beta1var[i+1] <- nll.weibull(theta)
}

plot(beta1, beta1var, type="l") #estimate optimal beta1 to be ~5.25

theta <- c(0,0,0)
for (i in 0:1000) {
  theta[3] <- 1.5 + i*0.005
  logsig[i+1] <- theta[3]
  logsigvar[i+1] <- nll.weibull(theta)
}

plot(logsig, logsigvar, type="l") #estimate optimal logsigma to be ~2.75

theta <- c(5.75,5.25,2.75)
weibull.optim <- optim(theta, nll.weibull, method="BFGS", hessian=TRUE)
theta.weibull <- weibull.optim$par
hess <- solve(weibull.optim$hessian)

#then to find the standard error, take the square root of the diagonal of the inverse hessian
se <- sqrt(diag(hess))

#95% CI = parameter estimate +- ~1.96*stderr
CI.norand <- c(theta.weibull[2]+qnorm(0.025)*se[2], theta.weibull[2]+qnorm(0.975)*se[2])

#appears to show that tumours appeared sooner in rats which received treatment
