#read in the data
rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#define negative log likelihood function to minimise
nll.llogistic <- function (theta) {
  #initialise variables from dataset
  t <- rats$time
  s <- rats$status
  rx <- rats$rx
  #define shape and scale of log-logistic distribution
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(theta[1] + theta[2]*rx) #scale
  #compute negative log likelihood
  loglik <- -sum(log(((k/lambda)*(t/lambda)^(k-1)/(1+(t/lambda)^k)^2)^s*(1/(1+(t/lambda)^k))^(1-s)))
  loglik
}

#to pick initial point, vary paramaters of theta individually and fix the others to find an approximate minimum

#initialize variables
theta <- c(0,0,0)
beta0 <- theta[1]
beta1 <- theta[2]
logsig <- theta[3]
beta0var <- 0
beta1var <- 0
logsigvar <- 0

for (i in 0:1000) {
  theta[1] <- 2.5 + i*0.005
  beta0[i+1] <- theta[1]
  beta0var[i+1] <- nll.llogistic(theta)
}
plot(beta0, beta0var, type="l") #estimate optimal beta0 to be ~ 5.75

theta <- c(5.75,0,0)
for (i in 0:1000) {
  theta[2] <- -2.5 + i*0.005
  beta1[i+1] <- theta[2]
  beta1var[i+1] <- nll.llogistic(theta)
}
plot(beta1, beta1var, type="l") #estimate optimal beta1 to be ~ -0.5

theta <- c(5.75,-0.5,0)
for (i in 0:1000) {
  theta[3] <- -2.5 + i*0.005
  logsig[i+1] <- theta[3]
  logsigvar[i+1] <- nll.llogistic(theta)
}
plot(logsig, logsigvar, type="l") #estimate optimal logsigma to be ~ -0.5

theta <- c(5.75,-0.5,-0.5)
nll.optim <- optim(theta, nll.llogistic, hessian=TRUE, method="BFGS")
theta.llog <- nll.optim$par #optimal theta
hess <- nll.optim$hessian

#then to find the standard error, take the sqrt of the diag of the inverse hessian
se <- sqrt(diag(solve(hess)))

#95% CI = parameter estimate +- ~1.96*stderr
CI <- c(theta.llog[2] + qnorm(0.025)*se[2], theta.llog[2] + qnorm(0.975)*se[2])

#appears to show that tumours appeared sooner in rats which received treatment
