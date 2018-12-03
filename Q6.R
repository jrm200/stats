fatigue<- read.table("http://people.bath.ac.uk/kai21/ASI/fatigue.txt")

gama <- 80
theta <- c(1,1,1)

nll <- function (theta) {
  s <- fatigue$s
  N <- fatigue$N
  eps <- lgamma(1+theta[3])
  n <- length(s)
  loglik <- -n*log(exp(theta[1])) - sum(theta[2]*log(s-gama)) - n*eps
  logN <- sum(log(N))
  Nlik <- abs(logN + loglik)
  Nlik
}

loglik <- optim(theta, nll, method="BFGS", hessian=TRUE)

hess <- solve(loglik$hessian)
se <- sqrt(diag(hess))
