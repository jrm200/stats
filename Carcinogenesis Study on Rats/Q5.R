lal <- function (theta, y, X, Z) {
  #compute the negative log likelihood of theta using the laplace approximation
  ratz <- optim(b, fn=lfyb, theta=theta, y=y, X=X, Z=Z, method="BFGS", hessian=TRUE)
  bhat <- ratz$par
  hess <- ratz$hessian
  lapprox <- -(length(bhat)/2)*log(2*pi) + (1/2)*sum(log(diag(hess))) + lfyb(b=bhat, y, theta, X, Z)
  #print(diag(hess))
  lapprox
}

log.posterior <- function (theta, times = y, X.mat=X, Z.mat=Z) {
  #log likelihood
  l.lik <- -lal(theta, times, X.mat, Z.mat) + dexp(x=theta[4], rate=5, log=TRUE)
  l.lik
}

n.rep <- 10
sigma.prop <- 0.5

MH <- function (theta0, sigma.prop, n.rep) {
 theta.vals <- matrix(0, n.rep, 4) #matrix to save generated values
 theta.vals[1,] <- theta0
 lp0 <- log.posterior(theta.vals[1,])
 alpha <- rep(0, n.rep)
 
 for (i in 2:n.rep) {
   current_theta <- theta.vals[i-1,]
   proposed_theta <- current_theta + rnorm(4, sd=sigma.prop)/10
   lp1 <- log.posterior(proposed_theta)
   acc <- exp(min(0, lp1-lp0))
   test.value <- runif(1)
   if (lp1 != -Inf && test.value <= acc) {
    theta.vals[i,] <- proposed_theta
    lp0 <- lp1
    alpha[i] <- 1
    } else {
    theta.vals[i,] <- current_theta
    lp1 <- lp0
   }
 }
 accept.rate <- sum(alpha)/n.rep
 list(theta=theta.vals, accept.rate=accept.rate)
}

mh <- MH(theta0 = thetahat$par, sigma.prop, n.rep)
theta <- mh$theta
accept.rate <- mh$accept.rate