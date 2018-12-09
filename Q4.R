rats<- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#theta = [beta0, beta1, log(sigma), log(sigmab)]
X <- model.matrix(~ rx + status, rats) #define X matrix
rats$litter <- factor(rats$litter)
Z <- model.matrix(~ litter - 1, rats) #define Z matrix
b <- rnorm(50,0,0.1) #set arbitrary b
y <- rats$time #define y vector
theta <- c(4.9831505, -0.2384417, -1.3324342, 0) #using optimal paramater estimates from Q1, & small log(sig.b)

lfyb <- function (b, y, theta, X, Z) {
  #function to compute the joint log density of y and b, which is the
  #(conditional log likelihood of y given b) + (likelihood of b)
  #first compute conditional density of y given b
  theta[theta<(-3)] <- -3 #constraints on theta
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

#status = 1
expr1 <- expression(k*(B0 + B1*rx + b.aux) + (y/exp(B0 + B1*rx + b.aux))^k)
nll.aux1 <- deriv(expr1, c("b.aux"), function.arg=c("b.aux", "B0", "B1", "k", "rx", "y"))

#status = 0
expr0 <- expression((y/exp(B0 + B1*rx + b.aux))^k)
nll.aux0 <- deriv(expr0, c("b.aux"), function.arg=c("b.aux", "B0", "B1", "k", "rx", "y"))

grad <- function (b, y, theta, X, Z) {
  #define variables
  theta[theta<(-3)] <- -3 #constraints on theta
  rx <- X[,2]
  s <- X[,3]
  B0 <- theta[1]
  B1 <- theta[2]
  k <- 1/exp(theta[3])
  b.aux <- Z%*%b
  #gradient if status = 1
  b.grad1 <- nll.aux1(b.aux, B0, B1, k, rx, y)
  g1 <- attr(b.grad1, "gradient")
  g1 <- as.numeric(g1)
  #gradient if status = 0
  b.grad0 <- nll.aux0(b.aux, B0, B1, k, rx, y)
  g0 <- attr(b.grad0, "gradient")
  g0 <- as.numeric(g0)
  #select contribution to gradient based on status
  g <- s*g1 +(1-s)*g0
  #gradient of negative log likelihood of b
  b.grad <- b/exp(2*theta[4])
  m <- matrix(g, nrow=3, ncol=50)
  g <- colSums(m) #add contributions from b terms
  g <- (g+b.grad)
  g
}

bhat <- optim(b, fn=lfyb, gr=grad, theta=theta, y=y, X=X, Z=Z, method="BFGS", hessian=TRUE)

lal <- function (theta, y, X, Z) {
  #constraints on theta
  #compute the negative log likelihood of theta using the laplace approximation
  b <- rep(0,50)
  lfyb.opt <- optim(b, fn=lfyb, gr=grad, theta=theta, y=y, X=X, Z=Z, method="BFGS", hessian=TRUE)
  b.opt <- lfyb.opt$par
  hess <- lfyb.opt$hessian
  lapprox <- -(length(b.opt)/2)*log(2*pi) + (1/2)*log(det(hess)) + lfyb(b=b.opt, y=y, theta=theta, X=X, Z=Z)
  lapprox
}

lal.opt <- optim(theta, lal, y=y, X=X, Z=Z, method="BFGS", hessian=TRUE)
thetahat <- lal.opt$par
se <- sqrt(diag(solve(lal.opt$hessian)))
CI1 <- c(thetahat[1] + qnorm(0.025)*se[1], thetahat[1] + qnorm(0.975)*se[1])
CI.rand <- c(thetahat[2] + qnorm(0.025)*se[2], thetahat[2] + qnorm(0.975)*se[2])
CI3 <- c(thetahat[3] + qnorm(0.025)*se[3], thetahat[3] + qnorm(0.975)*se[3])
CI4 <- c(thetahat[4] + qnorm(0.025)*se[4], thetahat[4] + qnorm(0.975)*se[4])
CI <- matrix(c(CI1, CI.rand, CI3, CI4), nrow=2, ncol=4)
