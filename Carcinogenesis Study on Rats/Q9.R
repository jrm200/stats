fatigue<- read.table("http://people.bath.ac.uk/kai21/ASI/fatigue.txt")

set.seed(7)

s <- fatigue$s
N <- fatigue$N
ro <- fatigue$ro
gama <- runif(26, min=0.01, max=s) #define random effect, gama between 0 and the value of s which it relates to

log.post <- function (theta, gama, s.=s, ro.=ro, N.=N) {
  #log density of N|gama
  k <- 1/exp(theta[3]) #shape
  lambda <- exp(theta[1])*(s. - gama)^theta[2] #scale
  ll0 <- sum(dweibull(N., shape=k, scale=lambda, log=TRUE)*(1-ro.))
  ll1 <- sum((-(N./lambda)^k)*ro.)
  ll <- ll0 + ll1
  if (is.nan(ll) || ll == -Inf) {
    ll <- -1e3
  }
  #log density of gama
  ll.gama <- sum(dweibull(gama, shape = 1/exp(theta[5]), scale = exp(theta[4]), log=TRUE))
  if (is.nan(ll.gama) || ll.gama == -Inf) {
    ll.gama <- -1e3
  }
  #log density of prior of sig.gama
  ll.prior <- dexp(exp(theta[5]), rate=5, log=TRUE)
  if (ll.prior == -Inf) {
    ll.prior <- -1e3
  }
  #sum of log densities
  lp <- ll + ll.gama + ll.prior
  lp
}


MH <- function (theta, sigma.prop, n.rep, s.=s, ro.=ro, N.=N) {
  theta.vals <- matrix(0, n.rep, 5) #matrix to save generated values
  theta.vals[1,] <- theta
  b <- runif(26, min=0.01, max=s.)
  b.vals <- matrix(0, n.rep, 26)
  b.vals[1,] <- b
  lp0 <- log.post(theta.vals[1,], b)
  accept.th <- 0
  accept.b <- 0
  
  for (i in 2:n.rep) {
    #update theta
    theta <- theta + rnorm(5, 0, sigma.prop[1:5])
    lp1 <- log.post(theta, b)
    if (runif(1) < exp(lp1 - lp0)){
      accept.th <- accept.th+1
      lp0 <- lp1
    }else{
      theta <- theta.vals[i-1,]
    }
    #update random effects
    b.step <- sigma.prop[6]
    b <- b + runif(26, -s., s.)*b.step #maybe tune
    while (length(b[b<0 | b>s.]) > 0) {
      b[b<0 | b>s.] <- b.vals[i-1,which(b<0 | b>s.)] + runif(length(b[b<0 | b>s.]), -s., s.)*b.step 
    }
    
    lp1 <- log.post(theta, b)
    if (runif(1) < exp(lp1 - lp0)){
      accept.b <- accept.b+1
      lp0 <- lp1
    }else{
      b <- b.vals[i-1,]
    }
    theta.vals[i,] <- theta
    b.vals[i,] <- b
  }
  accept.rate <- c(accept.th/n.rep, accept.b/n.rep)
  list(theta=theta.vals, accept.rate=accept.rate, gama=b.vals)
}

theta <- c(5, -2, 0, 4, -1.5)
sigma.prop <- c(0.055, 0.05, 0.05, 0.05, 0.05, 0.1)
n.rep <- 100000
mh <- MH(theta, sigma.prop, n.rep)
ar <- mh$accept.rate
print(ar)
lower <- n.rep/5+1

par(mfrow=c(3,2),mar=c(4,4,1,1))
for (i in 1:5){
  plot(mh$theta[lower:n.rep,i], type="l")
}
for (i in 1:7){
  plot(mh$gama[lower:n.rep,3*i], type="l")
}
