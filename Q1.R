ratlikelihood <- function (theta) {
  #create subset of data since for the purposes of calculating the optimal theta since when status = 0,
  #the survival function does not depend on theta. Need to optimise this function using optim.
  t <- rats[which(rats$status==1),,]$time
  
  #compute the negative log likelihood
  loglik <- -sum(log(dweibull(t, shape=1/exp(theta[3]), scale=exp(theta[1] + theta[2]))))
  loglik
}
