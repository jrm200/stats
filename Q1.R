ratlikelihood <- function (theta) {
  #set up variables of data to compute likelihood
  rats2 <- rats[which(rats$status==1),,]
  t <- rats2$time
  
  loglik <- sum(log(dweibull(t, shape=1/exp(theta[3]), scale=exp(theta[1] + theta[2]))))

  loglik
}