##Newton loop
b <- rnorm(50,0,0.01) #initial point
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
  
  g <- grad(theta,y,X,Z,b)
  
  if (max(abs(g)) <= (abs(ly)+0.001)*1e-10) { 
    break
  } else {
  #step <- -(1/Hess)%*%Grad
  step <- -solve(Hess,g)
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
