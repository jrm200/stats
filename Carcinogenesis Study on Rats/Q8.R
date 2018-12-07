gama <- 80
theta <- maxlik$par
s <- fatigue$s

N <- exp(theta[1])*(s-gama)^theta[2]*qweibull(0.1,shape=1)^exp(theta[3])

plot(s,N)

x <- seq(from=0,to=200, by= 0.01)
y <- dweibull(x, shape = 1, scale = exp(theta[1])*(81 - gama)^theta[2])
plot(x,y)