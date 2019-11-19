library(rjags)
X = cbind(rep(1,100),apply(as.matrix(rgamma(9,1)),1,function(x){rnorm(100,x)}))
B = as.matrix(c(10,rgamma(9,1)))
Y = X%*%B
Bhat = t((t(X)%*%Y))%*%solve(t(X)%*%X)
#Now let us do this with rjags
model_string <- "model{

#Start with the data

  for(i in 1:length(Y)) 
  {
    Y[i] ~ dnorm(mu[i],inv.var)
    mu[i] <- B[1]+X[i,1]*B[2]+X[i,2]*B[3]+X[i,3]*B[4]+X[i,4]*B[5]+X[i,5]*B[6]+X[i,6]*B[7]+X[i,7]*B[8]+X[i,8]*B[9]+X[i,9]*B[10]
  }

  for(j in 1:10) {
    B[j] ~ dgamma(1,1)
  }

  inv.var   ~ dgamma(0.01, 0.01)
  sigma     <- 1/sqrt(inv.var)
}
"

model <- jags.model(textConnection(model_string), data = list(Y=as.vector(Y),X=X))

update(model, 10000, progress.bar="none")
samp <- coda.samples(model, variable.names=c("B","sigma"), n.iter=100000, progress.bar="none")
summary(samp)
