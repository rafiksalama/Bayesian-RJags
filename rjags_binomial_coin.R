thetacoins = rbeta(10,100,100)
m = matrix(ncol=10,nrow=100,0)
for(i in 1:10){
  m[,i] = rbinom(100,1,prob = thetacoins[i])
}

#Now let us do this with rjags
model_string <- "
model{
  #Start with the data
  #column
  for(j in 1:length(Y[1,]))
  {
    #row
    for(i in 1:length(Y[,1]))  
    {
      Y[i,j] ~ dbern(theta[j])
    }
  }
  for(j in 1:length(Y[1,])){
    theta[j] ~ dbeta(a,b)
  }
  a ~ dunif(1,1000)
  b ~ dunif(1,1000)
}
"

model <- jags.model(textConnection(model_string), data = list(Y=m))

update(model, 10000, progress.bar="none")
samp <- coda.samples(model, variable.names=c("theta","a","b"), n.iter=20000, progress.bar="none")
summary(samp)
