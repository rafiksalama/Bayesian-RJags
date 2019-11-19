library(rjags)
data = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O/CapC-786O/CapC-786O_REdig_sorted_CC1_chr16:30076905-30077296_ALDOA.gff")
data = round(data[data[,1]=="chr16",6]/((data[data[,1]=="chr16",5]-data[data[,1]=="chr16",4])/1000))+1
max = which(data==max(data))
X = abs(c(1:length(data))-max)
range = c((max-750):(max+750))
X = X[range]+1
data = data[range]
pdata = data
Y = pdata

#Now let us do this with rjags
model_string <- "model{

#Start with the data

for(i in 1:length(Y)) 
{
#     nbin[i] ~ dnegbin(p,r)
#     Y[i] <- w[i]*norm[i]+(1-w[i])*nbin[i]
#     w[i] ~ dbeta(1,1)
    mu[i] <- B[1]+(B[2]*exp(-B[3]*X[i]))+(max(X)*exp(-B[4]*X[i]^2))
    Y[i] ~ dnorm(mu[i],inv.var)
}

#  p ~ dbeta(1.001,1.001)
#  r ~ dgamma(0.01,0.01)

  for(j in 1:4) {
    B[j] ~ dgamma(1,1)
  }

  inv.var   ~ dgamma(0.01, 0.01)
  sigma     <- 1/sqrt(inv.var)
}
"

model <- jags.model(textConnection(model_string), data = list(Y=as.vector(Y),X=X))

update(model, 10000, progress.bar="none")
samp <- coda.samples(model, variable.names=c("B","inv.var"), n.iter=20000, progress.bar="none")
summary(samp)
