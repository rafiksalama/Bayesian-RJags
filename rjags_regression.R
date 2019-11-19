library('coda')
library('rjags') # observe startup messages

#The model will start by setting the distribution of the data observed
y = myData$y
y = myData$s

Ntotal = length(y)
dataList = list( y=y , Ntotal=Ntotal)

modelString = "model{ 
  for ( i in 1: Ntotal) {
    y[i] ~ dbern( theta[s[i]] )
   }
  for ( s in 1 : Nsubj ) {
    theta[s] ~ dbeta( omega * (kappa-2)+1, (1-omega )*(kappa-2)+1)
  }
  omega ~ dbeta(1,1)
  kappa <- kappaMinusTwo + 2
  kappaMinusTwo ~ dgamma( 0.01, 0.01 )
}"
# mean = 0.01542, sd = 0.1232165 (generic vague ) 
writeLines( modelString, con='TEMPmodel.txt' )
