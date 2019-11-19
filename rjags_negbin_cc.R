library(rjags)
data = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O/CapC-786O_2/CapC-786O_2_REdig_sorted_CC1_chr10:99185458-99186090_PGAM1.gff")
chr=names(sort(table(data[,1]),decreasing = T)[1])
data = round(data[data[,1]==chr,6]/((data[data[,1]==chr,5]-data[data[,1]==chr,4])/1000))+1
max = which(data==max(data))
X = (c(1:length(data))-max)
window = min(min(max(X[X>0]),max(abs(X[X<0])))-1,500)
range = c((max-window):(max+window))
X = X[range]
data = data[range]
pdata = data
Y = pdata

#Now let us do this with rjags
model_string <- "model{

#Start with the data

  for(i in 1:length(Y)) 
  {
    mu[i] <- ((B[1]+(B[2]*exp(-B[3]*(i-c)^2))+(K*exp(-B[4]*(i-c)^2))))+(odds[i]*mu2[i])
    odds[i] = pi[i]/(1-pi[i])
    pi[i] ~ dbeta(.4,.4)
    p[i] = r[i]/(mu[i]+r[i])
    r[i] ~ dgamma(ar,br)
    Y[i] ~ dnegbin(p[i],r[i])
    mu2[i] ~ dgamma(am,bm)
  }
  
  for(j in 1:4){
    B[j] ~ dgamma(1,1)
  }
  
  #p2 ~ dbeta(1,1)
  #r2 ~ dgamma(am,bm)
  ar ~ dgamma(1,1)
  br ~ dgamma(1,1)
  am ~ dgamma(1,1)
  bm ~ dgamma(1,1)
}
"

model <- jags.model(textConnection(model_string), data = list(Y=as.vector(data),c=mid,K=max(data)))

update(model, 10000, progress.bar="none")
samp <- coda.samples(model, variable.names=c("br","ar","B","am","bm","pi","mu2"), n.iter=20000, progress.bar="none")
summary(samp)

#restrict length
pi1 = samp1[[1]][,grep("pi",colnames(samp1[[1]]))]
pi2 = samp[[1]][,grep("pi",colnames(samp[[1]]))]
#pi1 is larger than pi2, so we need to trim
pi1_start = ((ncol(pi1)-1001)/2)+1
pi2_start = ((ncol(pi2)-1001)/2)+1
pi1 = pi1[,pi1_start:(pi1_start+1000)]
pi2 = pi2[,pi2_start:(pi2_start+1000)]

#Now we need to multiply the two probability distributions
fitbeta <- function(x){
  m = mean(x)
  v = var(x)
  f = v+m^2-m
  A = -f*m/v
  B = f*(m-1)/v
  c(A,B)
}

v = vector()
for(i in 1:ncol(pi1))
{
  par1 = fitbeta(pi1[,i])
  par2 = fitbeta(pi2[,i])
  par = (par1+par2)/2
  v[i] = median(rbeta(1000000,par[1],par[2]))
  print(i)
}


pi1_1 = as.vector(apply(pi1,2,median))

pi2_2 = as.vector(apply(pi2,2,median))
odds_2 = pi2_2/(1-pi2_2)

