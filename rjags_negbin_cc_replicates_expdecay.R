library(rjags)
data1 = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O/CapC-786O/CapC-786O_REdig_sorted_CC1_chr1:8938714-8938951_ENO1.gff")
chr=names(sort(table(data1[,1]),decreasing = T)[1])
data1 = round(data1[data1[,1]==chr,6]/((data1[data1[,1]==chr,5]-data1[data1[,1]==chr,4])/1000))+1
max = which(data1==max(data1))
X = (c(1:length(data1))-max)
window = min(min(max(X[X>0]),max(abs(X[X<0])))-1,500)
range = c((max-window):(max+window))
X = X[range]
data1 = data1[range]

data2 = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O/CapC-786O_2/CapC-786O_2_REdig_sorted_CC1_chr1:8938714-8938951_ENO1.gff")
chr=names(sort(table(data2[,1]),decreasing = T)[1])
data2 = round(data2[data2[,1]==chr,6]/((data2[data2[,1]==chr,5]-data2[data2[,1]==chr,4])/1000))+1
max = which(data2==max(data2))
X = (c(1:length(data2))-max)
window = min(min(max(X[X>0]),max(abs(X[X<0])))-1,500)
range = c((max-window):(max+window))
X = X[range]
data2 = data2[range]

Y = cbind(data1,data2)
plot(X,Y[,1], type="l", ylim=c(0,500))
lines(X,Y[,2], type="l", col="red")

#Now let us do this with rjags
model_string <- "model{

#Start with the data

  for(i in 1:length(Y[,1])) 
  {
    mu[i] ~ dpois((B[1]+(K/exp(sqrt((i-c)^2)/B[2])))+(odds[i]*mu2[i]))
    odds[i] = pi[i]/(1-pi[i])
    pi[i] ~ dbeta(.5,.5)
    p[i] = r[i]/(mu[i]+r[i])
    r[i] ~ dgamma(ar,br)
    mu2[i] ~ dgamma(am,bm)
    Y[i,1] ~ dnegbin(p[i],r[i])
    Y[i,2] ~ dnegbin(p[i],r[i])
  }
  
  for(j in 1:2){
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

model <- jags.model(textConnection(model_string), data = list(Y=Y,c=window,K=max(Y)))

update(model, 10000, progress.bar="none")
samp <- coda.samples(model, variable.names=c("br","ar","B","am","bm","pi","mu2"), n.iter=20000, progress.bar="none")
summary(samp)

B = apply(samp[[1]][,1:2],2,median)
plot(X,Y[,1], type="l", ylim=c(0,500))
lines(X,Y[,2], type="l", col="red")
lines(X,max(Y[,1])/(exp(abs(X)/20))+10, type="l", col="blue",lwd=3)

#Bayes Factor
pi = samp[[1]][,grep("pi",colnames(samp[[1]]))]
pi = apply(pi,2,median)
Prior_param = seq(1,2,length=1000)
M1 = function(p,d){dbeta(d,p,p-1)}
M2 = function(p,d){dbeta(d,p-1,p)}
BF = vector()
for(i in 1:length(pi))
{
  BF[i] = 2*log(sum(dbeta(pi[i],Prior_param,Prior_param-1)*dgamma(Prior_param,1,1))/sum(dbeta(pi[i],Prior_param-1,Prior_param)*dgamma(Prior_param,1,1)))
}


points(X[BF>10],Y[BF>10,1], col="red")
points(X[BF>10],Y[BF>10,1], col="red")

# #restrict length
# pi1 = samp1[[1]][,grep("pi",colnames(samp1[[1]]))]
# pi2 = samp[[1]][,grep("pi",colnames(samp[[1]]))]
# #pi1 is larger than pi2, so we need to trim
# pi1_start = ((ncol(pi1)-1001)/2)+1
# pi2_start = ((ncol(pi2)-1001)/2)+1
# pi1 = pi1[,pi1_start:(pi1_start+1000)]
# pi2 = pi2[,pi2_start:(pi2_start+1000)]
# 
# #Now we need to multiply the two probability distributions
# fitbeta <- function(x){
#   m = mean(x)
#   v = var(x)
#   f = v+m^2-m
#   A = -f*m/v
#   B = f*(m-1)/v
#   c(A,B)
# }
# 
# v = vector()
# for(i in 1:ncol(pi1))
# {
#   par1 = fitbeta(pi1[,i])
#   par2 = fitbeta(pi2[,i])
#   par = (par1+par2)/2
#   v[i] = median(rbeta(1000000,par[1],par[2]))
#   print(i)
# }
# 
# 
# pi1_1 = as.vector(apply(pi1,2,median))
# 
# pi2_2 = as.vector(apply(pi2,2,median))
# odds_2 = pi2_2/(1-pi2_2)

