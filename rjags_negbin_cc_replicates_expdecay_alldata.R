library(rjags)
#let us start a superset from what we have
getMatrix <-function(rep1,rep2,mid)
{
  start = mid - 5E05
  end = mid + 5E05
  chr = names(sort(table(rep1[,1]),decreasing = T)[1])
  rep1 = rep1[rep1[,1]==chr,]
  r1 = rep1[rep1[,4]>start&rep1[,4]<end,]
  rep2 = rep2[rep2[,1]==chr,]
  r2 = rep2[rep2[,4]>start&rep2[,4]<end,]
  r = sort(unique(c(r1[,4],r2[,4])))
  m = matrix(ncol=2,nrow=length(r),1)
  m[match(r1[,4],r),1] = round(r1[,6])
  m[match(r2[,4],r),2] = round(r2[,6])
  row.names(m) <-r
  center = which(m[,1]==max(m[,1]))
  m[(center-500):(center+500),]
}

rep1 = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O/CapC-786O/CapC-786O_REdig_sorted_CC1_chr8:128882978-128883307_chr8+128883062+5.gff")
rep2 = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O/CapC-786O_2/CapC-786O_2_REdig_sorted_CC1_chr8:128882978-128883307_chr8+128883062+5.gff")
m1 = getMatrix(rep1,rep2,128882978)
rep1 = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O/CapC-786O/CapC-786O_REdig_sorted_CC1_chr10:99185458-99186090_PGAM1.gff")
rep2 = read.table("~/Documents/ResearchProjects/CC_Capture_Project/786O/CapC-786O_2/CapC-786O_2_REdig_sorted_CC1_chr10:99185458-99186090_PGAM1.gff")
m2 = getMatrix(rep1,rep2,99185458)
d = array(NA,dim=c(1001,2,2))
d[,,1] = m1
d[,,2] = m2
mid = rep(501,2)
max_point = apply(d[501,,],2,max)

# data1 = round(rep1[,6]/((rep1[,5]-rep1[,4])/1000))+1
# max = which(rep1[,4]<128882978&rep1[,4]>128881000)
# X = (c(1:length(data1))-max)
# window = min(min(max(X[X>0]),max(abs(X[X<0])))-1,500)
# range = c((max-window):(max+window))
# X = X[range]
# data1 = data1[range]
# rep1 = rep1[range,]
# 
# data2 = round(rep2[,6]/((rep2[,5]-rep2[,4])/1000))+1
# max = which(rep2[,4]<128882978&rep2[,4]>128881000)
# X = (c(1:length(data2))-max)
# window = min(min(max(X[X>0]),max(abs(X[X<0])))-1,500)
# range = c((max-window):(max+window))
# X = X[range]
# data2 = data2[range]
# rep2 = rep2[range,]

#We obviously have 

#clean data with less than 20 reads
#plot(X,Y[,1], type="l", ylim=c(0,500))
#lines(X,Y[,2], type="l", col="red")

#Now let us do this with rjags
model_string <- "model{

#Start with the data
for(j in 1:length(Y[1,1,]))
{
  for(i in 1:length(Y[,1,1]))
  {
        #Mean of my data is a poisson GMM, where parameters are distributed according
        #to a gamma prior distribution. Add an extra parameter for the Site.
        mu[j,i] ~ dpois((B[j,1]+(K[j]/exp(sqrt((i-c[j])^2)/B[j,2])))+(odds[j,i]*mu2[j,i]))
        #The mixture coefficient is distributed according to a beta distribution
        odds[j,i] = pi[j,i]/(1-pi[j,i])
        pi[j,i] ~ dbeta(.5,.5)
        p[i,1,j] = r[i,1,j]/(mu[j,i]+r[i,1,j])
        p[i,2,j] = r[i,2,j]/(mu[j,i]+r[i,2,j])
        mu2[j,i] ~ dgamma(am,bm)
        #let us assume that the mean is fixed
        #The scaling parameter is a convolution of two gamma distributions. The site
        #gamma distribution, and the replicates
        r[i,1,j] = site_r[j,i] + rep_r[1]
        r[i,2,j] = site_r[j,i] + rep_r[2]
        #site is sampled from a gamma distribution shrinking within the replicate
        #notice that for every site, we a new set of gamma parameters
        site_r[j,i] ~ dgamma(sar,sbr)
        #Replicate is sampled from prior distribution where we rather shrink 
        #across all sites for the same replicate
        #ar and br itself, should be sampled from a hyper prior 
        #The delaborte distribution, the two replicates are sampled from the same distribution
        Y[i,1,j] ~ dnegbin(p[i,1,j],r[i,1,j])
        Y[i,2,j] ~ dnegbin(p[i,2,j],r[i,2,j])
        #Split spread parameter to 2 convoluted sources, replicates, and sample
    }
    B[j,1] ~ dgamma(B1a,B1b)
    B[j,2] ~ dgamma(B2a,B2b)
}

  #two replicates
  for(x in 1:2)
  {
    rep_r[x] ~ dgamma(ar[x],br[x])
    ar[x] ~ dgamma(1,1)
    br[x] ~ dgamma(1,1)
  }
  
  B1a ~ dgamma(1,1)
  B1b ~ dgamma(1,1)
  B2a ~ dgamma(1,1) 
  B2b ~ dgamma(1,1) 
  sar ~ dgamma(1,1)
  sbr ~ dgamma(1,1)
  am ~ dgamma(1,1)
  bm ~ dgamma(1,1)
}
"

model <- jags.model(textConnection(model_string), data = list(Y=d,c=mid, K=max_point))

update(model, 10000)

samp <- coda.samples(model, variable.names=c("pi","r","rep_r","B"), n.iter=20000)

plotResult <- function(X=c(-500:500),Y,samp,samplenumber, bayesf, chr,coord)
{
  B = samp[[1]][,grep(paste("B\\[",samplenumber,",", sep=""),colnames(samp[[1]]))]
  B = apply(B,2,median)
  plot(X,Y[,1,samplenumber], type="l", ylim=c(0,200))
  lines(X,Y[,2,samplenumber], type="l", col="red")
  lines(X,(B[1]+(max(Y)/exp(sqrt(X^2)/B[2]))), type="l", col="blue",lwd=3)
  
  #Bayes Factor
  pi = samp[[1]][,grep(paste("pi\\[",samplenumber,",", sep=""),colnames(samp[[1]]))]
  pi = apply(pi,2,median)
  Prior_param = seq(1,2,length=1000)
  BF = vector()
  for(i in 1:length(pi))
  {
    BF[i] = (sum(dbeta(pi[i],Prior_param,Prior_param-1)*dgamma(Prior_param,1,1))/sum(dbeta(pi[i],Prior_param-1,Prior_param)*dgamma(Prior_param,1,1)))
  }
  #points(X[BF>=bayesf],apply(Y[BF>=bayesf,,samplenumber],1,max), col="red")
  cc = as.numeric(coord[BF>=bayesf])
  #plot(X,pi, type="l")
  
  points(X[BF>=bayesf],Y[BF>=bayesf,2,samplenumber], col="red")
  cbind(chr,cc,cc+1)
}

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

