rm(list=ls(all=TRUE))
set.seed(1)
library('mvtnorm')

setwd('U:/independent studies/kaplan amazon hydrology/hydrology data')
dat=read.csv('final edited.csv',as.is=T)
nsites=ncol(dat)-2
ntime=nrow(dat)

#create variance-covariance matrix
sig2.true=sig2=runif(nsites)
Sigma=diag(sig2)

p=5 #number of factors
Lambda=matrix(0,nsites,p)
delta=c(1.5,2,2,2,2)
for (i in 1:p){
  prec1=prod(delta[1:i])
  print(c(i,prec1))
  Lambda[,i]=rnorm(nsites,mean=0,sd=sqrt(1/prec1)) #this ignores phi's
}
boxplot(Lambda)
LtL=Lambda%*%t(Lambda)
Var1.true=Var1=LtL+Sigma

#create regression parameters (one for each location)
beta.true=beta=rnorm(nsites)

#create response variable
y=y.true=rmvnorm(ntime,mean=rep(0,nrow(Var1)),sigma=Var1)+
         matrix(beta,ntime,nsites,byrow=T)

#create some holes randomly
ind=sample(1:(nsites*ntime),size=(nsites*ntime)*0.1)
y[ind]=NA

#create data
y1=cbind(y,dat$ano,dat$mes)
colnames(y1)=colnames(dat)

#create NAs for impacted sites
tempo.pred=150:ntime
impacted.sites=c('w2360','w2370','w2371')
y1[tempo.pred,impacted.sites]=NA
image(is.na(y1))

setwd('U:/independent studies/kaplan amazon hydrology/git_IFM_GAP')
write.csv(y1,'fake data.csv',row.names=F)