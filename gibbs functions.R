fix.MH=function(lo,hi,old1,new1,jump){
  jold=pnorm(hi,mean=old1,sd=jump)-pnorm(lo,mean=old1,sd=jump)
  jnew=pnorm(hi,mean=new1,sd=jump)-pnorm(lo,mean=new1,sd=jump)
  list(ljold=log(jold),ljnew=log(jnew))
}
#----------------------------------------------------------------------------------------------
tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#-------------------------------
rmvnorm1=function (n, sigma, pre0.9_9994 = FALSE) 
{
#   retval <- chol(sigma, pivot = TRUE)
#   o <- order(attr(retval, "pivot"))
#   retval <- retval[, o]
  s. <- svd(sigma)
  if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
    warning("sigma is numerically not positive definite")
  }
  R = t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
  retval
}
#----------------------------
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<10000
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------
sample.sig2=function(param,not.nas,a.sig2,b.sig2,ntime,nsites,y){
  a1=(not.nas+2*a.sig2)/2
  media=matrix(param$betas,ntime,nsites,byrow=T)+param$etas%*%t(param$lambda)
  err3=apply((y-media)^2,2,sum,na.rm=T)

  b1=b.sig2+(err3/2)
  1/rgamma(nsites,a1,b1)
}
#----------------------------
sample.lambda=function(param,nsites,ntime,y){
  H=param$etas
  nfactors=length(param$tau)
  lambda=matrix(NA,nsites,nfactors)
  media=matrix(param$betas,ntime,nsites,byrow=T)
  err=y-media
  sig2=param$sig2
  
  for (i in 1:nsites){
    cond1=!is.na(y[,i])
    H1=H[cond1,]
    HtH=t(H1)%*%H1

    tmp=(1/param$phi[i,])*(1/param$tau)
    invD=diag(1/tmp)
    
    #to avoid numerical issues
    cond=invD>1e+13
    invD[cond]=1e+13
    
    invVar=invD+(1/sig2[i])*HtH
    Var=solve(invVar)
    # cond=rcond(invVar)<1e-13 #assess reciprocal of the condition number of a matrix
    # if (!cond) Var=solve(invVar)
    # if ( cond){ #matrix that  is really close to being singular
    #   #see matrix inverse document
    #   A=invD; invA=diag(tmp)
    #   U=(1/param$sig2[i])*t(H)
    #   C=diag(rep(1,ncol(U)))
    #   V=H
    #   Var=invA-invA%*%U%*%solve(C+V%*%invA%*%U)%*%V%*%invA
    #   #hist(Var-Varz)
    # }
    
    pmedia=(1/sig2[i])*t(H1)%*%err[cond1,i]
    lambda[i,]=rmvnorm1(1,sigma=Var)+t(Var%*%pmedia)
  }
  lambda
}
#----------------------------
sample.etas=function(param,ntime,nsites,y){
  nfactors=length(param$tau)    
  etas=matrix(NA,ntime,nfactors)
  medias=matrix(param$betas,ntime,nsites,byrow=T)
  err=y-medias
  sig2=param$sig2

  for (i in 1:ntime){
    cond=!is.na(y[i,])
    Lambda=param$lambda[cond,]
    if (sum(cond)==1) {
      dim(Lambda)=c(1,nfactors)
      invSigma=1/sig2[cond]
    }
    if (sum(cond)>1) invSigma=diag(1/sig2[cond])
    invVar=t(Lambda)%*%invSigma%*%Lambda+diag(rep(1,nfactors))
    Var=solve(invVar)
    pmedia=t(Lambda)%*%invSigma
    media=Var%*%pmedia%*%err[i,cond]
    etas[i,]=rmvnorm1(1,Var)+t(media)
  }
  etas
}
#----------------------------
sample.phi=function(param,nsites,v){
  a1=(v+1)/2
  nfactors=length(param$tau)
  tau=matrix(param$tau,nsites,nfactors,byrow=T)
  b1=(tau*(param$lambda^2)+v)/2
  tmp=rgamma(nsites*nfactors,a1,b1)
  matrix(tmp,nsites,nfactors)
}
#----------------------------
sample.delta1=function(param,nsites){
  nfactors=length(param$delta)
  a1=param$a1+((nsites*nfactors)/2)
  aux=param$phi*(param$lambda^2)
  soma=colSums(aux)
  
  tmp=rep(NA,nfactors)
  tmp[1]=1
  for (i in 2:nfactors){
    tmp[i]=prod(param$delta[2:i])
  }
  b1=1+(sum(tmp*soma)/2)
  rgamma(1,a1,b1)
}
#----------------------------
sample.deltah=function(param,nsites){
  nfactors=length(param$delta)
  aux=param$phi*(param$lambda^2)
  soma=colSums(aux)
  delta=param$delta
  
  for (i in 2:nfactors){
    a1=param$a2+(nsites*(i+1)/2)


    delta[i]=1 #this has the effect of eliminating the hth factor    
    p1=rep(NA,nfactors)
    for (j in 1:nfactors) p1[j]=prod(delta[1:j])
    
    soma1=sum(p1[i:nfactors]*soma[i:nfactors])
    b1=1+(soma1/2)
    delta[i]=rgamma(1,a1,b1)
  }
  delta[-1]
}
#----------------------------
sample.a1=function(param,jump){
  a1.old=param$a1
  a1.new=tnorm(1,lo=0,hi=Inf,mu=a1.old,sig=jump)
  pold=(a1.old-1)*log(param$delta[1])-lgamma(a1.old)+log(a1.old)-a1.old
  pnew=(a1.new-1)*log(param$delta[1])-lgamma(a1.new)+log(a1.new)-a1.new
  
  tmp=fix.MH(lo=0,hi=Inf,a1.old,a1.new,jump)
  k=acceptMH(pold+tmp$ljnew,pnew+tmp$ljold,a1.old,a1.new,F)
  list(a1=k$x,accept=k$x!=a1.old)
}
#----------------------------
sample.a2=function(param,jump){
  a2.old=param$a2
  a2.new=tnorm(1,lo=1,hi=Inf,mu=a2.old,sig=jump)
  nfactor=length(param$delta)
  pold=(a2.old-1)*sum(log(param$delta[-1]))-(nfactor-1)*lgamma(a2.old)+log(a2.old)-a2.old
  pnew=(a2.new-1)*sum(log(param$delta[-1]))-(nfactor-1)*lgamma(a2.new)+log(a2.new)-a2.new
  tmp=fix.MH(lo=1,hi=Inf,a2.old,a2.new,jump)
  
  k=acceptMH(pold+tmp$ljnew,pnew+tmp$ljold,a2.old,a2.new,F)
  list(a2=k$x,accept=k$x!=a2.old)
}
#------------------------------
calc.tau=function(delta){
  nfactors=length(delta)
  tau=rep(NA,nfactors)
  tau[1]=delta[1]
  for (i in 2:nfactors) tau[i]=prod(delta[1:i])
  tau
}
#-------------------------------
#remove redundant factors
sample.factor=function(i,param,neighb,nsites,ntime,v){
  # thresh=ifelse(1,i<500,0) 
  thresh=exp(-1-5*(10^(-4))*i)
  u=runif(1)
  nfactors=length(param$delta)
  if (u<thresh & nfactors!=2){ #I don't want to end up with a single factor
    tmp=colSums(abs(param$lambda)<neighb)
    cond=tmp==nrow(param$lambda)
    if (sum(cond)>0) {
      param$lambda=param$lambda[,!cond]
      param$etas=param$etas[,!cond]
      param$delta=param$delta[!cond]
      param$tau=calc.tau(param$delta)
      param$phi=param$phi[,!cond]
    }
    
    if (sum(cond)==0 & nfactors<(nsites-1)){ #just add another factor if we have not reached the maximum
      param$etas=cbind(param$etas,rnorm(ntime))
      param$delta=c(param$delta,rgamma(1,param$a2,1))
      param$tau=calc.tau(param$delta)
      param$phi=cbind(param$phi,rgamma(nsites,v/2,v/2))
      nfactors=ncol(param$phi)
      var1=(1/param$phi[,nfactors])*(1/param$tau[nfactors])
      param$lambda=cbind(param$lambda,rnorm(nsites,mean=0,sd=sqrt(var1)))
    }
  }
  list(delta=param$delta,tau=param$tau,phi=param$phi,lambda=param$lambda,etas=param$etas)
}
#-----------------------------------
#sample location specific seasonal parameters
sample.betas=function(param,nsites,y){
  invT=1
  betas=rep(NA,nsites)
  sig2=param$sig2
  
  for (i in 1:nsites){
    cond=!is.na(y[,i])
    n=sum(cond)
    prec=(n/sig2[i])+invT
    var1=1/prec
    pmedia=(1/sig2[i])*sum(y[cond,i]-param$etas[cond,]%*%t(t(param$lambda[i,])))
    betas[i]=rnorm(1,var1*pmedia,sqrt(var1))
  }
  betas
}
#-----------------------------------
make.predict=function(param,ntime,nsites,y){
  Lambda=param$lambda
  sig2=param$sig2
  
  Sigma=diag(sig2)
  Var1=Lambda%*%t(Lambda)+Sigma
  media=matrix(param$betas,ntime,nsites,byrow=T)
  res=y1=y
  for (i in 1:ntime){
    cond=is.na(y1[i,])
    if (!(sum(cond)%in%c(0,nsites))){ #either all data are present or all are missing
      mis=which(cond)
      obs=which(!cond)
      tmp=crmvnorm(media[i,],Var1,obs,mis,y1[i,])      
      res[i,cond]=tmp
    }
  }
  res
}
#-----------------------------------
#conditional multivariate normal
crmvnorm=function(media,Var1,obs,mis,y){
  Sigma.obs=Var1[obs,obs]
  if (length(obs)==1) Sigma.obs=matrix(Sigma.obs,1,1)
  invSigma.obs=solve(Sigma.obs)
  
  Sigma.mis=Var1[mis,mis]
  if (length(mis)==1) Sigma.mis=matrix(Sigma.mis,1,1)

  Sigma.mis.obs=Var1[mis,obs]
  if (length(mis)==1 | length(obs)==1) Sigma.mis.obs=matrix(Sigma.mis.obs,length(mis),length(obs))
  
  media1=media[mis]+Sigma.mis.obs%*%invSigma.obs%*%(y[obs]-media[obs])
  var1=Sigma.mis-Sigma.mis.obs%*%invSigma.obs%*%t(Sigma.mis.obs)
  t(rmvnorm1(1,var1))+media1
}