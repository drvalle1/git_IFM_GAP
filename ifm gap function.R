ifm.gap=function(y,ngibbs,ind.impact.sites,nadapt){
  not.nas=apply(!is.na(y),2,sum)
  nsites=ncol(y)
  ntime=nrow(y)
  nfactors=floor(5*log(nsites))
  
  etas=matrix(rnorm(ntime*nfactors),ntime,nfactors)
  phi=matrix(rgamma(nsites*nfactors,1,1),nsites,nfactors)
  delta=rep(1.1,nfactors)
  a1=10
  a2=5
  v=3 #changed this
  
  tau=rep(NA,nfactors)
  for (i in 1:nfactors) tau[i]=prod(delta[1:i])
  plot(tau,type='h')
  
  lambda=matrix(rnorm(nsites*nfactors),nsites,nfactors) #how each site is influenced by each factor
  betas=rep(0,nsites)
  #boxplot(lambda)
  
  a.sig2=b.sig2=0.1
  sig2=rep(0.1,nsites)
  
  #MCMC stuff
  vec.as=matrix(NA,ngibbs,2)
  vec.sig2=matrix(NA,ngibbs,nsites)
  vec.pred=matrix(NA,ngibbs,ntime*length(impacted.sites))
  vec.factors=matrix(NA,ngibbs,1)
  vec.prec=matrix(NA,ngibbs,nsites*nsites)
  vec.betas=matrix(NA,ngibbs,nsites)

  param=list(sig2=sig2,lambda=lambda,etas=etas,phi=phi,delta=delta,
             a1=a1,a2=a2,tau=tau,betas=betas)
  jump1=list(a1=10,a2=1)
  accept1=list(a1=0,a2=0)
  accept.output=50
  
  neighb=10^(-3) 
  options(warn=2) #this should be 2
  for (i in 1:ngibbs){
    tmp=length(param$delta)
    print(c(i,tmp))
    vec.factors[i]=tmp
    
    param$sig2=sample.sig2(param=param,not.nas=not.nas,a.sig2=a.sig2,
                           b.sig2=b.sig2,ntime=ntime,nsites=nsites,y=y)
    param$etas=sample.etas(param=param,ntime=ntime,nsites=nsites,y=y)
    param$lambda=sample.lambda(param=param,nsites=nsites,ntime=ntime,y=y)
    
    param$phi=sample.phi(param=param,nsites=nsites,v=v)
    #param$phi[]=1 #to avoid numerical issues in sample.lambda when param$phi[i,]*param$tau
    
    param$delta[1]=sample.delta1(param=param,nsites=nsites)
    param$delta[-1]=sample.deltah(param=param,nsites=nsites)
    param$tau=calc.tau(delta=param$delta)
    
    tmp=sample.a1(param=param,jump=jump1$a1)
    param$a1=tmp$a1
    accept1$a1=accept1$a1+tmp$accept
    #   param$a1=1
    
    tmp=sample.a2(param=param,jump=jump1$a2)
    param$a2=tmp$a2
    accept1$a2=accept1$a2+tmp$accept
    #   param$a2=2
    
    param$betas=sample.betas(param=param,nsites=nsites,y=y)
    # param$betas[]=0
    
    tmp=sample.factor(i=i,param=param,neighb=neighb,nsites=nsites,ntime=ntime,v=v)
    param$etas=tmp$etas
    param$delta=tmp$delta
    param$tau=tmp$tau
    param$phi=tmp$phi
    param$lambda=tmp$lambda
    
    if (i%%accept.output==0 & i<nadapt){
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
    
    tmp=make.predict(param=param,ntime=ntime,nsites=nsites,y=y)
    vec.pred[i,]=tmp[,ind.impact.sites]
    vec.as[i,]=c(param$a1,param$a2)
    vec.sig2[i,]=param$sig2
    vec.betas[i,]=param$betas
    vec.prec[i,]=solve(param$lambda%*%t(param$lambda)+diag(param$sig2))
  }
  list(pred=vec.pred,as=vec.as,sig2=vec.sig2,betas=vec.betas,
       prec=vec.prec,nfactors=vec.factors,
       lambda=param$lambda,delta=param$delta,tau=param$tau) 
  #notice that we just have the last iteration for
  #delta, tau, lambda
  #because size keeps on changing
}