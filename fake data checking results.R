compare=function(true,estim){
  rango=range(c(true,estim))
  plot(true,estim,ylim=rango,xlim=rango)
  lines(rango,rango)
}
ind=max(which(!is.na(res$betas[,1])))
compare(beta.true,res$betas[ind,])
compare(sig2.true,res$sig2[ind,])
#-----------------------------------------
tmp=res$as
plot(tmp[,1],type='l')
plot(tmp[,2],type='l')

plot(res$nfactors,type='l')
boxplot(res$lambda)

plot(log(res$tau),type='h')
plot(res$delta[-1],type='h')
#-----------------------------------------
#Look at predicted results for impacted sites after impact
pred=matrix(res$pred[ind,],ntime,length(impacted.sites))
compare(y.true[tempo.pred,ind.impact.sites],pred[tempo.pred,])
#-------------------------------
#Look at the same results separately for each impacted site
nimp=length(impacted.sites)
for (i in 1:nimp) {
  compare(y.true[tempo.pred,ind.impact.sites[i]],pred[tempo.pred,ind.impact.sites[i]])
}

#-------------------------------
#Look at induced covariance
var1=solve(matrix(res$prec[ind,],nsites,nsites))
compare(Var1.true,var1)

#-------------------------------
#Look at induced correlation
var1=solve(matrix(res$prec[ind,],nsites,nsites))
mat1=matrix(sqrt(diag(var1)),nsites,nsites)
mat2=matrix(sqrt(diag(var1)),nsites,nsites,byrow=T)
cor1=var1/(mat1*mat2)

mat1=matrix(sqrt(diag(Var1.true)),nsites,nsites)
mat2=matrix(sqrt(diag(Var1.true)),nsites,nsites,byrow=T)
cor1.true=Var1.true/(mat1*mat2)
compare(cor1.true,cor1)