rm(list=ls(all=TRUE))
set.seed(123120)

setwd('U:\\independent studies\\kaplan amazon hydrology\\git_IFM_GAP')
source('gibbs functions.R')
source('ifm gap function.R')
tmp=read.csv('fake data.csv',as.is=T)

impacted.sites=c('w2360','w2370','w2371')
ind.impact.sites=which(colnames(tmp)%in%impacted.sites)
ngibbs=1000

ind=which(colnames(tmp)%in%c('ano','mes'))
y=data.matrix(tmp[,-ind])

res=ifm.gap(y=y,ind.impact.sites=ind.impact.sites,ngibbs=ngibbs,nadapt=1000)

