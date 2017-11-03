rm(list=ls(all=TRUE))
set.seed(123120)

setwd('U:\\independent studies\\kaplan amazon hydrology\\git_IFM_GAP')
source('gibbs functions.R')
source('ifm gap function.R')

setwd('U:/independent studies/kaplan amazon hydrology/hydrology data')
tmp=read.csv('final edited.csv',as.is=T)

# estreito: em construcao
# lajeado (construcao 1998-2002), operacao 2001
# peixe angical: operacao 2006
# sao salvador: operacao 2009

#put NA's after 1997 for sites impacted by Lajeado upstream (sites 2360,2370,2371)
#lajeado started operating in 2001-2002 (http://dams-info.org/en/dams/view/lajeado/; http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1679-62252007000200005)
#construction period 1998-2002 (https://pt.wikipedia.org/wiki/Usina_Hidrel%C3%A9trica_Luiz_Eduardo_Magalh%C3%A3es)
nomes=c(paste('w',c(2360,2370,2371),sep=''),
        paste('w',c(2360,2370,2371),'_lag',sep=''))
cond=tmp$ano>1997
tmp[cond,nomes]=NA

#put NA's for sites likely to have been impacted by Tucurui
#not impacted are: 2885, 2370, 2371, 2360
#impacted: 2920, 2970, 2975, (perhaps 2910, 2905)
impacted=c(2920, 2970, 2975)

#remove lags of impacted sites
ind=which(colnames(tmp)%in%paste('w',impacted,'_lag',sep=''))
tmp=tmp[,-ind]

#plug in NA's for impacted period of impacted sites
impacted.sites=which(colnames(tmp)%in%paste('w',impacted,sep=''))
cond=tmp$ano >= 1979
tmp[cond,impacted.sites]=NA
image(data.matrix(tmp))

ind=which(colnames(tmp)%in%c('ano','mes'))
y=data.matrix(tmp[,-ind])

ngibbs=50000
res=ifm.gap(y=y,ind.impact.sites=impacted.sites,ngibbs=ngibbs,nadapt=1000)

setwd('U:/independent studies/kaplan amazon hydrology/gibbs/results')
write.csv(res$as,'vec as.csv',row.names=F)
write.csv(res$sig2,'vec sig2.csv',row.names=F)
write.csv(res$factors,'vec factors.csv',row.names=F)
write.csv(res$betas,'vec betas.csv',row.names=F)

ind=seq(from=ngibbs*(9/10),to=ngibbs,length.out=1000)
write.csv(res$prec[ind,],'vec precision.csv',row.names=F)
write.csv(res$pred[ind,],'vec pred.csv',row.names=F)