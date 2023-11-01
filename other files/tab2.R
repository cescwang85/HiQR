rm(list=ls())
library('stargazer')
setwd('~/Dropbox/Share/HQR_ADMM/sim')
Re<-readRDS('time-tab2.rds')
pp<-c(200,400,800,1200,1600,2000,2400)

dig=2
table<-function(obj){
  m<-format(round(apply(obj,1,mean,na.rm=TRUE),dig),nsmall =dig)
  s<-format(round(apply(obj,1,sd,na.rm=TRUE),dig),nsmall =dig)
  l<-rep('(',length(m));
  r<-rep(')',length(m));
  R1=paste(m,l,s,r,sep='')
  return(R1)
}
AA<-matrix(0,nrow=3,ncol=length(pp))
for ( k in 1:length(pp))
{AA[,k]=table(Re[,(1:10)+(k-1)*10])}

colnames(AA)=c('p=200','p=400','p=800','p=1200','p=1600','p=2000','p=2400');
rownames(AA)=c('glmnet','ncvreg','HiQR')
stargazer(AA)
