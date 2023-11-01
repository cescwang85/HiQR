rm(list=ls())
library('stargazer')
setwd('~/Dropbox/Share/HQR_ADMM/sim')
Re<-readRDS('time-tab1.rds')
pp<-c(100,200,400,800,1200)

dig=3
table<-function(obj){
  m<-format(round(apply(obj,1,mean,na.rm=TRUE),dig),nsmall =dig)
  s<-format(round(apply(obj,1,sd,na.rm=TRUE),dig),nsmall =dig)
  l<-rep('(',length(m));
  r<-rep(')',length(m));
  R1=paste(m,l,s,r,sep='')
  return(R1)
}
AA<-matrix(0,nrow=4,ncol=length(pp))
for ( k in 1:length(pp))
{AA[,k]=table(Re[,(1:10)+(k-1)*10])}
colnames(AA)=c('p=100','p=200','p=400','p=800','p=1200');
rownames(AA)=c('Naive', 'Woodbury','SVD','HiQR')
stargazer(AA)
