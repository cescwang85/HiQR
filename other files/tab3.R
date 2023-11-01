rm(list=ls())
library('stargazer')
setwd('~/Dropbox/Share/HQR_ADMM/sim')
Re<-readRDS('time-tab3.rds')
pp<-c(50,100,200)

dig=3
table<-function(obj){
  m<-format(round(apply(obj,1,mean,na.rm=TRUE),dig),nsmall =dig)
  s<-format(round(apply(obj,1,sd,na.rm=TRUE),dig),nsmall =dig)
  l<-rep('(',length(m));
  r<-rep(')',length(m));
  R1=paste(m,l,s,r,sep='')
  return(R1)
}
AA<-matrix(0,nrow=3,ncol=6)
for ( k in 1:length(pp))
{AA[k,]=table(Re[,(1:10)+(k-1)*10])}

rownames(AA)=c('p=50','p=100','p=200');
stargazer(AA)
BB<-matrix(0,nrow=6,ncol=4)
rownames(BB)<-c(1,1,2,2,3,4)
colnames(BB)<-c('Method','p=50','p=100','p=200')
BB[,1]=c('HiQR','FAMILY','HiQR','FAMILY','HiQR','HiQR')
BB[1,2:4]=AA[,1];
BB[2,2:4]=AA[,5];
BB[3,2:4]=AA[,2];
BB[4,2:4]=AA[,6];
BB[5,2:4]=AA[,3];
BB[6,2:4]=AA[,4];
stargazer(BB)
