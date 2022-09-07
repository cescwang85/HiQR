##check the hqr package with single lambda rm(list=ls())

setwd('~/Dropbox/ChengGroup/HQR_ADMM/qrlasso/')
set.seed(99)
library('tictoc')
library('hqr')
library('glmnet')
library('FAMILY')

x<-rnorm(5)
x2<-sqrt(sum(x^2))
lambda=1
y1<-l2_vec(x,lambda=lambda)
y2<-(x2>lambda)*(1-lambda/x2)*x
max(abs(y1-y2))

l1inf_vec_R<-function(x,lambda=0){
x1<-sort(abs(x[-1]),decreasing =TRUE);
if (lambda>=abs(x[1])+x1[1]) return(rep(0,length(x)));
lambda1=lambda-max((cumsum(x1)+lambda-abs(x[1]))/(2:length(x)))
lambda1=lambda1*(lambda1>0)
lambda1=lambda+(lambda1-lambda)*(lambda1<lambda)
y<-abs(x)
y<-(y-(lambda-lambda1))*(y<(lambda-lambda1))+(lambda-lambda1)
y[1]=(abs(x[1])-lambda1)*(abs(x[1])<lambda1)+lambda1
y=y*sign(x)
return(list(result=x-y,lambda1=lambda1))}



x<-rnorm(10)
x[1]=10
lambda =0.5
aa=l1inf_vec_R(x,lambda=lambda)$result
bb=as.vector(l1inf_vec(x,lambda=lambda))
cc<-l_inf_prox(x,lambda)
max(abs(aa-bb))
max(abs(bb-cc))

