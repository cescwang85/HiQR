##check the HiQR package with L1 penalty
#library("devtools")
#devtools::install_github("cescwang85/PIE")


rm(list=ls())
library(HiQR)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
set.seed(99)
library('tictoc')
library('glmnet')
library('Matrix')
library('FAMILY')
library('PIE')

p=100
n=200;
Omega=matrix(0,nrow=p,ncol=p);
Omega[6,6]=1
Omega[1,6]=2;Omega[6,10]=2;Omega=(Omega+t(Omega))/2;
beta<-rep(0,p);beta[c(1,6,10)]=1;
X=matrix(rnorm(n*p),nrow =n);
Y=as.vector(diag(X%*%Omega%*%t(X))+X%*%beta+rnorm(n));
##Centering
## PIE basd on the response
yOme<-PIE(X,Y)
##PIE based on the residual
hbeta<-as.vector(coef(cv.glmnet(X,Y,nfolds =5),s="lambda.min"))[-1];  
rOme<-PIE(X,Y-X%*%hbeta)
##Show the estimation
yOme[1:10,1:10]
rOme[1:10,1:10]