rm(list=ls())
set.seed(99)
setwd('~/Dropbox/Share/HQR_ADMM/sim/')
source('ef_1.R')
source('ADMM_0.1.19.R')
library(HiQR)
library('tictoc')
library('glmnet')
library('MASS')
library('Matrix')



p=50
n=500;
Sigma=toeplitz(0.5^(0:(p-1)))
Omega3<-matrix(0,nrow=p+1,ncol=p+1)
Omega3[2,6]=3;Omega3[6,6]=-2.5;Omega3[6,11]=4
Omega3=(Omega3+t(Omega3))/2;
Omega2=Omega3;
Omega2[1,6]=-1;Omega2[6,1]=-1;
Omega1=Omega2;
Omega1[1,2]=1;Omega1[2,1]=1
Omega1[1,11]=1;Omega1[11,1]=1
Omega3[c(1,2,6,11),c(1,2,6,11)] ##check matrix
X0<-mvrnorm(n,rep(0,p),Sigma)
X<-cbind(1,X0)
Y<-diag(X%*%Omega1%*%t(X))+rnorm(n);
obj2<-hiqr(X0,Y,type=1)

##update-D
p=50
A<-matrix(0,p,p)
g.list=c(list(A),list(A),list(A))
B<-matrix(rnorm(p^2),nrow=p)
lambda<-rep(0.1,3)
aaD<-update_D(y=Y,W=X,B=B,g.list =g.list,rho=1,lambda =lambda, norm = "l_inf")
aaE<-update_E(y=Y,W=X,B=B,g.list =g.list,rho=1,lambda =lambda, norm = "l_inf")
aaF<-update_F(y=Y,W=X,B=B,g.list =g.list,rho=1,lambda =lambda)
bbD<-proc(t(B),lambda =0.1,type=2);bbD=t(bbD)
bbE<-proc(B,lambda =0.1,type=2)
bbF<-proc(B,lambda =0.1,type=1)

max(abs(aaD-bbD))
max(abs(aaE-bbE))
max(abs(aaF-bbF))