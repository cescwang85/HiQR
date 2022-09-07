##check the hqr package with single lambda 
rm(list=ls())
library(HiQR)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
set.seed(99)
library('tictoc')
library('glmnet')
library('Matrix')
source('fun.R')

n=500
p=100
## Generate data
X0<-matrix(rnorm(n*(p-1)),nrow=n)
X<-cbind(1,X0)
Omega<-matrix(0,nrow=p,ncol=p)
Omega[1:3,1:3]=rnorm(9)
Omega=(Omega+t(Omega))/2
Y<-diag(X%*%Omega%*%t(X))+rnorm(n)


##Conducting Ridge regression for QR
lambda=10
tic('lasso-glmnet')
Xf1<-t(apply(X, 1, ww_new));
obj<-solve(t(Xf1)%*%Xf1/n+lambda*diag(ncol(Xf1)),t(Xf1)%*%Y/n)
toc()
A1=matrix(as.vector(obj),ncol=p)
A1=(A1+t(A1))/2
##New Ridge regression for QR
tic('qr')
obj_new=ridge_qr(X,Y,lambda =lambda)
toc()
A2=obj_new$Omega[[1]]

norm(A1-A2,type='F')/sqrt(p)
max(abs(A1-A2))





