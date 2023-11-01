##check the HiQR package with L1 penalty
rm(list=ls())
library(HiQR)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
set.seed(99)
library('tictoc')
library('glmnet')
library('Matrix')
source('fun.R')

n=200
p=50
## Generate data
X0<-matrix(rnorm(n*(p-1)),nrow=n)
X<-cbind(1,X0)
Omega<-matrix(0,nrow=p,ncol=p)
Omega[1:3,1:3]=rnorm(9)
Omega=(Omega+t(Omega))/2
Y<-diag(X%*%Omega%*%t(X))+rnorm(n)


##all-pairs LASSO
tic('HiQR')
obj0=hiqr(X0,Y)
toc()


##all-pairs LASSO
tic('HiQR')
obj1=hiqr(X0,Y,nlambda1 =25,nlambda2 =25, type =1)
toc()

##all-pairs LASSO
tic('HiQR')
obj2=hiqr(X0,Y,nlambda1 =25,nlambda2 =25,type =2)
toc()

##all-pairs LASSO
tic('HiQR')
obj3=hiqr(X0,Y,nlambda1 =25,nlambda2 =25,type =3)
toc()

obj3$Omega[[420]][1:10,1:10]
Omega[1:10,1:10]
