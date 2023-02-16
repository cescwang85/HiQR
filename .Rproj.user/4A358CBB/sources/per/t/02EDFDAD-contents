##check the HiQR package with L1 penalty
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


##Conducting all pairs LASSO
tic('lasso-glmnet')
Xfull<-t(apply(X, 1, ww_new));
obj<-glmnet(Xfull[,-1],Y,standardize =FALSE,nlambda=50,lambda.min.ratio =0.1);
toc()

lambda=obj$lambda

##all-pairs LASSO
tic('qr')
obj_new=qr1(X,Y,lambda =lambda)
toc()

##all-pairs LASSO
tic('HiQR')
obj_f=hiqr(X0,Y)
toc()

t(obj_new$rho) ## check varying rho
t(obj_new$niter) ##check the number of iterations

## check the estimation
NN=length(lambda)
for (k in 1:NN)
{
  A1=matrix(c(obj$a0[k], obj$beta[,k]),nrow=p)
  A1=(A1+t(A1))/2
A2=obj_f$Omega[[k]]
print(c(norm(A1-A2,type='F')/sqrt(p),max(abs(A1-A2))))
}



