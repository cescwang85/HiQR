##simulation for ridge regression with single lambda
rm(list=ls())
library(HiQR)
setwd('~/Dropbox/Share/HQR_ADMM/sim')
set.seed(99)
library('tictoc')
library('glmnet')
library('Matrix')
library('MASS')
ttime<-function(obj)
{t1<-proc.time();
obj
t2<-as.vector(proc.time()-t1);
return(sum(t2[1:2]))}

n=500
p=50
## Parameters
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
## Generate data
X0<-mvrnorm(n,rep(0,p),Sigma)
X<-cbind(1,X0)
Y1<-diag(X%*%Omega1%*%t(X))+rnorm(n);
Y2<-diag(X%*%Omega2%*%t(X))+rnorm(n);
Y3<-diag(X%*%Omega3%*%t(X))+rnorm(n);

##
lambda=1
Y=Y1

##Conducting Ridge regression for QR
tic('naive-inverse')
X2<-apply(X, 1, vec2inter); ##generate p^2*n interaction data matrix
obj<-solve(X2%*%t(X2)/n+lambda*diag(nrow(X2)),X2%*%Y/n)
t1=toc()
A1=matrix(as.vector(obj),ncol=p+1)
A1=(A1+t(A1))/2

tic('Woodbury-trick')
X2<-apply(X, 1, vec2inter); ##generate p^2*n interaction data matrix
b0=X2%*%Y/n;
a2=b0-X2%*%solve(t(X2)%*%X2/n+lambda*diag(ncol(X2)),t(X2)%*%b0/n);
a2=a2/lambda;
t2=toc()
A2=matrix(as.vector(a2),ncol=p+1)
A2=(A2+t(A2))/2

tic('svd trick')
X2<-apply(X, 1, vec2inter); ##generate p^2*n interaction data matrix
obj_svd=svd(X2/sqrt(n))
U<-obj_svd$u;
V<-obj_svd$v;
d<-obj_svd$d;
a3=U%*%((d/(d^2+lambda))*(t(V)%*%Y/sqrt(n)))
t3=toc()
A3=matrix(as.vector(a3),ncol=p+1)
A3=(A3+t(A3))/2
##New Ridge regression for QR
tic('qr')
obj_new=qr2(X,Y,lambda =lambda)
t0=toc()
A0=obj_new$Omega[[1]]

norm(A0-A1,type='F')
max(abs(A0-A1))

norm(A0-A2,type='F')
max(abs(A0-A2))

norm(A0-A3,type='F')
max(abs(A0-A3))

A0[1:5,1:5]



