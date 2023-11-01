##simulation for ridge regression with single lambda
rm(list=ls())
library(HiQR)
setwd('~/Dropbox/Share/HQR_ADMM/sim')
set.seed(99)
library('glmnet')
library('Matrix')
library('MASS')
ttime<-function(obj)
{t1<-proc.time();
obj
t2<-as.vector(proc.time()-t1);
return(sum(t2[1:2]))}

n=500
lambda=10
## Generate data
ff<-function(p){
## Parameters
Sigma=toeplitz(0.5^(0:(p-1)))
Omega<-matrix(0,nrow=p+1,ncol=p+1)
Omega[2,6]=3;Omega[6,6]=-2.5;Omega[6,11]=4
Omega[1,2]=1;Omega[2,1]=1
Omega[1,11]=1;Omega[11,1]=1

X0<-mvrnorm(n,rep(0,p),Sigma)
X<-cbind(1,X0)
Y<-diag(X%*%Omega%*%t(X))+rnorm(n);

##
runtime<-rep(0,4);
##naive-inverse
if (p<=100){
runtime[1]=ttime({
X2<-apply(X, 1, vec2inter); ##generate p^2*n interaction data matrix
obj<-solve(X2%*%t(X2)/n+lambda*diag(nrow(X2)),X2%*%Y/n)
})}

## Woodbury-trick
if (p<=1000){
runtime[2]=ttime({
X2<-apply(X, 1, vec2inter); ##generate p^2*n interaction data matrix
b0=X2%*%Y/n;
a2=b0-X2%*%solve(t(X2)%*%X2/n+lambda*diag(ncol(X2)),t(X2)%*%b0/n);
a2=a2/lambda;
})}


## svd trick
if (p<=1000){
runtime[3]=ttime({
X2<-apply(X, 1, vec2inter); ##generate p^2*n interaction data matrix
obj_svd=svd(X2/sqrt(n))
U<-obj_svd$u;
V<-obj_svd$v;
d<-obj_svd$d;
a3=U%*%((d/(d^2+lambda))*(t(V)%*%Y/sqrt(n)))})}
##New Ridge regression for QR
runtime[4]=ttime({obj_new=qr2(X,Y,lambda =lambda)})
print(c(p,runtime))
return(runtime)
}

pp<-c(100,200,400,800,1200)
Re<-sapply(rep(pp,each=10),ff)
saveRDS(Re,file='time-tab1.rds')

