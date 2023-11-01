##Check the QR-LASSO 
## ADMM不是一个高精度算法，精度不能设置的过小；尝试了变ADMM系数，带来的提升非常有限，简单起见规定rho=5

rm(list=ls())
setwd('~/Dropbox/ChengGroup/HQR_ADMM/qrlasso/')
set.seed(123)
library('tictoc')
library('glmnet')
library('Rcpp')
library('Matrix')
library('hqr')
source('fun.R')
n=200
p=50
## Generate data
X0<-matrix(rnorm(n*(p-1)),nrow=n)
X<-cbind(1,X0)
Omega<-matrix(0,nrow=p,ncol=p)
Omega[1:3,1:3]=rnorm(9)
Y<-diag(X%*%Omega%*%t(X))+rnorm(n)
H0=(X%*%t(X))^2/n

lam=0.2

tic('qr-R')
err_abs=1e-4
err_rel=1e-2
D0=gram(X,Y)/n; 
K=500
ee1=rep(0,K) ##dual residual s^{k+1}=rho (z^{k+1}-z^k)
ee2=rep(0,K)  ## primal residual r^{k+1}=x^{k+1}-z^{k+1}
err_pri=rep(0,K)
err_dual=rep(0,K)
ee3=rep(0,K)
A=matrix(0,p,p) 
B=A;
U=A;
rho=5
H=solve(H0+rho*diag(n))

for (k in 1:K)
{Dk=D0+rho*(B-U)
  #H=HU%*%((1/(Hd+rho))*t(HU))
  ## A-update
  w1=as.vector(H%*%qrow(X,Dk)/n)
  A=(Dk-gram(X,w1))/rho
 # A=(A+t(A))/2
  ##B-update
  oldB=B
  B=soft(A+U,lam/rho)
 # print(B[1,1])
  ##U-update
  U=A-B+U
  ee1[k]=rho*norm(B-oldB,type='F')
  ee2[k]=norm(A-B,type='F')
  err_pri[k]=p*err_abs+err_rel*max(norm(A,type='F'),norm(B,type='F'));
  err_dual[k]=p*err_abs+err_rel*norm(U,type='F')
## Update rho if necessary
#if (ee2[k]>10*ee1[k]) {rho=2*rho;print(c(k,rho))}
#  if (ee1[k]>10*ee2[k]) {rho=rho/2;print(c(k,rho))}
if ((ee1[k]<err_dual[k])&&(ee2[k]<err_pri[k])) break;
  }
toc()
print(k)
plot(log(ee1),col='red')
lines(log(ee2),col='blue')

## check for single lam
tic('qr-Rcpp')
 obj=qr_lasso(X,Y,lambda =lam)
 toc()
 obj$niter
max(abs(obj$Omega[[1]]-B)) 

## Check for a seq lambda
tic('glmnet')
Xfull<-t(apply(X0, 1, ww));
obj1<-glmnet(Xfull,Y,standardize =FALSE);
toc()
tic('qr-Rcpp')
obj2=qr_lasso(X,Y,lambda =obj1$lambda)
toc()

lambda=obj1$lambda
NN=length(lambda)
for (k in 1:NN)
{A1=inww(c(obj1$a0[k], obj1$beta[,k]))
A2=obj2$Omega[[NN+1-k]]
print(norm(A1-A2,type='F'))
}

##check
#B[1:5,1:5]
#B_qr[1:5,1:5]
#B_glmnet[1:5,1:5]

