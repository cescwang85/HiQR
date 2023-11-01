##toy example to compare HiQR and FAMILY with FAMILY's lambdas
##FAMILY问题1. F[-1,-1]作者本意是去除第一个元素，实际效果去除了第一行第一列 
rm(list=ls())
library(HiQR)
set.seed(99)
library('tictoc')
library('glmnet')
library('MASS')
library('Matrix')
library('FAMILY')

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



alphas<- seq(0.1,0.9,length =9)
lambdas<- seq(0.5,1,length =25)
nlambdas<-length(lambdas)
lambda1=as.vector(lambdas[nlambdas:1]%*%t(alphas))
lambda2=as.vector(lambdas[nlambdas:1]%*%t(1-alphas))*sqrt(p)
tic('l_inf+FAMILY')
fit1<- FAMILY(X0, X0, Y, lambdas , alphas, norm ='l_inf',  e.abs = 1e-4, e.rel = 1e-4, verbose =FALSE)
toc()
tic('l_inf+HiQR')
obj1<-hiqr(X0,Y,type=12,lambda =lambda1,lambda2 =lambda2,err_abs = 1e-4,err_rel = 1e-4)
toc()



fit<-fit1
obj<-obj1
for (i in 1:length(alphas)){
  for (j in 1:length(lambdas)){
   #print(lambda1[(nlambdas+1-j)+(i-1)*nlambdas]-alphas[i]*lambdas[j])  
    A1=as.matrix(fit$Estimate[[i]][[j]]$finB)  ##
    A1=(A1+t(A1))/2
    A2=as.matrix(obj$Omega[[(nlambdas+1-j)+(i-1)*nlambdas]])
    aa=norm(A1-A2,type='F')
    print(round(c(fit$Estimate[[i]][[j]]$iters,aa,csi(Omega1,A1),csi(Omega1,A2)),3))
  }}

tic('l2+FAMILY')
fit2<- FAMILY(X0, X0, Y, lambdas , alphas, norm ='l2',e.abs = 1e-4, e.rel = 1e-3,verbose =FALSE)
toc()
fit2$time
tic('l2+HiQR')
obj2<-hiqr(X0,Y,type=13,lambda =lambda1,lambda2 =lambda2,err_abs = 1e-4,err_rel = 1e-3)
toc()

fit<-fit2
obj<-obj2
for (i in 1:length(alphas)){
  for (j in 1:length(lambdas)){
    #print(lambda1[(nlambdas+1-j)+(i-1)*nlambdas]-alphas[i]*lambdas[j])  
    A1=as.matrix(fit$Estimate[[i]][[j]]$finB)  ##
    A1=(A1+t(A1))/2
    A2=as.matrix(obj$Omega[[(nlambdas+1-j)+(i-1)*nlambdas]])
    aa=norm(A1-A2,type='F')
    print(round(c(fit$Estimate[[i]][[j]]$iters,aa,csi(Omega1,A1),csi(Omega1,A2)),3))
  }}

