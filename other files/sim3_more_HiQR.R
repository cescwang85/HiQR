##toy example to compare HiQR and FAMILY with HiQR's lambdas
rm(list=ls())
library(HiQR)
setwd('~/Dropbox/Share/HQR_ADMM/sim')
set.seed(99)
library('glmnet')
library('Matrix')
library('MASS')
library('tictoc')
library("lattice")
source('ef_1.R')
source('ADMM_0.1.19.R')

p=50
n=100;
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


nalpha<-9;
nlambda<-50;
alphas<- seq(0.1,1,length =nalpha)



tic('l_inf+HiQR')
obj1<-hiqr(X0,Y,type=13,nlambda =nlambda,nalpha=nalpha)
toc()


lambdas=matrix(obj1$lambda1,ncol=nalpha,byrow=TRUE)[,1]/min(alphas)
lambdas_more=matrix(obj1$lambda2,ncol=nalpha,byrow=TRUE)[,1]/(1-min(alphas))
lambdas=sort(lambdas)
lambdas_more=sort(lambdas_more)
tic('l_inf+FAMILY')
fit2<-new_FAMILY(X0, X0, Y, lambdas=lambdas,lambdas_more=lambdas_more,alphas=alphas, 
                 norm ='l2', verbose =FALSE)
toc()


dd<-function(fit,obj){
  re<-matrix(0,nalpha,nlambda)
for (i in 1:length(alphas)){
  for (j in 1:nlambda){
   #print(obj1$lambda1[i+(j-1)*nalpha]-alphas[i]*lambdas[nlambda+1-j])  
    A1=as.matrix(fit$Estimate[[i]][[nlambda+1-j]]$finB)  ##
    A1=(A1+t(A1))/2
    A2=as.matrix(obj$Omega[[i+(j-1)*nalpha]])
    re[i,j]=norm(A1-A2,type='F')
  }}
  return(as.vector(re))
  }
summary(dd(fit2,obj1))

