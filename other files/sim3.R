##time performance of HiQR and FAMILY
rm(list=ls())
library(HiQR)
setwd('~/Dropbox/Share/HQR_ADMM/sim')
set.seed(99)
library('glmnet')
library('Matrix')
library('MASS')
library("lattice")
source('ef_1.R')
source('ADMM_0.1.19.R')

ttime<-function(obj)
{t1<-proc.time();
obj
t2<-as.vector(proc.time()-t1);
return(sum(t2[1:2]))}

n=500
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
  Y=as.vector(diag(X%*%Omega%*%t(X))+rnorm(n))

nalpha=10
nlambda =10

dd<-function(fit,obj){
  re<-matrix(0,nalpha,nlambda);
  re_fit<-matrix(0,nalpha,nlambda);
  re_obj<-matrix(0,nalpha,nlambda);
  for (i in 1:length(alphas)){
    for (j in 1:nlambda){
      #print(obj1$lambda1[i+(j-1)*nalpha]-alphas[i]*lambdas[nlambda+1-j])  
      A1=as.matrix(fit$Estimate[[i]][[nlambda+1-j]]$finB)  ##
      A1=(A1+t(A1))/2
      A2=as.matrix(obj$Omega[[i+(j-1)*nalpha]])
      re[i,j]=norm(A1-A2,type='F');
      re_fit[i,j]=fit$Estimate[[i]][[j]]$iters;
    }}
  return(c(sum(obj$niter),sum(re_fit), summary(as.vector(re))))
}

alphas<- seq(0.1,1,length =nalpha)

runtime<-rep(0,6);
## type=12 (l1+l_inf)
runtime[1]=ttime({obj12=hiqr(X0,Y,type=12,nlambda =nlambda,lambda.min.ratio =0.5,nalpha=nalpha,err_rel = 1e-4)})
## type=13 (l1+l2)
runtime[2]=ttime({obj13=hiqr(X0,Y,type=13,nlambda =nlambda,lambda.min.ratio =0.5,nalpha=nalpha,err_rel = 1e-4)})
## type=14 (l1+l1/l_inf)
runtime[3]=ttime({obj14=hiqr(X0,Y,type=14,nlambda =nlambda,lambda.min.ratio =0.5,nalpha=nalpha,err_rel = 1e-4)})
## type=15 (l1+nuclear norm)
runtime[4]=ttime({obj15=hiqr(X0,Y,type=15,nlambda =nlambda,lambda.min.ratio =0.5,nalpha=nalpha,err_rel = 1e-4)})

## type=12 (l1+l_inf)+FAMILY
lambdas=matrix(obj12$lambda1,ncol=nalpha,byrow=TRUE)[,1]/min(alphas)
lambdas_more=matrix(obj12$lambda2,ncol=nalpha,byrow=TRUE)[,1]/(1-min(alphas))
lambdas=sort(lambdas)
lambdas_more=sort(lambdas_more)
runtime[5]=ttime({fit12<-new_FAMILY(X0, X0, Y, lambdas=lambdas,lambdas_more=lambdas_more,
                                    alphas=alphas, norm ='l_inf',e.abs =1e-4,e.rel =1e-4,verbose =FALSE)})
print(dd(fit12,obj12))
## type=13 (l1+l2)+FAMILY
lambdas=matrix(obj13$lambda1,ncol=nalpha,byrow=TRUE)[,1]/min(alphas)
lambdas_more=matrix(obj13$lambda2,ncol=nalpha,byrow=TRUE)[,1]/(1-min(alphas))
lambdas=sort(lambdas)
lambdas_more=sort(lambdas_more)
runtime[6]=ttime({fit13<-new_FAMILY(X0, X0, Y, lambdas=lambdas,lambdas_more=lambdas_more,
                                    alphas=alphas, norm ='l2',e.abs =1e-4,e.rel =1e-4, verbose =FALSE)})
print(dd(fit13,obj13))
print(c(p,runtime))
return(runtime)
}

pp<-c(50,100,200)
Re<-sapply(rep(pp,each=10),ff)
Re
saveRDS(Re,file='time-tab3.rds')




