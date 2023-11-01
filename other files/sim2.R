##time performance of HiQR for LASSO problem
rm(list=ls())
library(HiQR)
setwd('~/Dropbox/Share/HQR_ADMM/sim')
set.seed(99)
library('glmnet')
library('Matrix')
library('MASS')
library("lattice")
library('ncvreg')

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

runtime<-rep(0,3);
## glmnet
if (p<1200){
runtime[1]=ttime({Xfull<-t(apply(X, 1, vec2inter));
                 obj<-glmnet(Xfull[,-1],Y,standardize =FALSE,nlambda=50,lambda.min.ratio =sqrt(log(p)/n));})
lambda=obj$lambda
## ncvreg
runtime[2]=ttime({Xfull<-t(apply(X, 1, vec2inter));
obj1<-ncvreg(Xfull[,-1],Y,penalty ='lasso',lambda=lambda,family='gaussian',returnX =FALSE)})}
## HiQR
runtime[3]=ttime({obj5=hiqr(X0,Y,type=1,nlambda=50,err_abs =min(1e-4,0.1/p),err_rel =min(1e-4,0.1/p),rho=sqrt(p),maxIter =500,lambda.min.ratio =sqrt(log(p)/n))})
print(c(p,runtime,mean(obj5$niter)))
return(runtime)
}
pp<-c(200,400,800,1200,1600,2000,2400)
Re<-sapply(rep(pp,each=10),ff)
Re
saveRDS(Re,file='time-tab2.rds')




