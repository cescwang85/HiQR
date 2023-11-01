##time performance of HiQR 
rm(list=ls())
library(HiQR)
setwd('~/Dropbox/Share/HQR_ADMM/sim')
set.seed(123)
library('glmnet')
library('Matrix')
library('MASS')
library("lattice")
library('ncvreg')
library('tictoc')

ttime<-function(obj)
{t1<-proc.time();
obj
t2<-as.vector(proc.time()-t1);
return(sum(t2[1:2]))}

n=500
## Generate data
p=1000
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

    runtime[1]=ttime({Xfull<-t(apply(X, 1, vec2inter));
    obj<-cv.glmnet(Xfull[,-1],Y,standardize =FALSE,nlambda=50,lambda.min.ratio=log(p)/n)})
  plot(obj)
    lambda=obj$lambda
    ## ncvreg
    runtime[2]=ttime({Xfull<-t(apply(X, 1, vec2inter));
    obj1<-ncvreg(Xfull[,-1],Y,penalty ='lasso',lambda=lambda,family='gaussian',returnX =FALSE)})
  ## HiQR
  runtime[3]=ttime({obj2=hiqr(X0,Y,type=1,nlambda=50,err_abs =1e-05,err_rel = 1e-4,lambda.min.ratio =sqrt(log(p)/n),rho=2)})
  print(c(p,runtime))
sum(obj2$niter)
# check the estimation
NN=length(lambda)
for (k in 1:NN)
{
  A1=matrix(c(obj$a0[k], obj$beta[,k]),nrow=p+1)
  A1=(A1+t(A1))/2
  A2=matrix(obj1$beta[,k],nrow =p+1)
  A2=(A2+t(A2))/2
A3=as.matrix(obj2$Omega[[k]])
print(c(max(abs(A1-A2)),max(abs(A1-A3))))
}

alphas<-0.99999
#alphas=0.9999
lambdas<- lambda
nlambdas<-length(lambdas)
lambda1=as.vector(lambdas[nlambdas:1]%*%t(alphas))
lambda2=as.vector(lambdas[nlambdas:1]%*%t(1-alphas))*sqrt(p)
tic('l_inf+FAMILY')
fit1<- FAMILY(X0, X0, Y, lambdas , alphas, norm ='l_inf', e.abs = 1e-5, e.rel = 1e-5,verbose =FALSE)
toc()
tic('l_inf+HiQR')
obj1<-hiqr(X0,Y,type=14,lambda =lambda1,lambda2 =lambda2,err_abs = 1e-5,err_rel = 1e-5)
toc()

tic('l2+FAMILY')
fit2<- FAMILY(X0, X0, Y, lambdas , alphas, norm ='l2',e.abs = 1e-5, e.rel = 1e-5,verbose =FALSE)
toc()
tic('l2+HiQR')
obj2<-hiqr(X0,Y,type=13,lambda =lambda1,lambda2 =lambda2,err_abs = 1e-5,err_rel = 1e-5)
toc()
i=1
  for (j in 1:length(lambdas)){
    A0=matrix(c(obj$a0[j], obj$beta[,j]),nrow=p+1)
    A0=(A0+t(A0))/2
    
    #print(lambda1[(nlambdas+1-j)+(i-1)*nlambdas]-alphas[i]*lambdas[j])  
    A1=as.matrix(fit2$Estimate[[i]][[j]]$finB)  ##
    print(fit2$Estimate[[i]][[j]]$iters)
    A1=(A1+t(A1))/2
    A2=as.matrix(obj2$Omega[[(nlambdas+1-j)+(i-1)*nlambdas]])
    A1[1,1]=A2[1,1]
    aa=norm(A1-A2,type='F')
    bb=norm(A0-A2,type='F')
   print(round(c(aa,bb,csi(Omega,A1),csi(Omega,A2)),3))
  }



