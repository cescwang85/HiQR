##探索变系数rho对结果的影响
rm(list=ls())
set.seed(99)
setwd('~/Dropbox/Share/HQR_ADMM/sim')
library(HiQR)
library('tictoc')
library('glmnet')
library('MASS')
library('Matrix')
library('FAMILY')
source('hiqr_R_version.R')
p=200
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


##Step 1 检查确认R版本和Rcpp版本结果一致
# obj=hiqr(X0,Y,type=1,nlambda=50,err_abs =1e-05,err_rel = 1e-5,lambda.min.ratio =sqrt(log(p)/n),rho=5,maxIter =1000)
# lambda=obj$lambda
# for (k in 1:length(lambda)){print(nnzero(obj$Omega[[k]]))}
#obj_R=qr1_R(X,Y,lambda=lambda,err_abs =1e-05,err_rel = 1e-5,rho=5)

# max(abs(obj$lambda-obj_R$lambda))
# max(abs(obj$niter-obj_R$niter))
# for (k in 1:length(lambda)){
# print(max(abs(obj$Omega[[k]]-obj_R$Omega[[k]])))
#}

## Step 2 对于LASSO检查变系数带来的影响
lam=2


loss<-function(X,Y,B,lambda1,lambda2=0)
{
  re=mean((Y-diag(X%*%B%*%t(X)))^2)/2+lambda1*(sum(abs(B))-abs(B[1,1]))#+lambda2*max(svd(B)$d);
  return(re)
}

maxIter =500
aa=qr1_single(X,Y,lambda =lam,err_abs =0,err_rel = 0,maxIter =maxIter,rho=sqrt(p),rho_vary =0)
f1<-rep(0,maxIter)
for (k in 1:maxIter)
{f1[k]=loss(X,Y,aa$Omega[[k]],lambda1=lam,lambda2=0)
}
plot(1:maxIter,log(f1-min(f1)),type='l')


plot(1:maxIter,log(aa$abs_err),type='l')
lines(1:maxIter,log(aa$pri_err),type='l',col='red')
lines(1:maxIter,log(aa$dual_err),type='l',col='blue')


## Step 3 对于Sparse Low Rank 检查变系数带来的影响
lambda2=3
bb=qr3_rank_single(X,Y,lambda1 =lam,lambda2=lambda2,err_abs =0,err_rel = 0,maxIter =maxIter,rho=1,rho_vary =0)
f2<-rep(0,maxIter)
for (k in 1:maxIter)
{f2[k]=loss(X,Y,bb$Omega[[k]],lambda1=lam,lambda2=lambda2)
}
min(f2)
plot(1:maxIter,log(f2-min(f2)),type='l')
bb$rho




