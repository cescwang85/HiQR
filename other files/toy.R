##check the glmnet for very large p
rm(list=ls())
library(HiQR)
set.seed(99)
library('tictoc')
library('glmnet')
library('Matrix')


n=500
p=500
## Generate data
X0<-matrix(rnorm(n*(p-1)),nrow=n)
X<-cbind(1,X0)
Omega<-matrix(0,nrow=p,ncol=p)
Omega[1:3,1:3]=rnorm(9)
Omega=(Omega+t(Omega))/2
Y<-diag(X%*%Omega%*%t(X))+rnorm(n)
##Conducting all pairs LASSO fail (n=p=1000)
# tic('lasso-glmnet')
# Xfull<-t(apply(X, 1, vec2inter));
# obj<-glmnet(Xfull[,-1],Y,standardize =FALSE);
# toc()

# ## Naive Inverse 
lambda=10
# 
 X2<-apply(X, 1, vec2inter); ##generate p^2*n interaction data matrix
##naive-inverse（fail n=500, p=200） 
#obj<-solve(X2%*%t(X2)/n+lambda*diag(nrow(X2)),X2%*%Y/n)

## Woodbury-trick
    b0=X2%*%Y/n;
    a2=b0-X2%*%solve(t(X2)%*%X2/n+lambda*diag(ncol(X2)))%*%(t(X2)%*%b0/n);
    a2=a2/lambda;


## svd trick
    obj_svd=svd(X2/sqrt(n))
    U<-obj_svd$u;
    V<-obj_svd$v;
    d<-obj_svd$d;
    a3=U%*%((d/(d^2+lambda))*(t(V)%*%Y/sqrt(n)))
##New Ridge regression for QR
obj_new=qr2(X,Y,lambda =lambda)

AA=as.matrix(obj_new$Omega[[1]])

A2=matrix(as.vector(a2),nrow =p)
max(abs(AA-A2))


A3=matrix(as.vector(a3),nrow =p)
max(abs(AA-A2))
max(abs(AA-A3))

    