#ospq_lasso, the code is based on the official code of lasso Matlab
rm(list=ls())
set.seed(123)
library('osqp')
library('glmnet')
library('tictoc')
##define a new lasso solver based on ospq
ospq_lasso<-function(X,Y,lambda=1)
{
  n=nrow(X);p=ncol(X);
  mX<-colMeans(X);mY=mean(Y);
  X<-scale(X,scale =FALSE);
  Y<-scale(Y,scale=FALSE);
  ##OSQP data
  P=matrix(0,nrow=n+2*p,ncol=n+2*p)
  diag(P)=c(rep(0,p),rep(1/n,n),rep(0,p))
  P=Matrix(P, sparse = TRUE)  
  q=rep(0,n+2*p)
  A=matrix(0,nrow=n+2*p,ncol=n+2*p)
  A[1:n,]=cbind(X,-diag(n),matrix(0,n,p))
  A[(n+1):(n+p),1:p]=diag(p); A[(n+1):(n+p),(n+p+1):(n+2*p)]=-diag(p);
  A[(n+p+1):(n+2*p),1:p]=diag(p); A[(n+p+1):(n+2*p),(n+p+1):(n+2*p)]=diag(p);
  A=Matrix(A, sparse = TRUE)  
  l<-c(Y,rep(-Inf,p),rep(0,p));
  u<-c(Y,rep(0,p),rep(Inf,p));
  
  # Solve
  settings <- osqpSettings(verbose =FALSE,warm_start =TRUE)
  model <- osqp(P, q, A, l, u, settings)
  #res <- model$Solve()
  # Define new vector
  Beta_hat<-matrix(0,nrow=p+1,ncol=length(lambda))
  for (k in 1:length(lambda)){
    q_new <- lambda[k]*c(rep(0,n+p),rep(1,p))
    # Update model and solve again
    model$Update(q = q_new)
    res <- model$Solve()
    Beta_hat[-1,k]=res$x[1:p];
    Beta_hat[1,k]=mY-t(mX)%*%Beta_hat[-1,k];
  }
  return(list(beta=Beta_hat,lambda=lambda))
}

p=500
n=1000
beta_true<-rbinom(p,1,0.05)*rnorm(p)/sqrt(p)
X<-matrix(rnorm(n*p),nrow=n)*matrix(rbinom(n*p,1,0.15),nrow=n)
X<-X*matrix(rbinom(n*p,1,0.25),nrow=n)
Y=as.vector(X%*%beta_true+rnorm(n))
tic()
S=cov(X)
dim(S)
toc()
tic('glmnet')
obj<-glmnet(X,Y,standardize =FALSE,intercept =TRUE,nlambda =50,lambda.min.ratio = 0.01)
toc()
tic('ospq')
obj_ospq<-ospq_lasso(X,Y,lambda=obj$lambda)  
toc()
max(abs(obj$beta-obj_ospq$beta[-1,]))
max(abs(obj$a0-obj_ospq$beta[1,]))
