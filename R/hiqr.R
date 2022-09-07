#' Main function for the High dimensional Quadratic Regression (HiQR)
#' @param \code{X} covariates data matrix of dimension n*p.
#' @param Y response variable with length n.
#' @param lambda1 user supplied tuning parameter for the \eqn{\ell_1} penalty; Default is NULL and the program compute its own
#' sequence based on \code{nlambda}.
#' @param  nlambda1 the length of the tuning parameter sequence which is available when lambda is NULL. Default is 50.
#' @param  lambda.min.ratio1 smallest value for lambda, as a fraction of lambda.max which is available when lambda is NULL. Default is 0.1.
#' @param err the precision used to stop the convergence. Default is 1e-4.
#' @param maxIter Maximum number of iterations. Default is 1000.
#' @param rho step parameter for the ADMM. Default is 1.
#' @return A sparse asymmetric p*p matrix for the interactions.
#' @examples
#' rm(list = ls())
#'library('PIE')
#'library('glmnet')
#'set.seed(99)
#'p=100
#'n=200;
#'Omega=matrix(0,nrow=p,ncol=p);
#'Omega[6,6]=1
#'Omega[1,6]=2;Omega[6,10]=2;Omega=(Omega+t(Omega))/2;
#'beta<-rep(0,p);beta[c(1,6,10)]=1;
#'X=matrix(rnorm(n*p),nrow =n);
#'Y=as.vector(diag(X%*%Omega%*%t(X))+X%*%beta+rnorm(n));
##Centering
## PIE basd on the response
#'yOme<-PIE(X,Y)
##PIE based on the residual
#'hbeta<-as.vector(coef(cv.glmnet(X,Y,nfolds =5),s="lambda.min"))[-1];  
#'rOme<-PIE(X,Y-X%*%hbeta)
#'yOme[1:10,1:10]
#'rOme[1:10,1:10]
hiqr<-function(X,Y,lambda1=NULL,nlambda1=10,lambda.min.ratio1=0.1,lambda2=NULL,nlambda2=10,lambda.min.ratio2=0.1,type=NULL,err_abs=1e-4,err_rel=1e-3,maxIter=200,rho=5)
{
  X<-cbind(1,as.matrix(X));
  Y<-as.vector(Y);
  n=nrow(X);
  p=ncol(X);
  An<-t(X)%*%diag(Y-mean(Y))%*%X/n;
  if (is.null(lambda1)){
 lambda1<-max(abs(An))*exp(seq(from=log(lambda.min.ratio1),to=0,length.out=nlambda1));
    nlambda1=length(lambda1);}
   if (is.null(type))  {obj<-qr1(X,Y,lambda =lambda1,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho);return(obj)}
  if (is.null(lambda2)){
    if (type==1) {    
      lmax=sqrt(max(apply(An[,-1]*An[,-1],2,sum)))/2;
      lambda2<-lmax*exp(seq(from=log(lambda.min.ratio2),to=0,length.out=nlambda2));
      nlambda2=length(lambda2); }
    if (type==2) {    
      lmax=max(apply(abs(An[,-1]),2,sum))/2;
      lambda2<-lmax*exp(seq(from=log(lambda.min.ratio2),to=0,length.out=nlambda2));
      nlambda2=length(lambda2); }
    if (type==3) {    
      lmax=max(apply(abs(An[2:p,2:p]),1,max)+abs(An[1,2:p]))/2;
      lambda2<-lmax*exp(seq(from=log(lambda.min.ratio2),to=0,length.out=nlambda2));
      nlambda2=length(lambda2); }
  }
  obj<-qr2(X,Y,lambda1=lambda1,lambda2=lambda2,type=type,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho);
  return(obj);
}
