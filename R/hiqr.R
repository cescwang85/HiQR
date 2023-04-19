#' Main function for the High dimensional Quadratic Regression (HiQR)
#' @param \code{X}, covariates data matrix of dimension \eqn{n \times p}.
#' @param \code{Y}, response variable with length \eqn{n}.
#' @param \code{type}, The penalty to use: 
#' 1(Default)--\eqn{\ell_1} penalty(LASSO);
#' 2--squared \eqn{\ell_2} penalty(Ridge Regression);
#' 5--nuclear norm penalty(Reduced Rank Regression);
#' 12--\eqn{\ell_1+\ell_\infty} penalties;
#' 13-\eqn{\ell_1+\ell_2} penalties;
#' 14--\eqn{\ell_1+\ell_1/\ell_\infty} penalties;
#' 15--\eqn{\ell_1+}nuclear norm penalties.
#' @param \code{lambda1}, user supplied tuning parameter for the first penalty; Default is NULL and the program compute its own
#' sequence based on \code{nlambda1}.
#' @param  \code{nlambda1}, the length of the tuning parameter sequence which is available when \code{lambda1} is NULL. Default is 50.
#' @param  \code{lambda.min.ratio1}, smallest value for lambda, as a fraction of lambda.max which is available when \code{lambda1} is NULL. Default is 0.1.
#' @param \code{lambda2}, user supplied tuning parameter for the second penalty(if exists); Default is NULL and the program compute its own
#' sequence based on \code{nlambda1}.
#' @param  \code{nlambda2}, the length of the tuning parameter sequence which is available when \code{lambda2} is NULL. Default is 50.
#' @param  \code{lambda.min.ratio2}, smallest value for lambda, as a fraction of lambda.max which is available when \code{lambda2} is NULL. Default is 0.1.
#' @param \code{err_abs},\code{err_rel}, the absolute and relative precision used to stop the convergence. Default are 1e-4 and 1e-3.
#' @param \code{maxIter}, Maximum number of iterations. Default is 200.
#' @param \code{rho}, step parameter for ADMM. Default is 5.
#' @return A list with components
#' \item{Omega}{a list with length \code{nlambda1}(for single penalty) or \code{nlambda1}*\code{nlambda2}(for two penalties) of sparse \eqn{(p+1) \times (p+1)} matrices. 
#'For two penalties, the \eqn{(i-1)}*\code{nlambda2}\eqn{+j} element of the list is for (\code{lambda1}[i],\code{lambda2}[j]}).
#' \item{lambda1}{the used lambda1 for the solution path.}
#' \item{lambda2}{the used lambda2 for the solution path.}
#' \item{niter}{the number of iterations for each element which is a \code{nlambda1} vector(for single penalty) or a \code{nlambda1}*\code{nlambda2} matrix (for two penalties).}
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
hiqr<-function(X,Y,type=NULL,lambda1=NULL,nlambda1=50,lambda.min.ratio1=0.1,lambda2=NULL,nlambda2=50,lambda.min.ratio2=0.1,err_abs=1e-4,err_rel=1e-3,maxIter=200,rho=5)
{
  X<-cbind(1,as.matrix(X));
  Y<-as.vector(Y);
  n=nrow(X);
  p=ncol(X);
  An<-t(X)%*%diag(Y-mean(Y))%*%X/n;
  if (type==1){
    if (is.null(lambda1)){lambda1<-max(abs(An))*exp(seq(from=log(lambda.min.ratio1),to=0,length.out=nlambda1));}  
    obj<-qr1(X,Y,lambda=lambda1,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho);  
  }
  if (type==2){
    if (is.null(lambda1)){lambda1<-5*max(abs(An))*exp(seq(from=log(lambda.min.ratio1),to=0,length.out=nlambda1));}  
    obj<-qr2(X,Y,lambda=lambda1);  
  }
  if (type==5){
    if (is.null(lambda1)){lambda1<-max(eigen(An)$values)*exp(seq(from=log(lambda.min.ratio1),to=0,length.out=nlambda1));}  
    obj<-qr1_rank(X,Y,lambda=lambda1,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho); 
  }
  if (type==12){
    if (is.null(lambda1)){lambda1<-max(abs(An))*exp(seq(from=log(lambda.min.ratio1),to=0,length.out=nlambda1));}  
    if (is.null(lambda2)){lambda2<-sqrt(max(apply(An[,-1]*An[,-1],2,sum)))/2*exp(seq(from=log(lambda.min.ratio2),to=0,length.out=nlambda2));}  
    obj<-qr3(X,Y,lambda1=lambda1,lambda2=lambda2,type=2,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho); 
  }
  if (type==13){
    if (is.null(lambda1)){lambda1<-max(abs(An))*exp(seq(from=log(lambda.min.ratio1),to=0,length.out=nlambda1));}  
    if (is.null(lambda2)){lambda2<-max(apply(abs(An[,-1]),2,sum))/2*exp(seq(from=log(lambda.min.ratio2),to=0,length.out=nlambda2));}  
    obj<-qr3(X,Y,lambda1=lambda1,lambda2=lambda2,type=3,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho); 
  }
  if (type==14){
    if (is.null(lambda1)){lambda1<-max(abs(An))*exp(seq(from=log(lambda.min.ratio1),to=0,length.out=nlambda1));}  
    if (is.null(lambda2)){lambda2<-max(apply(abs(An[2:p,2:p]),1,max)+abs(An[1,2:p]))/2*exp(seq(from=log(lambda.min.ratio2),to=0,length.out=nlambda2));}  
    obj<-qr3(X,Y,lambda1=lambda1,lambda2=lambda2,type=4,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho); 
  }
  if (type==15){
    if (is.null(lambda1)){lambda1<-max(abs(An))*exp(seq(from=log(lambda.min.ratio1),to=0,length.out=nlambda1));}  
    if (is.null(lambda2)){lambda2<-max(eigen(An)$values)*exp(seq(from=log(lambda.min.ratio2),to=0,length.out=nlambda2));}  
    obj<-qr3_rank(X,Y,lambda1=lambda1,lambda2=lambda2,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho); 
  }
  return(obj);
}
