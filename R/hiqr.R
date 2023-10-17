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
#' @param \code{lambda}, user supplied tuning parameter for the first penalty; Default is NULL and the program compute its own
#' sequence based on \code{nlambda} and \code{lambda.min.ratio}.
#' @param \code{lambda2}, user supplied tuning parameter for the second penalty(if exists) which has the same length as \code{lambda};
#' Default is NULL and the program compute its own sequence based on \code{nlambda}, \code{lambda.min.ratio}, \code{nalpha} and \code{alpha.min.ratio}.
#' @param  \code{nlambda}, the length of the tuning parameter sequence which is available when \code{lambda} is NULL. Default is 50.
#' @param  \code{lambda.min.ratio}, smallest value for lambda, as a fraction of lambda.max which is available when \code{lambda} is NULL. Default is 0.25.
#' @param  \code{nalpha}, the length of the second parameter to control the magnitude of penalties on individual covariates versus interaction terms.
#' which is available when \code{lambda} is NULL. Default is 5.
#' @param  \code{alpha.min.ratio}, smallest value for \code{alpha}. Default is 0.1.
#' @param \code{err_abs},\code{err_rel}, the absolute and relative precision used to stop the convergence. Default are 1e-4 and 1e-3.
#' @param \code{maxIter}, Maximum number of iterations. Default is 200.
#' @param \code{rho}, step parameter for ADMM. Default is 5.
#' @param \code{rho_vary}, whether varying penalty parameter \eqn{\rho} for ADMM. Default is 0 (no varying).
#' @return A list with components
#' \item{Omega}{a list with length \code{nlambda}(for single penalty) or \code{nlambda}*\code{nalpha}(for two penalties) of sparse \eqn{(p+1) \times (p+1)} matrices. 
#'For two penalties, the element of the list is corresponding to (\code{lambda1}[k],\code{lambda2}[k]).}
#' \item{\code{lambda1}}{the used \code{lambda1} for the solution path.}
#' \item{\code{lambda2}}{the used \code{lambda2} for the solution path.}
#' \item{niter}{the number of iterations for each element which is a \code{nlambda1} vector(for single penalty) or a \code{nlambda1}*\code{nlambda2} vector (for two penalties).}
#' @examples
#' rm(list = ls())
#'library('HiQR')
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
#' obj1<-hiqr(X,Y,type=1) ##LASSO
#' obj2<-hiqr(X,Y,type=2) ## Ridge Regression 
#'obj5<-hiqr(X,Y,type=5)## Reduced Rank Regression
#'obj12<-hiqr(X,Y,type=12) ## $\ell_1+\ell_2$ norm
#'obj13<-hiqr(X,Y,type=13) ## $\ell_1+\ell_\infty$ norm
#'obj14<-hiqr(X,Y,type=14) ## $\ell_1+\ell_1/\ell_\infty$ norm
#'obj15<-hiqr(X,Y,type=15) ## $\ell_1+\ell_*$ norm
hiqr<-function(X,Y,type=1,lambda=NULL,lambda2=NULL,nlambda=50,lambda.min.ratio=0.25,nalpha=5,alpha.min.ratio=0.1,err_abs=1e-4,err_rel=1e-3,maxIter=200,rho=5,rho_vary=0)
{
  X<-cbind(1,as.matrix(X));
  Y<-as.vector(Y);
  n=nrow(X);
  p=ncol(X);
  An<-t(X)%*%diag(Y-mean(Y))%*%X/n;
  if (type==1){
    if (is.null(lambda)){lambda<-max(abs(An))*exp(seq(from=0,to=log(lambda.min.ratio),length.out=nlambda));}  
    obj<-qr1(X,Y,lambda=lambda,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho,rho_vary=rho_vary);  
  }
  if (type==2){
    if (is.null(lambda)){lambda<-5*max(abs(An))*exp(seq(from=0,to=log(lambda.min.ratio),length.out=nlambda));}  
    obj<-qr2(X,Y,lambda=lambda);  
  }
  if (type==5){
    if (is.null(lambda)){lambda<-max(svd(An)$d)*exp(seq(from=0,to=log(lambda.min.ratio),length.out=nlambda));}  
    obj<-qr1_rank(X,Y,lambda=lambda,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho,rho_vary=rho_vary); 
  }
  if (type==12){lambda1=lambda;
    if (is.null(lambda)){
    lambda<-exp(seq(from=0,to=log(lambda.min.ratio),length.out=nlambda));
    alpha<-seq(from=alpha.min.ratio,to=1,length.out=nalpha);
    lambdafull<-rep(lambda,each=length(alpha));
    alphafull<-rep(alpha,times=length(lambda));
    lambda1<-max(abs(An))*(lambdafull*alphafull);
    lambda2<-sqrt(max(apply(An[,-1]*An[,-1],2,sum)))/2*(lambdafull*(1-alphafull));}  
    obj<-qr3(X,Y,lambda1=lambda1,lambda2=lambda2,type=2,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho,rho_vary=rho_vary); 
  }
  if (type==13){lambda1=lambda;
    if (is.null(lambda)){
      lambda<-exp(seq(from=0,to=log(lambda.min.ratio),length.out=nlambda));
      alpha<-seq(from=alpha.min.ratio,to=1,length.out=nalpha);
      lambdafull<-rep(lambda,each=length(alpha));
      alphafull<-rep(alpha,times=length(lambda));
       lambda1<-max(abs(An))*(lambdafull*alphafull);
    lambda2<-mean(apply(abs(An[,-1]),2,sum))/2*(lambdafull*(1-alphafull));} ##using the mean of all rows
    obj<-qr3(X,Y,lambda1=lambda1,lambda2=lambda2,type=3,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho,rho_vary=rho_vary); 
  }
  if (type==14){lambda1=lambda;
    if (is.null(lambda)){
      lambda<-exp(seq(from=0,to=log(lambda.min.ratio),length.out=nlambda));
      alpha<-seq(from=alpha.min.ratio,to=1,length.out=nalpha);
      lambdafull<-rep(lambda,each=length(alpha));
      alphafull<-rep(alpha,times=length(lambda));
      lambda1<-max(abs(An))*(lambdafull*alphafull);
      lambda2<-mean(apply(abs(An[2:p,2:p]),1,max)+abs(An[1,2:p]))/2*(lambdafull*(1-alphafull));}   ##using the mean of all rows
    obj<-qr3(X,Y,lambda1=lambda1,lambda2=lambda2,type=4,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho,rho_vary=rho_vary); 
  }
  if (type==15){lambda1=lambda;
    if (is.null(lambda)){
      lambda<-exp(seq(from=0,to=log(lambda.min.ratio),length.out=nlambda));
      alpha<-seq(from=alpha.min.ratio,to=1,length.out=nalpha);
      lambdafull<-rep(lambda,each=length(alpha));
      alphafull<-rep(alpha,times=length(lambda));
      lambda1<-max(abs(An))*(lambdafull*alphafull);
      lambda2<-max(svd(An)$d)*(lambdafull*(1-alphafull));}  
    obj<-qr3_rank(X,Y,lambda1=lambda1,lambda2=lambda2,err_abs=err_abs,err_rel=err_rel,maxIter=maxIter,rho=rho,rho_vary=rho_vary); 
  }
  return(obj);
}
