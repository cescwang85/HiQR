#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Element-wise soft-thresholding function for a vector
//' @description Soft-thresholding function for a vector which is the solution of
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_1.}
//' @param \code{x}, a vector.
//' @param \code{lambda}, a scalar.
//' @return a vector after threholding.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::vec soft_vec(arma::vec x,double lambda){
  lambda=std::abs(lambda);
  arma::vec y=(x>=lambda)%(x-lambda)+(x<=(-lambda))%(x+lambda); 
  return y;}

//' @title Element-wise soft-thresholding function for a matrix
//' @description Soft-thresholding function for a matrix which is the solution of
//' \deqn{\argmin_Y ||Y-X||_2^2/2+\lambda ||Y||_1.}
//' @param \code{X}, a matrix.
//' @param \code{lambda}, scalar.
//' @param  \code{k}, Not penalized the first element. Default is \code{1}(TRUE).
//' 
//' @return a matrix after threholding.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat soft(arma::mat X,double lambda,int k=1){
  lambda=std::abs(lambda);
  arma::mat Y=(X>=lambda)%(X-lambda)+(X<=(-lambda))%(X+lambda); 
  if (k==1){Y(0,0)=X(0,0);} 
  return Y;}



//' @title Proximal projection of a \eqn{\ell_\infty} penalty. 
//' @description Proximal projection solution of a vector which is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_\infty.}
//' @param \code{x}, a vector.
//' @param \code{lambda}, a scalar.
//' @return a vector after proximal projection.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::vec inf_vec(arma::vec x,double lambda=0){
  lambda=std::abs(lambda);
  arma::vec x1=arma::abs(x);
  arma::uvec x_order=arma::sort_index(x1); 
  x1=arma::sort(x1);
  int p=x.n_elem;
  arma::vec y1=x1;
arma::vec ak=arma::sum(x1)-arma::cumsum(x1)-arma::linspace(p-1,0,p)%x1;
    if (lambda>=ak[0]) {y1=arma::zeros(p,1); return y1;}  
      else {int k=arma::max(arma::linspace(1,p,p)%(ak>lambda));
        for(arma::uword i=k; i<x.n_elem; ++i) {y1(i) = x1[k-1]+(ak[k-1]-lambda)/(p-k);}
      y1(x_order)=y1;
      y1=arma::sign(x)%y1;
      return(y1);}}

//' @title Proximal projection of a \eqn{\ell_\infty} penalty for a matrix 
//' @description Proximal projection solution of a Matrix where each column is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_\infty.}
//' @param \code{X}, a Matrix.
//' @param \code{lambda}, a scalar.
//' @param \code{k}, Not penalized the first \code{k} columns. Default is 1.
//' @return a matrix.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat inf_mat(arma::mat X,double lambda=0,int k=1){
  arma::mat Y=X;
  int p=X.n_cols;
  for(arma::uword i=k; i<p; ++i) {Y.col(i)=inf_vec(X.col(i),lambda);}
  return(Y);}


//' @title Proximal projection of a \eqn{\ell_2} penalty (Group LASSO). 
//' @description Proximal projection solution of a vector which is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_2.}
//' @param \code{x}, a vector.
//' @param \code{lambda}, a scalar.
//' @return a vector.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::vec l2_vec(arma::vec x,double lambda=0){
  lambda=std::abs(lambda);
  arma::vec y=x;
  double x2=arma::norm(x, 2);
if (x2>0) y=(x2>lambda)*(1-lambda/x2)*x;
    return(y);
  }

//' @title Proximal projection of a \eqn{\ell_2} penalty for a matrix 
//' @description Proximal projection solution of a Matrix where each column is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_2.}
//' @param \code{X}, a Matrix.
//' @param \code{lambda}, a scalar.
//' @param \code{k}, Not penalized the first \code{k} columns. Default is 1.
//' @return A Matrix.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat l2_mat(arma::mat X,double lambda=0,int k=1){
  arma::mat Y=X;
  int p=X.n_cols;
  for(arma::uword i=k; i<p; ++i) {Y.col(i)=l2_vec(X.col(i),lambda);}
  return(Y);}




//' @title Proximal projection of a hybrid \eqn{\ell_1/\ell_\infty} penalty for a vector 
//' @description Proximal projection solution of a vector which is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda  \max\{|y_1|, \sum_{i=2}^p |y_i|\}.}
//' @param \code{x}, a vector.
//' @param \code{lambda}, a scalar.
//' @return a vector.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::vec l1inf_vec(arma::vec x,double lambda=0){
  lambda=std::abs(lambda);
  int p=x.n_elem;
  arma::vec x1=x.subvec(1,p-1);
  x1=arma::sort(arma::abs(x1),"descend");
  if (lambda>=std::abs(x(0))+x1(0)) return arma::zeros(p);
  arma::vec y=arma::cumsum(x1)+(lambda-std::abs(x(0)))*arma::ones(p-1);
  double lambda1=lambda-max(y/arma::linspace(2,p,p-1));
  lambda1=lambda1*(lambda1>0);
  lambda1=lambda+(lambda1-lambda)*(lambda1<lambda);
  y=soft_vec(x,lambda-lambda1);
  y(0)=(x(0)>lambda1)*(x(0)-lambda1)+(x(0)<(-lambda1))*(x(0)+lambda1); 
  return(y);}


//' @title Proximal projection of a hybrid \eqn{\ell_1/\ell_\infty} penalty for a matrix 
//' @description Proximal projection solution of a Matrix where each column is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda  \max\{|y_1|, \sum_{i=2}^p |y_i|\}.}
//' @param \code{X}, a matrix.
//' @param \code{lambda}, a scalar.
//' @param \code{k}, Not penalized the first \code{k} columns. Default is 1.
//' @return a matrix.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat l1inf_mat(arma::mat X,double lambda=0,int k=1){
  arma::mat Y=X;
  int p=X.n_cols;
  for(arma::uword i=k; i<p; ++i) {Y.col(i)=l1inf_vec(X.col(i),lambda);}
  return(Y);}

//' @title Proximal projection of a nuclear norm penalty.
//' @description Proximal projection solution of a matrix with a nuclear norm penalty: 
//' \deqn{\argmin_Y ||Y-X||_2^2/2+\lambda \|Y\|_*.}
//' @param \code{X}, a matrix.
//' @param \code{lambda}, a scalar.
//' @return a matrix.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat nuclear_mat(arma::mat X,double lambda){
  lambda=std::abs(lambda);
  arma::vec d;
  arma::mat u;
  arma::mat v;
  svd(u,d,v,X);
  arma::mat D=X.zeros();
  D.diag()=soft_vec(d,lambda);
  arma::mat Y=u*D*v.t(); 
  return Y;}

//' @title Proximal projection of a Matrix with penalty. 
//' @description Proximal projection solution of a Matrix where each column is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda \cdot p(y).}
//' or proximal projection of nuclear norm.
//' @param \code{X}, a matrix.
//' @param \code{lambda}, a scalar.
//' @param \code{type}, The penalty to use. \eqn{1} (Default) is the \eqn{\ell_1} penalty; 
//' \eqn{2} is the \eqn{\ell_\infty} penalty; \eqn{3} is the \eqn{\ell_2} penalty;
//'  \eqn{4} is the hybrid \eqn{\ell_1/\ell_\infty} penalty.
//' @param \code{k}, Not penalized the first element (for \eqn{\ell_1} penalty) or the first \code{k} columns (for other penalties). Default is 1.
//' Not work for the nuclear norm. 
//' @return a matrix.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat proc(arma::mat X,double lambda=0,int type=1,int k=1){
  arma::mat Y=X;
  if (type==1) Y=soft(X,lambda,k);
  if (type==2) Y=inf_mat(X,lambda,k);
  if (type==3) Y=l2_mat(X,lambda,k);
  if (type==4) Y=l1inf_mat(X,lambda,k);
  return(Y);}


//' @title Weighted Gram matrix
//' @description Calculate the weighted Gram matrix:\eqn{X^\top diag(w)X}.
//' @param \code{X},  a \eqn{n \times p} matrix.
//' @param \code{w}, a \eqn{n} dimensional weight vector.
//' @return a \eqn{p \times p} matrix.
//' 
//' @export
// [[Rcpp::export]]
arma::mat gram(arma::mat X,arma::vec w){
  arma::mat Y=X.each_col()%w;
  Y=X.t()*Y;
  return Y;}

//' @title Quadratic form for each sample 
//' @description Calculate the quadratic form for each sample with a Weight matrix \code{W}: 
//' \deqn{ diag(X W X^\top)=(X_1^\top W X_1,\cdots,X_n^\top W X_n)^\top.}
//' @param \code{X},  a \eqn{n \times p} matrix.
//' @param \code{W}, a \eqn{p \times p} weight matrix.
//' @return a \eqn{n} dimensional vector.
//' 
//' @export
// [[Rcpp::export]]
arma::vec qrow(arma::mat X,arma::mat W){
  arma::vec B=arma::sum(((X*W)%X),1);
  return B;}

//' @title Quadratic regression with squared \eqn{\ell_2} penalty (Ridge regression)  
//' @description Algorithm for High dimensional Quadratic Regression(HiQR) with squared \eqn{\ell_2} penalty:
//' \deqn{\argmin_Y \frac{1}{2n}\sum_{i=1}^n (Y_i-X_i^\top \Omega X_i)^2+\lambda \|\Omega\|_2^2.}
//' @param \code{X}, a \eqn{n \times p} data matrix.
//' @param \code{Y}, a \eqn{n} dimensional response vector.  
//' @param \code{lambda}, user supplied tuning parameter; 
//' @return A list with components
//' \item{Omega}{a list of sparse \eqn{p \times p} matrices corresponding to lambda.}
//' \item{lambda}{the used lambda.}
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List qr2(arma::mat X,arma::vec Y,arma::vec lambda){
  int n=X.n_rows;
  lambda=arma::sort(lambda,"descend");
  int nlambda=lambda.size();
  /*Centering int*/
  arma::mat D=gram(X,Y/n);
  arma::mat H0=X*X.t();
  arma::mat H2=H0%H0/n;
  arma::mat eigvec;
  arma::vec eigval;
  eig_sym(eigval,eigvec,H2); 
  arma::vec w0=eigvec.t()*qrow(X,D/n);
  arma::vec w=w0;
 Rcpp::List Omega_all(nlambda);
  double lam;
for (int k=0;k<nlambda;++k) {
    lam=lambda(k);
    w=eigvec*((1/(eigval+lam))%w0);
    Omega_all(k)=arma::mat(D/lam-gram(X,w/lam));
  }
  return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                            Rcpp::Named("lambda") =lambda); }


//' @title Quadratic regression with \eqn{\ell_1} penalty 
//' @description ADMM algorithm for high dimensional Quadratic regression with a \eqn{\ell_1} norm penalty:
//' \deqn{\argmin_Y \frac{1}{2n}\sum_{i=1}^n (Y_i-X_i^\top \Omega X_i)^2+\lambda \|\Omega\|_1.} 
//' @param \code{X}, a \eqn{n*p} input data matrix.
//' @param \code{Y}, a \eqn{n} response vector.  
//' @param \code{lambda}, user supplied tuning parameter.
//' @param \code{err_abs}, \code{err_rel}   the precision used to stop the convergence of ADMM. 
//' @param \code{maxIter}, Maximum number of iterations. Default is 1000.
//' @param \code{rho}, initial step parameter for ADMM.
//' @return A list with components
//' \item{Omega}{a list of sparse p*p matrices corresponding to lambda.}
//' \item{lambda}{the used lambda for the solution path.}
//' \item{niter}{the number of iterations for each element of lambda.}
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List qr1(arma::mat X,arma::vec Y,arma::vec lambda,double err_abs=10^(-4),double err_rel=10^(-3),int maxIter=200,double rho=1){
  int n=X.n_rows;
  int p=X.n_cols;
  lambda=arma::sort(lambda,"descend");
  int nlambda=lambda.size();
  /*Centering int*/
  arma::mat D0=gram(X,Y/n);
  arma::mat Dk=D0;
  arma::mat H0=X*X.t();
  arma::mat H2=H0%H0/n;
  arma::mat H=inv_sympd(H2+rho*arma::eye(n,n));
  Rcpp::List Omega_all(nlambda);
  arma::vec niter=lambda;
  arma::vec rholist=lambda;
  /*Initialization*/
  arma::mat A=arma::zeros(p,p);
  arma::mat U=arma::zeros(p,p);
  arma::mat B=A;
  arma::mat old_B=B;
  arma::vec w=Y;
  double lam;
  double ee_pri=1;
  double ee_dual=1;
for (int k=0;k<nlambda;++k) {
   lam=lambda(k);
    int i=0;
    while ((i<maxIter)||(i==0))
    {
      Dk=D0+rho*(B-U);
      /* A-update*/
      w=H*qrow(X,Dk/n);
      A=(Dk-gram(X,w))/rho;
      /*B-update*/
      old_B=B;
      B=soft(A+U,lam/rho);
      /*U-update*/
        U=A-B+U;
i=i+1;
/*Stop Rule*/
ee_dual=rho*norm(B-old_B,"fro"); /*dual residual*/
ee_pri=norm(A-B,"fro"); /*primal residual*/
double err_pri_new=p*err_abs+err_rel*std::max(norm(A,"fro"),norm(B,"fro"));
double err_dual_new=p*err_abs+err_rel*norm(U,"fro");
if ((ee_dual<err_dual_new)&&(ee_pri<err_pri_new)) break;
/*Varying Penalty Parameter*/
/*if (ee_pri>10*ee_dual) {rho=2*rho;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
/*if (ee_dual>10*ee_pri) {rho=rho/2;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
}
Omega_all(k)=arma::sp_mat(B);
niter(k)=i;
rholist(k)=rho;
}
  return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                            Rcpp::Named("lambda") =lambda,
                            Rcpp::Named("rho") =rholist,
                            Rcpp::Named("niter") =niter); }

//' @title Quadratic regression with \eqn{\ell_1} penalty and an additional penalty
//' @description ADMM algorithm for high dimensional Quadratic regression with \eqn{\ell_1} and another penalty:
//' \deqn{\argmin_Y \frac{1}{2n}\sum_{i=1}^n (Y_i-X_i^\top \Omega X_i)^2+\lambda_1 \|\Omega\|_1+\lambda_2 p(\Omega)+\lambda_2 p(\Omega^\top).} 
//' @param \code{X}, a \eqn{n\times p} data matrix.
//' @param \code{Y}, a \eqn{n} dimensional response vector.  
//' @param \code{lambda1}, user supplied tuning parameters.
//' @param \code{lambda2}, user supplied tuning parameters.
//' @param \code{type} The additional penalty to use for the quadratic regression.
//'  \eqn{2} is the \eqn{\ell_\infty} penalty; \eqn{3} is the \eqn{\ell_2} penalty (Group LASSO); \eqn{4} is the hybrid \eqn{\ell_1/\ell_\infty} penalty.
//' @param \code{err_abs}, \code{err_rel},   the precision used to stop the convergence of ADMM. 
//' @param \code{maxIter}, maximum number of iterations. Default is 200.
//' @param \code{rho}, initial step parameter for ADMM.
//' 
//' @return A list with components
//' \item{Omega}{a list of sparse \eqn{p \times p} matrices corresponding to \code{lambda1}[k] and \code{lambda2}[k].}
//' \item{lambda1}{the used lambda1 for the solution path.}
//' \item{lambda2}{the used lambda2 for the solution path.}
//' \item{niter}{the number of iterations.}
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List qr3(arma::mat X,arma::vec Y,arma::vec lambda1,arma::vec lambda2,int type=2,double err_abs=10^(-4),double err_rel=10^(-3),int maxIter=200,double rho=5){
  int n=X.n_rows;
  int p=X.n_cols;
  int nlambda=lambda1.size();
  /*Preparing*/
  arma::mat D0=gram(X,Y/n);
  arma::mat Dk=D0;
  arma::mat H0=X*X.t();
  arma::mat H2=H0%H0/n;
  arma::mat H=inv_sympd(H2+rho*arma::eye(n,n));
  
  Rcpp::List Omega_all(nlambda);
  arma::vec niter=arma::zeros(nlambda);
  arma::vec rholist=arma::zeros(nlambda);

  arma::vec w=Y;
  double lam1;
  double lam2;
  /*Intialization*/
  arma::mat B1;
  arma::mat B2;
  arma::mat B3;
  arma::mat B4;
  arma::mat old_B;
  arma::mat B=arma::zeros(p,p);
  arma::mat U1=arma::zeros(p,p);
  arma::mat U2=arma::zeros(p,p);
  arma::mat U3=arma::zeros(p,p);
  arma::mat U4=arma::zeros(p,p);
  
  for (int k=0;k<nlambda;++k) {
    lam1=lambda1(k);
    lam2=lambda2(k);
    double ee_pri=1;
    double ee_dual=1;
    int i=0;
    while ((i<maxIter)||(i==0))
    {
      /* Dual-update*/
      Dk=D0+rho*(B-U1);
      w=H*qrow(X,Dk/n);
      B1=Dk/rho-gram(X,w/rho);
      B2=soft(B-U2,lam1/rho);
      B3=proc(B-U3,lam2/rho,type,1);
      B4=proc((B-U4).t(),lam2/rho,type,1);B4=B4.t();
      /*Prime-update*/
      old_B=B;
      B=(B1+B2+B3+B4)/4;
      /*U-update*/
      U1=U1+B1-B;
      U2=U2+B2-B;
      U3=U3+B3-B;
      U4=U4+B4-B;
      i=i+1;
      
      /*Stop Rule*/
      ee_dual=rho*norm(B-old_B,"fro"); /*dual residual*/
      ee_pri=norm(B1-B,"fro")+norm(B2-B,"fro")+norm(B3-B,"fro")+norm(B4-B,"fro"); /*primal residual*/
      double err_pri_new=p*err_abs+err_rel*norm(B,"fro");
      double err_dual_new=p*err_abs+err_rel*norm(U1,"fro");
      if ((ee_dual<err_dual_new)&&(ee_pri<err_pri_new)) break;
      /*Varying Penalty Parameter*/
      /* if (ee_pri>10*ee_dual) {rho=2*rho;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
       /*if (ee_dual>10*ee_pri) {rho=rho/2;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
    }
    niter(k)=i;
    rholist(k)=rho;
    Omega_all[k]=arma::sp_mat(B2);
  }
  return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                            Rcpp::Named("lambda1") =lambda1,
                            Rcpp::Named("lambda2") =lambda2,
                            Rcpp::Named("rho") =rholist,
                            Rcpp::Named("niter") =niter); }




//' @title Quadratic regression with nuclear norm penalty 
//' @description ADMM algorithm for high dimensional Quadratic regression with a \eqn{\ell_1} norm penalty:
//' \deqn{\argmin_Y \frac{1}{2n}\sum_{i=1}^n (Y_i-X_i^\top \Omega X_i)^2+\lambda \|\Omega\|_*.} 
//' @param \code{X}, a \eqn{n*p} input data matrix.
//' @param \code{Y}, a \eqn{n} response vector.  
//' @param \code{lambda}, user supplied tuning parameter; 
//' @param \code{err_abs}, \code{err_rel}   the precision used to stop the convergence of ADMM. 
//' @param \code{maxIter}, Maximum number of iterations. Default is 1000.
//' @param \code{rho}, initial step parameter for ADMM.
//' @return A list with components
//' \item{Omega}{a list of sparse p*p matrices corresponding to lambda.}
//' \item{lambda}{the used lambda for the solution path.}
//' \item{niter}{the number of iterations for each element of lambda.}
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List qr1_rank(arma::mat X,arma::vec Y,arma::vec lambda,double err_abs=10^(-4),double err_rel=10^(-3),int maxIter=200,double rho=1){
         int n=X.n_rows;
         int p=X.n_cols;
         lambda=arma::sort(lambda,"descend");
         int nlambda=lambda.size();
         /*Centering int*/
         arma::mat D0=gram(X,Y/n);
         arma::mat Dk=D0;
         arma::mat H0=X*X.t();
         arma::mat H2=H0%H0/n;
         arma::mat H=inv_sympd(H2+rho*arma::eye(n,n));
         Rcpp::List Omega_all(nlambda);
         arma::vec niter=lambda;
         arma::vec rholist=lambda;
         /*Intialization*/
         arma::mat A=arma::zeros(p,p);
         arma::mat U=arma::zeros(p,p);
         arma::mat B=A;
         arma::mat old_B=B;
         arma::vec w=Y;
         double lam;
         double ee_pri=1;
         double ee_dual=1;
         
         for (int k=0;k<nlambda;++k) {
           lam=lambda(k);
           int i=0;
           while ((i<maxIter)||(i==0))
           {
             Dk=D0+rho*(B-U);
             /* A-update*/
             w=H*qrow(X,Dk/n);
             A=(Dk-gram(X,w))/rho;
             /*B-update*/
             old_B=B;
             B=nuclear_mat(A+U,lam/rho);
             /*U-update*/
             U=A-B+U;
             i=i+1;
             /*Stop Rule*/
             ee_dual=rho*norm(B-old_B,"fro"); /*dual residual*/
             ee_pri=norm(A-B,"fro"); /*primal residual*/
             double err_pri_new=p*err_abs+err_rel*std::max(norm(A,"fro"),norm(B,"fro"));
             double err_dual_new=p*err_abs+err_rel*norm(U,"fro");
             if ((ee_dual<err_dual_new)&&(ee_pri<err_pri_new)) break;
             /*Varying Penalty Parameter*/
             /*if (ee_pri>10*ee_dual) {rho=2*rho;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
             /*if (ee_dual>10*ee_pri) {rho=rho/2;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
           }
           Omega_all(k)=arma::sp_mat(B);
           niter(k)=i;
           rholist(k)=rho;
         }
         return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                                   Rcpp::Named("lambda") =lambda,
                                   Rcpp::Named("rho") =rholist,
                                   Rcpp::Named("niter") =niter); }

//' @title Quadratic regression with \eqn{\ell_1} penalty and nuclear norm penalty (Sparse and Low rank)
//' @description ADMM algorithm for high dimensional Quadratic regression with \eqn{\ell_1} and nuclear norm penalty:
//' \deqn{\argmin_Y \frac{1}{2n}\sum_{i=1}^n (Y_i-X_i^\top \Omega X_i)^2+\lambda_1 \|\Omega\|_1+\lambda_2 \|\Omega\|_*.} 
//' @param \code{X}, a \eqn{n\times p} data matrix.
//' @param \code{Y}, a \eqn{n} dimensional response vector.  
//' @param \code{lambda1}, user supplied tuning parameters.
//' @param \code{lambda2}, user supplied tuning parameters.
//' @param \code{err_abs}, \code{err_rel},   the precision used to stop the convergence of ADMM. 
//' @param \code{maxIter}, maximum number of iterations. Default is 200.
//' @param \code{rho}, initial step parameter for ADMM.
//' 
//' @return A list with components
//' \item{Omega}{a list of sparse \eqn{p \times p} matrices corresponding to lambda.}
//' \item{lambda1}{the used lambda1 for the solution path.}
//' \item{lambda2}{the used lambda2 for the solution path.}
//' \item{niter}{the number of iterations for each element of (lambda1,lambda2).}
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List qr3_rank(arma::mat X,arma::vec Y,arma::vec lambda1,arma::vec lambda2,double err_abs=10^(-4),double err_rel=10^(-3),int maxIter=200,double rho=5){
  int n=X.n_rows;
  int p=X.n_cols;
  int nlambda=lambda1.size();
    /*Preparing*/
  arma::mat D0=gram(X,Y/n);
  arma::mat Dk=D0;
  arma::mat H0=X*X.t();
  arma::mat H2=H0%H0/n;
  arma::mat H=inv_sympd(H2+rho*arma::eye(n,n));
  Rcpp::List Omega_all(nlambda);
  arma::vec niter=arma::zeros(nlambda);
  arma::vec rholist=arma::zeros(nlambda);
  arma::vec w=Y;
  double lam1;
  double lam2;
  /*Intialization*/
  arma::mat B1;
  arma::mat B2;
  arma::mat B3;
  arma::mat old_B;
  arma::mat B=arma::zeros(p,p);
  arma::mat U1=arma::zeros(p,p);
  arma::mat U2=arma::zeros(p,p);
  arma::mat U3=arma::zeros(p,p);
  
  for (int k=0;k<nlambda;++k) {
    lam1=lambda1(k);
    lam2=lambda2(k);
    double ee_pri=1;
    double ee_dual=1;
    int i=0;
  while ((i<maxIter)||(i==0))
{
 /* Dual-update*/
    Dk=D0+rho*(B-U1);
    w=H*qrow(X,Dk/n);
    B1=Dk/rho-gram(X,w/rho);
    B2=soft(B-U2,lam1/rho);
    B3=nuclear_mat(B-U3,lam2/rho);
/*Prime-update*/
    old_B=B;
    B=(B1+B2+B3)/3;
/*U-update*/
    U1=U1+B1-B;
    U2=U2+B2-B;
    U3=U3+B3-B;
    i=i+1;
/*Stop Rule*/
    ee_dual=rho*norm(B-old_B,"fro"); /*dual residual*/
    ee_pri=norm(B1-B,"fro")+norm(B2-B,"fro")+norm(B3-B,"fro"); /*primal residual*/
    double err_pri_new=p*err_abs+err_rel*norm(B,"fro");
    double err_dual_new=p*err_abs+err_rel*norm(U1,"fro");
    if ((ee_dual<err_dual_new)&&(ee_pri<err_pri_new)) break;
/*Varying Penalty Parameter*/
/* if (ee_pri>10*ee_dual) {rho=2*rho;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
/* if (ee_dual>10*ee_pri) {rho=rho/2;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
  }
                    niter(k)=i;
                    rholist(k)=rho;
                    Omega_all[k]=arma::sp_mat(B2);
                }
                return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                                          Rcpp::Named("lambda1") =lambda1,
                                          Rcpp::Named("lambda2") =lambda2,
                                          Rcpp::Named("rho") =rholist,
                                          Rcpp::Named("niter") =niter); }