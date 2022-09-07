//This is the main function for the penalized quadratic regression with \ell_1 penalty 

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Element-wise Soft Thresholding Function for Vector
//' @description Soft thresholding function for a vector which is the solution of the  
//' \deqn{\argmin_x ||y-x||_2^2/2+\lambda ||x||_1.}
//' @param \code{y} a vector
//' @param \code{lambda} a scalar
//' @return Vector after threholding
//' 
//' @export
//'
// [[Rcpp::export]]
arma::vec soft_vec(arma::vec y,double lambda){
  lambda=std::abs(lambda);
  arma::vec x=(y>=lambda)%(y-lambda)+(y<=(-lambda))%(y+lambda); 
  return x;}

//' @title Element-wise Soft Thresholding Function for a Matrix
//' @description Soft thresholding function for Matrix which is the solution of the  
//' \deqn{\argmin_X ||X-A||_2^2/2+a ||X||_1.}
//' @param \code{A} Matrix
//' @param \code{a} scalar
//' @return Matrix after threholding
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat soft(arma::mat A,double a){
  a=std::abs(a);
  arma::mat B=(A>=a)%(A-a)+(A<=(-a))%(A+a); 
  B(0,0)=A(0,0); //Not penalized for the constant
  return B;}



//' @title Proximal projection of a Vector with  \eqn{\ell_\infty} penalty. 
//' @description Proximal projection solution of a vector which is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_\infty.}
//' @param \code{x} a vector
//' @param \code{lambda} a scalar
//' @return A vector which is proximal projection solution.
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

//' @title Proximal projection of a Matrix with  \eqn{\ell_\infty} penalty. 
//' @description Proximal projection solution of a Matrix where each column is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_\infty.}
//' @param \code{X} a Matrix
//' @param \code{lambda} a scalar
//' @param \code{k} Not penalized the first \code{k} columns. Default is 1.
//' @return A Matrix.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat inf_mat(arma::mat X,double lambda=0,int k=1){
  arma::mat Y=X;
  int p=X.n_cols;
  for(arma::uword i=k; i<p; ++i) {Y.col(i)=inf_vec(X.col(i),lambda);}
  return(Y);}


//' @title Proximal projection of a Vector with  \eqn{\ell_2} penalty (Group LASSO). 
//' @description Proximal projection solution of a vector which is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_2.}
//' @param \code{x} a vector
//' @param \code{lambda} a scalar
//' @return A vector which is proximal projection solution.
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

//' @title Proximal projection of a Matrix with  \eqn{\ell_2} penalty. 
//' @description Proximal projection solution of a Matrix where each column is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_2.}
//' @param \code{X} a Matrix
//' @param \code{lambda} a scalar
//' @param \code{k} Not penalized the first \code{k} columns. Default is 1.
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




//' @title Proximal projection of a Vector with a hybrid \eqn{\ell_1/\ell_\infty} penalty. 
//' @description Proximal projection solution of a vector which is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda  \max\{|y_1|, \sum_{i=2}^p |y_i|\}.}
//' @param \code{x} a vector
//' @param \code{lambda} a scalar
//' @return A vector which is proximal projection solution.
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


//' @title Proximal projection of a Matrix with a hybrid \eqn{\ell_1/\ell_\infty} penalty. 
//' @description Proximal projection solution of a Matrix where each column is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda  \max\{|y_1|, \sum_{i=2}^p |y_i|\}.}
//' @param \code{X} a Matrix
//' @param \code{lambda} a scalar
//' @param \code{k} Not penalized the first \code{k} columns. Default is 1.
//' @return A Matrix.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat l1inf_mat(arma::mat X,double lambda=0,int k=1){
  arma::mat Y=X;
  int p=X.n_cols;
  for(arma::uword i=k; i<p; ++i) {Y.col(i)=l1inf_vec(X.col(i),lambda);}
  return(Y);}

//' @title Proximal projection of a Matrix with different penalties. 
//' @description Proximal projection solution of a Matrix where each column is the solution:
//' \deqn{\argmin_y ||y-x||_2^2/2+\lambda \cdot p(y).}
//' @param \code{X} a Matrix
//' @param \code{lambda} a scalar
//' @param \code{type} The penalty to use for the columns of matrix \code{X}.
//'  \eqn{1} is the \eqn{\ell_\infty} penalty; \eqn{2} is the \eqn{\ell_2} penalty (Group LASSO); \eqn{3} is the hybrid \eqn{\ell_1/\ell_\infty} penalty.
//' @param \code{k} Not penalized the first \code{k} columns. Default is 1.
//' @return A Matrix.
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat proc(arma::mat X,double lambda=0,int type=1,int k=1){
  arma::mat Y=X;
  if (type==1) Y=inf_mat(X,lambda,k);
  if (type==2) Y=l2_mat(X,lambda,k);
  if (type==3) Y=l1inf_mat(X,lambda,k);
  return(Y);}


//' @title Weighted Gram Matrix
//' @description \code{gram} computes the weighted Gram matrix: 
//' \deqn{ X^\top diag(w) X.}
//' @param \code{X}  a \eqn{n \times p} Matrix
//' @param \code{w} a \eqn{n} dimensional weight vector
//' @return Weighted Gram Matrix
//' 
//' @export
// [[Rcpp::export]]
arma::mat gram(arma::mat X,arma::vec w){
  arma::mat B=X.each_col()%w;
  B=X.t()*B;
  return B;}

//' @title Quadratic of Each Row 
//' @description \code{qrow} computes the quadratic of Each Row : 
//' \deqn{ \diag(X W X^\top)=(X_1^\top W X_1,\cdots,X_n^\top W X_n)^\top}
//' @param \code{X}  a \eqn{n*p} Matrix
//' @param \code{W} a \eqn{p*p} weight Matrix
//' @return Quadratic of Each Row 
//' 
//' @export
// [[Rcpp::export]]
arma::vec qrow(arma::mat X,arma::mat W){
  arma::vec B=arma::sum(((X*W)%X),1);
  return B;}



//' @title Quadratic regression with squared \eqn{\ell_2} penalty (Ridge regression)  
//' @description Algorithm for high dimensional Quadratic regression with squared \eqn{\ell_2} penalty
//' @param \code{X} a \eqn{n*p} input data matrix.
//' @param \code{Y} a \eqn{n} response vector.  
//' @param \code{lambda} user supplied tuning parameter; 
//' @return A list with components
//' \item{Omega}{a list of sparse p*p matrices corresponding to lambda.}
//' \item{lambda}{the used lambda.}
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List ridge_qr(arma::mat X,arma::vec Y,arma::vec lambda){
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
//' @description ADMM algorithm for high dimensional Quadratic regression with \eqn{\ell_1} penalty
//' @param \code{X} a \eqn{n*p} input data matrix.
//' @param \code{Y} a \eqn{n} response vector.  
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
if (ee_pri>10*ee_dual) {rho=2*rho;H=inv_sympd(H2+rho*arma::eye(n,n));}
if (ee_dual>10*ee_pri) {rho=rho/2;H=inv_sympd(H2+rho*arma::eye(n,n));}
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
//' @description ADMM algorithm for high dimensional Quadratic regression with \eqn{\ell_1} and other penalties.
//' @param \code{X}, a \eqn{n*p} input data matrix.
//' @param \code{Y}, a \eqn{n} response vector.  
//' @param \code{lambda}, user supplied tuning parameter; 
//' @param \code{type} The additional penalty to use for the quadratic regression.
//'  \eqn{1} is the \eqn{\ell_\infty} penalty; \eqn{2} is the \eqn{\ell_2} penalty (Group LASSO); \eqn{3} is the hybrid \eqn{\ell_1/\ell_\infty} penalty.
//' @param \code{err_abs}, \code{err_rel}   the precision used to stop the convergence of ADMM. 
//' @param \code{maxIter}, maximum number of iterations. Default is 200.
//' @param \code{rho}, initial step parameter for ADMM.
//' 
//' @return A list with components
//' \item{Omega}{a list of sparse p*p matrices corresponding to lambda.}
//' \item{lambda}{the used lambda for the solution path.}
//' \item{niter}{the number of iterations for each element of lambda.}
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List qr2(arma::mat X,arma::vec Y,arma::vec lambda1,arma::vec lambda2,int type=1,double err_abs=10^(-4),double err_rel=10^(-3),int maxIter=200,double rho=5){
  int n=X.n_rows;
  int p=X.n_cols;
  lambda1=arma::sort(lambda1,"descend");
  lambda2=arma::sort(lambda2,"descend");
  int nlambda1=lambda1.size();
  int nlambda2=lambda2.size();
  /*Preparing*/
  arma::mat D0=gram(X,Y/n);
  arma::mat Dk=D0;
  arma::mat H0=X*X.t();
  arma::mat H2=H0%H0/n;
  arma::mat H=inv_sympd(H2+rho*arma::eye(n,n));
  
  Rcpp::List Omega_all(nlambda1*nlambda2);
  arma::mat niter=arma::zeros(nlambda1,nlambda2);
  arma::mat rholist=arma::zeros(nlambda1,nlambda2);

  arma::vec w=Y;
  double lam1;
  double lam2;
  /*Intialization*/
  arma::mat B1;
  arma::mat B2;
  arma::mat B3;
  arma::mat B4;
  arma::mat old_B;


  Rcpp::List Omega_part(nlambda2);
  Omega_part(0)=arma::zeros(p,p);
  for (int k1=0;k1<nlambda1;++k1) {
    lam1=lambda1(k1);
    arma::mat B=Omega_part(0);
    arma::mat U1=arma::zeros(p,p);
    arma::mat U2=U1;
    arma::mat U3=U1;
    arma::mat U4=U1;
    
    double ee_pri=1;
    double ee_dual=1;
    for (int k2=0;k2<nlambda2;++k2) {
      lam2=lambda2(k2);
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
      if (ee_pri>10*ee_dual) {rho=2*rho;H=inv_sympd(H2+rho*arma::eye(n,n));}
      if (ee_dual>10*ee_pri) {rho=rho/2;H=inv_sympd(H2+rho*arma::eye(n,n));}
    }
    Omega_part[k2]=B2;
    niter(k1,k2)=i;
    rholist(k1,k2)=rho;
    Omega_all[k1*nlambda2+k2]=arma::sp_mat(B2);
  }
  }
  return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                            Rcpp::Named("lambda1") =lambda1,
                            Rcpp::Named("lambda2") =lambda2,
                            Rcpp::Named("rho") =rholist,
                            Rcpp::Named("niter") =niter); }



