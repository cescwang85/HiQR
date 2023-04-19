#' Transform a long vector into a matrix 
#' @param \code{aa}, a \eqn{p(p+3)/2} dimensional vector.
#' @return a \eqn{(p+1)\times (p+1)} sparse matrix.
longvec2mat<-function(aa)
{pp=length(aa);
p=(sqrt(8*pp+1)-1)/2;  ## p(p+1)/2=pp
p=p-1
Omega<-matrix(0,p,p);
id<-as.vector(tril(matrix(1:p^2,p)));
id<-id[id>0];
Omega[id]<-aa[-(1:(p+1))];
Omega1<-matrix(0,p+1,p+1);
Omega1[2:(p+1),2:(p+1)]=Omega;
Omega1[,1]=aa[1:(p+1)]
Omega1=Matrix((Omega1+t(Omega1))/2, sparse = TRUE) 
return(Omega1)
}