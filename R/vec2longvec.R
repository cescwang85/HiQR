#' Calculate a interaction vector without duplication for a vector 
#' @param \code{a}, a \eqn{p} dimensional vector.
#' @return a \eqn{p(p+3)/2} dimensional vector
vec2longvec<-function(a)
{a<-as.vector(a);
p<-length(a);
id<-as.vector(tril(matrix(1:p^2,p))); ##pp=p*(p+3)/2
id<-id[id>0];
A<-a%*%t(a);
return(c(a,as.vector(A[id])))}