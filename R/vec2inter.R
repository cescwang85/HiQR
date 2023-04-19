#' Calculate a p2 interaction vector for a vector 
#' @param \code{a}, a \eqn{p} dimensional vector.
#' @return a \eqn{p^2} dimensional vector
vec2inter<-function(a)
{a<-as.vector(a);
return(c(as.vector(a%*%t(a))))}