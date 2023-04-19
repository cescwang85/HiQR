#' Transform a long vector into a matrix 
#' @param \code{A}, a matrix(True parameter).
#' @param \code{B}, a matrix(Estimation).
#' @return a 4 dimensional vector with nnzero of \code{B}, Sensitivity, Specificity and also the Frobenious norm of A-B.
measure<-function(A,B)
{
  B1=B[A==0];
  B2=B[A!=0]
  sens=mean(B2!=0)
  spec=mean(B1==0)
  size=nnzero(B)
  return(c(size,sens,spec,norm(A-B,type='F')))
}