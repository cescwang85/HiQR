#' critical success index (CSI) or threat score (TS)
#' @param \code{A}, a matrix(True parameter).
#' @param \code{B}, a matrix(Estimation).
#' @return critical success index (CSI).
csi<-function(A,B)
{A=as.matrix(A);B=as.matrix(B);
A[1,1]=0; B[1,1]=0; ##ignore the constant
size=nnzero(A*B)/(nnzero(A)+nnzero(B)-nnzero(A*B))
return(size)
}