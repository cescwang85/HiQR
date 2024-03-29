% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{soft_vec}
\alias{soft_vec}
\title{Element-wise soft-thresholding function for a vector}
\usage{
soft_vec(x, lambda = 0)
}
\arguments{
\item{\code{x}, }{a vector.}

\item{\code{lambda}, }{a scalar.}
}
\value{
a vector after thresholding.
}
\description{
Soft-thresholding function for a vector which is the solution of
\deqn{\argmin_y ||y-x||_2^2/2+\lambda ||y||_1.}
}
