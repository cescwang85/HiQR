% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{l1inf_vec}
\alias{l1inf_vec}
\title{Proximal projection of a hybrid \eqn{\ell_1/\ell_\infty} penalty for a vector}
\usage{
l1inf_vec(x, lambda = 0)
}
\arguments{
\item{\code{x}, }{a vector.}

\item{\code{lambda}, }{a scalar.}
}
\value{
a vector.
}
\description{
Proximal projection solution of a vector which is the solution:
\deqn{\argmin_y ||y-x||_2^2/2+\lambda  \max\{|y_1|, \sum_{i=2}^p |y_i|\}.}
}
