# Introduction for R package HiQR
We develop an Efficient admm algorithm via the QUAdratic Loss (HiQR) for quadratic regression with penalty. The computation complexity for each iteration of the algorithm is linear in both the sample size (n) and the number of parameters (p^2).  


This is my first R package and welcome any comments or suggestions.

# References 
Still in develop

# Getting Started
These instructions will give you a toy example for implementing the package.

## Prerequisites
What things you need to install the software and how to install them.  The key functions of the package is writing in C++ supported by the great Rcpp package. So, make sure your OS can complies C++ code. For example,  you should install Rtools under Windows and Xcode under MacOS.  After that, the following R packages are also necessary.

```
install.packages("MASS")
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("Matrix")
install.packages("devtools")
```
## Install 

```
library("devtools")
devtools::install_github("cescwang85/HiQR")
```

## Toy example 

```
rm(list = ls())
set.seed(123)
library('MASS')
library('Rcpp')
library('Matrix')
library('HiQR')
n=200
p=100
Omega<-toeplitz(0.5^(1:p-1))
X=mvrnorm(n,rep(0,p),solve(Omega))
```

