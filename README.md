# Introduction for R package HiQR
We develop an efficient algorithm for High dimensional Quadratic Regression (HiQR):
$$\arg\min \frac{1}{2n} \sum_{i=1}^n (y_i-x_i^{T} B x_i)^2+f(B).$$

 Several penalty functions are included: 
1. (type=1) LASSO ($\ell_1$ norm): $f(B)=\|B\|_1$;
2. (type=2) Ridge regression (squared $\ell_2$ norm): $f(B)=\|B\|_2^2$;
3. (type=5) Reduced rank regression (nuclear norm): $f(B)=\|B\|_*$;
4. (type=12) $\ell_1+\ell_2$ norm: 
$$f(B)=\lambda_1 \|B\|_1+\lambda_2 \sum_{k=2}^p \|B_{\cdot,k}\|_2+\lambda_2 \sum_{k=2}^p \|B_{k,\cdot}\|_2;$$
5. (type=13) $\ell_1+\ell_\infty$ norm:
$$f(B)=\lambda_1 \|B\|_1+\lambda_2 \sum_{k=2}^p \|B_{\cdot,k}\|_\infty+\lambda_2 \sum_{k=2}^p \|B_{k,\cdot}\|_\infty;$$
6. (type=14) $\ell_1+\ell_1/\ell_\infty$ norm: 
$$f(B)=\lambda_1 \|B\|_1+\lambda_2 \sum_{k=2}^p \max\{|B_{1,k}|, \|B_{-1,k}\|_1\}+\lambda_2 \sum_{k=2}^p \max\{|B_{k,1}|, \|B_{k,-1}\|_1\};$$
7. (type=15) $\ell_1+\ell_*$ norm: 
$$f(B)=\lambda_1 \|B\|_1+\lambda_2 \|B\|_*.$$


# References 
Cheng Wang, Haozhe Chen and Binyan Jiang. "HiQR: An efficient algorithm for high dimensional quadratic regression with penalties", 2023.  

## Prerequisites
What things you need to install the software and how to install them.  The key functions of the package is writing in C++ supported by the great Rcpp package. So, make sure your OS can complies C++ code. For example,  you should install Rtools under Windows and Xcode under MacOS.  After that, the following R packages are also necessary.

```
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
X=matrix(rnorm(n*p),nrow=n)
Y=rnorm(n)


## LASSO
obj2<-hiqr(X,Y,type=1)
## Ridege 
obj2<-hiqr(X,Y,type=2)
## Ridege 
obj5<-hiqr(X,Y,type=5)
## $\ell_1+\ell_2$ norm:
obj12<-hiqr(X,Y,type=12)
## $\ell_1+\ell_\infty$ norm:
obj13<-hiqr(X,Y,type=13)
## $\ell_1+\ell_1/\ell_\infty$ norm: 
obj14<-hiqr(X,Y,type=14)
## $\ell_1+\ell_*$ norm: 
obj15<-hiqr(X,Y,type=15)
```

