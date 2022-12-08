##check the HiQR package with L1 penalty
#library("devtools")
#devtools::install_github("cescwang85/PIE")


rm(list=ls())
library(HiQR)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
set.seed(99)
library('tictoc')
library('glmnet')
library('Matrix')
library('FAMILY')
library('PIE')

p=50
n=200;
Omega=matrix(0,nrow=p,ncol=p);
Omega[6,6]=1
Omega[1,6]=2;Omega[6,10]=2;Omega=(Omega+t(Omega))/2;
beta<-rep(0,p);beta[c(1,6,10)]=1;
X=matrix(rnorm(n*p),nrow =n);
Y=as.vector(diag(X%*%Omega%*%t(X))+X%*%beta+rnorm(n));
## PIE basd on the response
yOme<-PIE(X,Y)
##PIE based on the residual
#all-pairs LASSO
tic('HiQR')
obj0=hiqr(X,Y)
toc()


##LASSO+l_infty penalty
tic('HiQR')
obj1=hiqr(X,Y,nlambda1 =25,nlambda2 =25, type =1)
toc()

##LASSO+group LASSO
tic('HiQR')
obj2=hiqr(X,Y,nlambda1 =25,nlambda2 =25,type =2)
toc()

##LASSO+hybrid \eqn{\ell_1/\ell_\infty} penalty
tic('HiQR')
obj3=hiqr(X,Y,nlambda1 =25,nlambda2 =25,type =3)
toc()

obj3$Omega[[420]][1:10,1:10]
Omega[1:10,1:10]


library(FAMILY)
library(pROC)
library(pheatmap)

#####################################################################################
#####################################################################################
############################# EXAMPLE - CONTINUOUS RESPONSE #########################
#####################################################################################
#####################################################################################

############################## GENERATE DATA ########################################

#Generate training set of covariates X and Z
p=30
set.seed(1)
X.tr<- matrix(rnorm(10*100),ncol =p, nrow = 100)
#Scale appropiately
X.tr<- scale(X.tr, scale = FALSE)


#Generate full matrix of Covariates
w.tr<- c()
X1<- cbind(1,X.tr)
for(i in 1:16){
  for(j in 1:11){
    w.tr<- cbind(w.tr,X1[,j]*X1[,i])
  }
}
B<- matrix(0,ncol = 16,nrow = 11)
rownames(B)<- c("inter" , paste("X",1:(nrow(B)-1),sep = ""))
colnames(B)<- c("inter" , paste("X",1:(ncol(B)-1),sep = ""))

# First, we simulate data as follows:
# The first five features in X, and the first five features in Z, are non-zero.
# And given the non-zero main effects, all possible interactions are involved.
# We call this "high strong heredity"
B_high_SH<- B
B_high_SH[1:6,1:6]<- 1
#View true coefficient matrix
pheatmap(as.matrix(B_high_SH), scale="none",
         cluster_rows=FALSE, cluster_cols=FALSE)

Y_high_SH <- as.vector(w.tr%*%as.vector(B_high_SH))+rnorm(100,sd = 2)


############################## FIT SOME MODELS ########################################

#Define alphas and lambdas
#Define 3 different alpha values
#Low alpha values penalize groups more
#High alpha values penalize individual Interactions more
alphas<- c(0.01,0.5,0.99)
lambdas<- seq(0.1,1,length = 50)

#high Strong heredity with l2 norm
fit_high_SH<- FAMILY(X.tr, X.tr, Y_high_SH, lambdas ,
                     alphas, quad = TRUE,iter=500, verbose = TRUE )

#Low alpha, higher penalty on groups
plot(fit_high_SH$Estimate[[ 1 ]][[ 25 ]])
#Medium alpha, equal penalty on groups and individual interactions
#plot(fit_high_SH$Estimate[[ 2 ]][[ 25  ]])
#High alpha, more penalty on individual interactions
#plot(fit_high_SH$Estimate[[ 3 ]][[ 40 ]])

