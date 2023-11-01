##simulation for all-pairs lasso
rm(list=ls())
library(HiQR)
setwd('~/Dropbox/Share/HQR_ADMM/sim')
set.seed(99)
library('glmnet')
library('Matrix')
library('MASS')
library("lattice")


ttime<-function(obj)
{t1<-proc.time();
obj
t2<-as.vector(proc.time()-t1);
return(sum(t2[1:2]))}

n=500
## Generate data
p=200
## Parameters
Sigma=toeplitz(0.5^(0:(p-1)))
Omega3<-matrix(0,nrow=p+1,ncol=p+1)
Omega3[2,6]=3;Omega3[6,6]=-2.5;Omega3[6,11]=4
Omega3=(Omega3+t(Omega3))/2;
Omega2=Omega3;
Omega2[1,6]=-1;Omega2[6,1]=-1;
Omega1=Omega2;
Omega1[1,2]=1;Omega1[2,1]=1
Omega1[1,11]=1;Omega1[11,1]=1
Omega3[c(1,2,6,11),c(1,2,6,11)] ##check matrix
Omegalist=list(Omega1,Omega2,Omega3)

##
for (mk in 1:3) {
Omega=Omegalist[[mk]]

X0<-mvrnorm(n,rep(0,p),Sigma)
X<-cbind(1,X0)
Y=as.vector(diag(X%*%Omega%*%t(X))+rnorm(n))


nalpha=50
runtime<-rep(0,7);
## type=1 (all-pairs LASSO)
runtime[1]=ttime({obj1=hiqr(X0,Y,type=1)})
## type=2 (Ridge regression)
runtime[2]=ttime({obj2=hiqr(X0,Y,type=2)})
## type=5 (Reduced Rank)
runtime[3]=ttime({obj5=hiqr(X0,Y,type=5)})
## type=12 (l1+l_inf)
runtime[4]=ttime({obj12=hiqr(X0,Y,type=12,nalpha=nalpha)})
## type=13 (l1+l2)
runtime[5]=ttime({obj13=hiqr(X0,Y,type=13,nalpha=nalpha)})
## type=14 (l1+l1/l_inf)
runtime[6]=ttime({obj14=hiqr(X0,Y,type=14,nalpha=nalpha)})
## type=15 (l1+nuclear norm)
runtime[7]=ttime({obj15=hiqr(X0,Y,type=15,nalpha=nalpha)})
runtime

lambda<-exp(seq(from=0,to=log(0.25),length.out=50));
alpha<-seq(from=0.1,to=1,length.out=nalpha);

nlambda<-length(lambda)*length(alpha)
re<-matrix(0,nrow=4,ncol=nlambda)
for (k in 1:nlambda)
{A1=obj12$Omega[[k]];
A2=obj13$Omega[[k]];
A3=obj14$Omega[[k]];
A4=obj15$Omega[[k]];

re[1,k]=csi(Omega,A1);
re[2,k]=csi(Omega,A2);
re[3,k]=csi(Omega,A3);
re[4,k]=csi(Omega,A4);
}
grid <- expand.grid(x=lambda, y=alpha)

grid$z <- re[1,]
pdf(paste(n,'-',p,'-',mk,'-','12.pdf',sep=''))
plot1<-levelplot(z ~ x * y, grid, xlab="lambda", ylab="alpha",col.regions = gray(100:0/100))
print(plot1)
dev.off()

grid$z <- re[2,]
pdf(paste(n,'-',p,'-',mk,'-','13.pdf',sep=''))
plot2<-levelplot(z ~ x * y, grid, xlab="lambda", ylab="alpha",col.regions = gray(100:0/100))
print(plot2)
dev.off()

grid$z <- re[3,]
pdf(paste(n,'-',p,'-',mk,'-','14.pdf',sep=''))
plot3<-levelplot(z ~ x * y, grid, xlab="lambda", ylab="alpha",col.regions = gray(100:0/100))
print(plot3)
dev.off()

grid$z <- re[4,]
pdf(paste(n,'-',p,'-',mk,'-','15.pdf',sep=''))
plot4<-levelplot(z ~ x * y, grid, xlab="lambda", ylab="alpha",col.regions = gray(100:0/100))
print(plot4)
dev.off()
}