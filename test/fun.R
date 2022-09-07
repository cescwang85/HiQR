##Given covariate a, generate a vector of c(a,a_1a_1,a_1a_2,\cdots,a_pa_p)
ww<-function(a)
{a<-as.vector(a);
m<-length(a);
id<-as.vector(tril(matrix(1:m^2,m))); ##pp=p*(p+3)/2
id<-id[id>0];
A<-a%*%t(a);
return(c(a,as.vector(A[id])))}

inww<-function(aa)
{pp=length(aa);
p=(sqrt(8*pp+1)-1)/2;  ## p(p+1)/2=pp
p=p-1
Omega<-matrix(0,p,p);
id<-as.vector(tril(matrix(1:p^2,p)));
id<-id[id>0];
Omega[id]<-aa[-(1:(p+1))];
Omega1<-matrix(0,p+1,p+1);
Omega1[2:(p+1),2:(p+1)]=Omega;
Omega1[,1]=aa[1:(p+1)]
Omega1=Matrix((Omega1+t(Omega1))/2, sparse = TRUE) 
return(Omega1)
}

ww_new<-function(a)
{a<-as.vector(a);
return(c(as.vector(a%*%t(a))))}