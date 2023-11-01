rm(list=ls())
library(HiQR)
setwd('~/Dropbox/Share/HQR_ADMM/sim/')
source('ef_1.R')
source('ADMM_0.1.19.R')
#the l_infinity prox function.
x<-rnorm(1000)
lambda=sum(abs(x))/2

y1=as.vector(inf_vec(x,lambda))
y2=l_inf_prox(x,lambda)
max(abs(y1-y2))

##Soft

y3=soft_vec(x,0.3)
y4=soft_vec(x,0.3)
max(abs(y3-y4))