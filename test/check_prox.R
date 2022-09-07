rm(list=ls())
library(HiQR)



#The function I use for the l_infinity prox function.

l_inf_prox<- function(y,lam,na.rm = TRUE){
  had.na<- FALSE
  if(na.rm == TRUE){
    ind.na <- which(is.na(y))
    if(length(ind.na) != 0){
      had.na<- TRUE
      y<- y[-ind.na]
    } 
  }
  
  sorted.y<- sort(abs(y),TRUE)
  ans<- (cumsum(sorted.y)-lam)/(1:length(y))
  
  if(all(ans<=0)){
    ret.obj<- rep(0,length(y)) 
    
    if(had.na == TRUE){
      ret.obj<- append(ret.obj,NaN, after = ind.na-1 )
    }
    return(ret.obj)
  }
  statement<- ans > sorted.y
  if(all(statement==FALSE)){
    c<- tail(ans,1)
    ret.obj<- sign(y)*(pmin(c,abs(y)))
    
    if(had.na == TRUE){
      ret.obj<- append(ret.obj,NaN, after = ind.na-1 )
    }
    return(ret.obj)
    
  }else{
    ind<- which(statement==TRUE)[1]   
    c<- ans[ind-1]
    ret.obj<- sign(y)*(pmin(c,abs(y)))
    
    if(had.na == TRUE){
      ret.obj<- append(ret.obj,NaN, after = ind.na-1 )
    }
    return(ret.obj) 
  }
}

#####################################################################################
x<-rnorm(1000)
x<-rep(0,10)
lambda=5
y1=l_inf_prox(x,lambda)
y3=inf_vec(x,lambda)
max(abs(y1-y3))

