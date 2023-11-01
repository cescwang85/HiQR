## qr1 function in R
qr1_R<-function(X,Y,lambda,err_abs=10^(-4),err_rel=10^(-3),maxIter=200,rho=1){
  n=nrow(X);
  p=ncol(X);
  lambda=sort(lambda, decreasing =TRUE);
  nlambda=length(lambda);
  ##Centering int##
  D0=gram(X,Y/n);
  Dk=D0;
  H0=X%*%t(X);
  H2=H0*H0/n;
  H=solve(H2+rho*diag(n));
    Omega_all=NULL;
    niter=lambda;
    rholist=lambda;
  ##Initialization##
    A=matrix(0,p,p);
    U=matrix(0,p,p);
    B=A;
    old_B=B;
    w=Y;
    ee_pri=1;
    ee_dual=1;
      for (k in 1:nlambda) {
        lam=lambda[k];
        i=0;
        while ((i<maxIter)||(i==0))
        {
          Dk=D0+rho*(B-U);
          ## A-update##
            w=H%*%qrow(X,Dk/n);
            A=(Dk-gram(X,w))/rho;
          ##B-update##
              old_B=B;
              B=soft(A+U,lam/rho);
          ##U-update##
              U=A-B+U;
              i=i+1;
          ##Stop Rule##
              ee_dual=rho*norm(B-old_B,type='F'); #dual residual
              ee_pri=norm(A-B,type='F');     #primal residual
              err_pri_new=p*err_abs+err_rel*max(norm(A,type='F'),norm(B,type='F'));
              err_dual_new=p*err_abs+err_rel*norm(U,type='F')*rho;
                      if ((ee_dual<err_dual_new)&&(ee_pri<err_pri_new)) break;
                      # /*Varying Penalty Parameter*/
                      #   /*if (ee_pri>10*ee_dual) {rho=2*rho;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
                      #   /*if (ee_dual>10*ee_pri) {rho=rho/2;H=inv_sympd(H2+rho*arma::eye(n,n));}*/
        }
        Omega_all=c(Omega_all,list(B));
        niter[k]=i;
        rholist[k]=rho;
      }
      return(list(Omega=Omega_all,lambda=lambda,rho=rholist,niter =niter))} 

## qr1 for single lambda function in R
qr1_single<-function(X,Y,lambda,err_abs=10^(-4),err_rel=10^(-3),maxIter=200,rho=1,rho_vary=1){
  n=nrow(X);
  p=ncol(X);
  ##Centering int##
  D0=gram(X,Y/n);
  Dk=D0;
  H0=X%*%t(X);
  H2=H0*H0/n;
  H=solve(H2+rho*diag(n));
  Omega_all=NULL;
  rholist=rep(0,maxIter)
  err_dual_list=rep(0,maxIter)
  err_pri_list=rep(0,maxIter)
  err_abs_list=rep(0,maxIter)
  ##Initialization##
  A=matrix(0,p,p);
  U=matrix(0,p,p);
  B=A;
  old_B=B;
  w=Y;
  ee_pri=1;
  ee_dual=1;
  
    lam=lambda;
    i=0;
    while ((i<maxIter)||(i==0))
    {
      Dk=D0+rho*(B-U);
      ## A-update##
      w=H%*%qrow(X,Dk/n);
      A=(Dk-gram(X,w))/rho;
      ##B-update##
      old_B=B;
      B=soft(A+U,lam/rho);
      ##U-update##
      U=A-B+U;
      i=i+1;
      ##Stop Rule##
      ee_dual=rho*norm(B-old_B,type='F'); #dual residual
      ee_pri=norm(A-B,type='F');     #primal residual
      err_pri_new=p*err_abs+err_rel*max(norm(A,type='F'),norm(B,type='F'));
      err_dual_new=p*err_abs+err_rel*norm(U,type='F')*rho;
      if ((ee_dual<err_dual_new)&&(ee_pri<err_pri_new)) break;
      # /*Varying Penalty Parameter*/
      if (rho_vary>0){
     if (ee_pri>10*ee_dual) {rho=2*rho;H=solve(H2+rho*diag(n));}
    if (ee_dual>10*ee_pri) {rho=rho/2;H=solve(H2+rho*diag(n));}
      }
      Omega_all=c(Omega_all,list(B));
      err_dual_list[i]=ee_dual;
      err_pri_list[i]=ee_pri;
      err_abs_list[i]=norm(B-old_B,type='M')
      rholist[i]=rho;
    }
  return(list(Omega=Omega_all,lambda=lambda,rho=rholist,dual_err=err_dual_list,pri_err=err_pri_list,abs_err=err_abs_list))} 

qr3_rank_single<-function(X,Y,lambda1,lambda2,err_abs=10^(-4),err_rel=10^(-3),maxIter=200,rho=5,rho_vary=1){
  n=nrow(X);
  p=ncol(X);
D0=gram(X,Y/n);
  Dk=D0;
H0=X%*%t(X);
H2=H0*H0/n;
H=solve(H2+rho*diag(n));
Omega_all=NULL;
rholist=rep(0,maxIter)
    w=Y;
    B=matrix(0,p,p);
    U1=B;
    U2=B;
    U3=B;

  lam1=lambda1;
  lam2=lambda2;
  ee_pri=1;
  ee_dual=1;
    i=0;
      while ((i<maxIter)||(i==0))
      {
          Dk=D0+rho*(B-U1);
          w=H%*%qrow(X,Dk/n);
          B1=Dk/rho-gram(X,w/rho);
          B2=soft(B-U2,lam1/rho);
          B3=nuclear_mat(B-U3,lam2/rho);
        
            old_B=B;
            B=(B1+B2+B3)/3;
            
              U1=U1+B1-B;
              U2=U2+B2-B;
              U3=U3+B3-B;
              i=i+1;
         
                ee_dual=rho*norm(B-old_B,type ='F'); #dual residual*/
                  ee_pri=norm(B1-B,type ='F')+norm(B2-B,type ='F')+norm(B3-B,type ='F'); #primal residual*/
                  err_pri_new=p*err_abs+err_rel*norm(B,type ='F');
                  err_dual_new=p*err_abs+err_rel*norm(U1,type ='F')*rho;
                    if ((ee_dual<err_dual_new)&&(ee_pri<err_pri_new)) break;
                    #/*Varying Penalty Parameter*/
                  if (rho_vary>0){
                      if (ee_pri>10*ee_dual) {rho=2*rho;H=solve(H2+rho*diag(n));}
                      if (ee_dual>10*ee_pri) {rho=rho/2;H=solve(H2+rho*diag(n));}
                  }  
                  Omega_all=c(Omega_all,list(B));
                  rholist[i]=rho;
      }
    return(list(Omega=Omega_all,rho=rholist))} 



             
             
             