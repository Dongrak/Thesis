afttest_omni=function(path,b,std,Time,Delta,Covari,tol){
  # path=200;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  # path=200;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  # path=200;b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parameters
  
  Covari=matrix(Covari,nrow=n)
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  pi_i_z=list(NA)
  for(i in 1:n){
    pi_i_z[[i]]=apply(apply(Covari,2,function(x){(x<=((x[order(x)])[i]))*1}),1,prod)
  }
  pi_i_z=as.list(data.frame(t(matrix(unlist(pi_i_z),nrow=n))))
  
  N_i_t=list(NA)
  for(j in 1:n){
    N_i_t[[j]]=(e_i_beta>=e_i_beta[j])*Delta[j]
  }
  #N_i_t
  
  Y_i_t=list(NA)
  for(j in 1:n){
    Y_i_t[[j]]=(e_i_beta<=e_i_beta[j])*1
  }
  #Y_i_t
  
  N_d_t=Reduce('+',N_i_t)
  #N_d_t
  
  S_0_t=Reduce('+',Y_i_t)
  #S_0_t
  
  S_1_t=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,
                          as.list(data.frame(t(Covari))),SIMPLIFY=FALSE))
  #S_1_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t 
  
  Mhat_i_t=mapply("-",N_i_t,lapply(lapply(
    Y_i_t,'*',dLambdahat_0_t),cumsum),SIMPLIFY=FALSE)
  #Mhat_i_t
  
  obs_path=Reduce('+',mapply(function(x,y){x%*%t(y)},
                             Mhat_i_t,pi_i_z,SIMPLIFY=FALSE))/sqrt(n)
  #obs_path
  
  S_pi_t.z=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,pi_i_z,SIMPLIFY=FALSE))
  #S_pi_t.z
  
  dMhat_i_t=lapply(Mhat_i_t,function(x){diff(c(0,x))})
  #dMhat_i_t
  
  #-----------------------------------------------------------
  #----------------------kernel Smoothing---------------------
  #-----------------------------------------------------------
  
  #-----------------------------g0----------------------------
  Ghat_0_t=1-exp(-Lambdahat_0_t)
  #Ghat_0_t
  
  dGhat_0_t=diff(c(0,Ghat_0_t))
  #dGhat_0_t
  
  ghat_0_t=(ksmooth(e_i_beta,dGhat_0_t,"normal",
                    bandwidth = 1.06*sd(dGhat_0_t)*n^(-0.2),x.points=e_i_beta)$y)
  #ghat_0_t
  
  ghat_t.z=list(NA)
  for(j in 1:p){
    ghat_t.z[[j]]=Reduce('+',lapply(mapply('*',pi_i_z,Covari[,j],SIMPLIFY=FALSE),
                                    function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n
  }
  #ghat_t.z
  
  #-----------------------------f0----------------------------
  Fhat_0_e=1-cumprod(1-(Delta/S_0_t))
  #Fhat_0_e
  
  dFhat_0_e=diff(c(0,Fhat_0_e))
  #dFhat_0_e
  
  Condi.Ehat=(cumsum(e_i_beta*dFhat_0_e))/(1-Fhat_0_e)
  #Condi.Ehat
  
  rhat_i=Delta*e_i_beta+(1-Delta)*Condi.Ehat
  rhat_i[is.nan(rhat_i)]=0
  #rhat_i
  
  den.f=ksmooth(rhat_i,dFhat_0_e,"normal",
                bandwidth = 1.06*sd(dFhat_0_e)*n^(-0.2),x.points=rhat_i)
  #den.f
  
  fhat_0_t=predict(loess(den.f$y~den.f$x),e_i_beta)
  #fhat_0_t
  
  fhat_t.z=list(NA)
  for(j in 1:p){
    fhat_t.z[[j]]=Reduce('+',lapply(mapply('*',pi_i_z,Delta*Covari[,j],SIMPLIFY=FALSE),
                                    function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n
  }
  #fhat_t.z
  
  #-----------------------------------------------------------
  #--------Find Beta_hat_star by using optimize function------
  #-----------------------------------------------------------
  U_beta=function(beta_U){
    #beta_U=b;
    
    Time_U=Time;Delta_U=Delta;Covari_U=Covari;
    
    e_i_beta_U=as.vector(log(Time_U)+Covari_U%*%beta_U)
    
    order_resid_U=order(e_i_beta_U)
    
    Time_U=Time_U[order_resid_U]
    Covari_U=matrix(Covari_U[order_resid_U,],nrow=n)
    Delta_U=Delta_U[order_resid_U]
    e_i_beta_U=e_i_beta_U[order_resid_U]
    
    N_i_t_U=list(NA)
    for(j in 1:n){
      N_i_t_U[[j]]=(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]
    }
    #N_i_t_U
    
    dN_i_t_U=lapply(N_i_t_U,function(x){diff(c(0,x))})
    #dN_i_t_U
    
    Y_i_t_U=list(NA)
    for(j in 1:n){
      Y_i_t_U[[j]]=(e_i_beta_U<=e_i_beta_U[j])*1
    }
    #Y_i_t_U
    
    S_0_t_U=Reduce('+',Y_i_t_U)
    #S_0_t_U
    
    S_1_t_U=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t_U,as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))
    #S_1_t_U
    
    U_inf_U=apply(S_0_t_U*Reduce('+',mapply(function(x,y){x%*%t(y)},dN_i_t_U,
                                            as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))-S_1_t_U*Reduce('+',dN_i_t_U),2,sum)/n
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  # tolerance=tol+1 #initial value
  
  # while (tolerance>tol){
  
  #------------------------SAMPLE PATH-----------------------
  
  app_path=list(NA)
  
  path_check=ceiling(path/2)
  
  for (k in 1:path){
    
    if(k%%path_check==0) {
      cat("Sample Path",k,"\n")
    }
    
    phi_i=rnorm(n)
    #phi_i
    
    U_pi_phi_t.z=apply(S_0_t*Reduce('+',mapply('*',mapply(function(x,y){x%*%t(y)},dMhat_i_t,
                                                          pi_i_z,SIMPLIFY=FALSE),phi_i,SIMPLIFY=FALSE))-S_pi_t.z*Reduce('+',mapply('*',
                                                                                                                                   dMhat_i_t,phi_i,SIMPLIFY=FALSE)),2,cumsum)/n
    #U_pi_phi_t.z
    
    U_phi_inf=apply(S_0_t*Reduce('+',mapply('*',mapply(function(x,y){x%*%t(y)},dMhat_i_t,
                                                       as.list(data.frame(t(Covari))),SIMPLIFY=FALSE),phi_i,SIMPLIFY=FALSE))-S_1_t*
                      Reduce('+',mapply('*',dMhat_i_t,phi_i,SIMPLIFY=FALSE)),2,sum)/n
    #U_phi_inf
    
    if(p==1){
      beta_hat_s_list=optimize(function(BETA){sum((U_beta(BETA)-U_phi_inf)^2)},
                               c(b-2*std,b+2*std),tol = 1e-16)
      #beta_hat_s_list
      
      beta_hat_s=beta_hat_s_list$minimum
      #beta_hat_s
      
      tolerance=beta_hat_s_list$objective
      #tolerance
    }
    if(p>1){
      beta_hat_s_list=optim(b,function(BETA){sum((U_beta(BETA)-U_phi_inf)^2)})
      #beta_hat_s_list
      
      beta_hat_s=beta_hat_s_list$par
      #beta_hat_s
      
      tolerance=beta_hat_s_list$value
      #tolerance
    }
    # }
    
    e_i_beta_s=as.vector(log(Time)+Covari%*%beta_hat_s)
    
    order_resid_s=order(e_i_beta_s)
    
    Delta_s=Delta[order_resid_s]
    e_i_beta_s=e_i_beta_s[order_resid_s]
    
    N_i_t_s=list(NA)
    for(j in 1:n){
      N_i_t_s[[j]]=(e_i_beta_s>=e_i_beta_s[j])*Delta_s[j]
    }
    #N_i_t_s
    
    Y_i_t_s=list(NA)
    for(j in 1:n){
      Y_i_t_s[[j]]=(e_i_beta_s<=e_i_beta_s[j])*1
    }
    #Y_i_t_s
    
    N_d_t_s=Reduce('+',N_i_t_s)
    #N_d_s_s_t_s
    
    S_0_t_s=Reduce('+',Y_i_t_s)
    #S_0_t_s
    
    J_t_s=(S_0_t_s>0)*1
    #J_t_s
    
    dN_d_t_s=diff(c(0,N_d_t_s))
    #dN_d_t_s
    
    Lambdahat_0_t_s=cumsum((J_t_s/S_0_t_s)*dN_d_t_s)
    #Lambdahat_0_t_s
    
    F.T.=(1/sqrt(n))*U_pi_phi_t.z
    S.T.=sqrt(n)*Reduce('+',mapply('*',mapply('+',fhat_t.z,mapply(function(x){apply(x,2,
      cumsum)},lapply(ghat_t.z,'*',dLambdahat_0_t),SIMPLIFY=FALSE),SIMPLIFY=FALSE),
      (b-beta_hat_s),SIMPLIFY=FALSE))
    T.T.=apply((S_pi_t.z*diff(c(0,Lambdahat_0_t-Lambdahat_0_t_s))),2,cumsum)/sqrt(n)
    
    app_path[[k]]=F.T.-S.T.-T.T.
    #app_path
  }
  
  std.boot=matrix(apply(mapply(function(x){as.vector(x)},app_path),1,sd),nrow=n)
  # std.boot
  
  app_std.path=lapply(app_path,function(x){x/std.boot})
  # app_std.path
  
  obs_std.path=obs_path/std.boot
  # obs_std.path
  
  #-----------------------MAXIMUM VALUE-----------------------
  max_app_path=unlist(lapply(app_path,function(x){max(abs(x))}))
  # max_app_path
  
  max_obs_path=max(abs(obs_path))
  # max_obs_path
  
  max_app_std.path=unlist(lapply(app_std.path,function(x){max(abs(x))}))
  # max_app_std.path
  
  max_obs_std.path=max(abs(obs_std.path))
  # max_obs_std.path
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_app_path>max_obs_path)*1==1))/path
  # p_value
  
  std.p_value=length(which((max_app_std.path>max_obs_std.path)*1==1))/path
  # std.p_value
  
  result=list(Time,Delta,Covari,e_i_beta,std.boot,
              app_path,app_std.path,
              obs_path,obs_std.path,
              p_value,std.p_value)
  
  names(result)=c("Time","Delta","Covari","Resid","std.boot",
                  "app_path","app_std.path",
                  "obs_path","obs_std.path",
                  "p_value","std.p_value")
  # result
  
  return(result)
}
#afttest_omni()

aa=afttest_omni(200,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
aa$p_value
aa$std.p_value
plotting_omni(aa,"rank",50)

bb=afttest_omni(200,beta_hat_ln_cox,std_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
bb$p_value
bb$std.p_value
plotting_omni(bb,"rank",50)

afttest_form=function(path,b,std,Time,Delta,Covari,tol,form=1){
  # path=200;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;form=1
  # path=200;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;form=1
  # path=200;b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parameters
  
  Covari=matrix(Covari,nrow=n)
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  pi_i_z=list(NA)
  Covari_form=Covari[,form]
  for(i in 1:n){
    pi_i_z[[i]]=(Covari_form<=((Covari_form[order(Covari_form)])[i]))*1
  }
  pi_i_z=as.list(data.frame(t(matrix(unlist(pi_i_z),nrow=n))))
  
  N_i_t=list(NA)
  for(j in 1:n){
    N_i_t[[j]]=(e_i_beta>=e_i_beta[j])*Delta[j]
  }
  #N_i_t
  
  Y_i_t=list(NA)
  for(j in 1:n){
    Y_i_t[[j]]=(e_i_beta<=e_i_beta[j])*1
  }
  #Y_i_t
  
  N_d_t=Reduce('+',N_i_t)
  #N_d_t
  
  S_0_t=Reduce('+',Y_i_t)
  #S_0_t
  
  S_1_t=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,
                          as.list(data.frame(t(Covari))),SIMPLIFY=FALSE))
  #S_1_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t 
  
  Mhat_i_t=mapply("-",N_i_t,lapply(lapply(
    Y_i_t,'*',dLambdahat_0_t),cumsum),SIMPLIFY=FALSE)
  #Mhat_i_t
  
  Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
  #Mhat_i_inf
  
  obs_path=Reduce('+',mapply('*',pi_i_z,Mhat_i_inf,SIMPLIFY=FALSE))/sqrt(n)
  #obs_path
  
  S_pi_t.z=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,pi_i_z,SIMPLIFY=FALSE))
  #S_pi_t.z
  
  dMhat_i_t=lapply(Mhat_i_t,function(x){diff(c(0,x))})
  #dMhat_i_t
  
  #-----------------------------------------------------------
  #----------------------kernel Smoothing---------------------
  #-----------------------------------------------------------
  
  #-----------------------------g0----------------------------
  Ghat_0_t=1-exp(-Lambdahat_0_t)
  #Ghat_0_t
  
  dGhat_0_t=diff(c(0,Ghat_0_t))
  #dGhat_0_t
  
  ghat_0_t=(ksmooth(e_i_beta,dGhat_0_t,"normal",
                    bandwidth = 1.06*sd(dGhat_0_t)*n^(-0.2),x.points=e_i_beta)$y)
  #ghat_0_t
  
  ghat_t.z=list(NA)
  for(j in 1:p){
    ghat_t.z[[j]]=Reduce('+',lapply(mapply('*',pi_i_z,Covari[,j],SIMPLIFY=FALSE),
                                    function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n
  }
  #ghat_t.z
  
  #-----------------------------f0----------------------------
  Fhat_0_e=1-cumprod(1-(Delta/S_0_t))
  #Fhat_0_e
  
  dFhat_0_e=diff(c(0,Fhat_0_e))
  #dFhat_0_e
  
  Condi.Ehat=(cumsum(e_i_beta*dFhat_0_e))/(1-Fhat_0_e)
  #Condi.Ehat
  
  rhat_i=Delta*e_i_beta+(1-Delta)*Condi.Ehat
  rhat_i[is.nan(rhat_i)]=0
  #rhat_i
  
  den.f=ksmooth(rhat_i,dFhat_0_e,"normal",
                bandwidth = 1.06*sd(dFhat_0_e)*n^(-0.2),x.points=rhat_i)
  #den.f
  
  fhat_0_t=predict(loess(den.f$y~den.f$x),e_i_beta)
  #fhat_0_t
  
  fhat_t.z=list(NA)
  for(j in 1:p){
    fhat_t.z[[j]]=Reduce('+',lapply(mapply('*',pi_i_z,Delta*Covari[,j],SIMPLIFY=FALSE),
                                    function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n
  }
  #fhat_t.z
  
  fhat_inf.z=lapply(fhat_t.z,function(x){x[n,]})
  #fhat_inf.z
  
  #-----------------------------------------------------------
  #--------Find Beta_hat_star by using optimize function------
  #-----------------------------------------------------------
  U_beta=function(beta_U){
    #beta_U=b;
    
    Time_U=Time;Delta_U=Delta;Covari_U=Covari;
    
    e_i_beta_U=as.vector(log(Time_U)+Covari_U%*%beta_U)
    
    order_resid_U=order(e_i_beta_U)
    
    Time_U=Time_U[order_resid_U]
    Covari_U=matrix(Covari_U[order_resid_U,],nrow=n)
    Delta_U=Delta_U[order_resid_U]
    e_i_beta_U=e_i_beta_U[order_resid_U]
    
    N_i_t_U=list(NA)
    for(j in 1:n){
      N_i_t_U[[j]]=(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]
    }
    #N_i_t_U
    
    dN_i_t_U=lapply(N_i_t_U,function(x){diff(c(0,x))})
    #dN_i_t_U
    
    Y_i_t_U=list(NA)
    for(j in 1:n){
      Y_i_t_U[[j]]=(e_i_beta_U<=e_i_beta_U[j])*1
    }
    #Y_i_t_U
    
    S_0_t_U=Reduce('+',Y_i_t_U)
    #S_0_t_U
    
    S_1_t_U=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t_U,as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))
    #S_1_t_U
    
    U_inf_U=apply(S_0_t_U*Reduce('+',mapply(function(x,y){x%*%t(y)},dN_i_t_U,
                                            as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))-S_1_t_U*Reduce('+',dN_i_t_U),2,sum)/n
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  # tolerance=tol+1 #initial value
  
  # while (tolerance>tol){
  
  #------------------------SAMPLE PATH-----------------------
  
  app_path=list(NA)
  
  path_check=ceiling(path/2)
  
  for (k in 1:path){
    
    if(k%%path_check==0) {
      cat("Sample Path",k,"\n")
    }
    
    phi_i=rnorm(n)
    #phi_i
    
    U_pi_phi_inf.z=apply(S_0_t*Reduce('+',mapply('*',mapply(function(x,y){x%*%t(y)},dMhat_i_t,
                                                            pi_i_z,SIMPLIFY=FALSE),phi_i,SIMPLIFY=FALSE))-S_pi_t.z*Reduce('+',mapply('*',
                                                                                                                                     dMhat_i_t,phi_i,SIMPLIFY=FALSE)),2,sum)/n
    #U_pi_phi_inf.z
    
    U_phi_inf=apply(S_0_t*Reduce('+',mapply('*',mapply(function(x,y){x%*%t(y)},dMhat_i_t,
                                                       as.list(data.frame(t(Covari))),SIMPLIFY=FALSE),phi_i,SIMPLIFY=FALSE))-S_1_t*
                      Reduce('+',mapply('*',dMhat_i_t,phi_i,SIMPLIFY=FALSE)),2,sum)/n
    #U_phi_inf
    
    if(p==1){
      beta_hat_s_list=optimize(function(BETA){sum((U_beta(BETA)-U_phi_inf)^2)},
                               c(b-5*std,b+5*std),tol = 1e-16)
      #beta_hat_s_list
      
      beta_hat_s=beta_hat_s_list$minimum
      #beta_hat_s
      
      tolerance=beta_hat_s_list$objective
      #tolerance
      
    }
    
    if(p>1){
      beta_hat_s_list=optim(b,function(BETA){sum((U_beta(BETA)-U_phi_inf)^2)})
      #beta_hat_s_list
      
      beta_hat_s=beta_hat_s_list$par
      #beta_hat_s
      
      tolerance=beta_hat_s_list$value
      #tolerance
    }
    # }
    
    e_i_beta_s=as.vector(log(Time)+Covari%*%beta_hat_s)
    
    order_resid_s=order(e_i_beta_s)
    
    Delta_s=Delta[order_resid_s]
    e_i_beta_s=e_i_beta_s[order_resid_s]
    
    N_i_t_s=list(NA)
    for(j in 1:n){
      N_i_t_s[[j]]=(e_i_beta_s>=e_i_beta_s[j])*Delta_s[j]
    }
    #N_i_t_s
    
    Y_i_t_s=list(NA)
    for(j in 1:n){
      Y_i_t_s[[j]]=(e_i_beta_s<=e_i_beta_s[j])*1
    }
    #Y_i_t_s
    
    N_d_t_s=Reduce('+',N_i_t_s)
    #N_d_s_s_t_s
    
    S_0_t_s=Reduce('+',Y_i_t_s)
    #S_0_t_s
    
    J_t_s=(S_0_t_s>0)*1
    #J_t_s
    
    dN_d_t_s=diff(c(0,N_d_t_s))
    #dN_d_t_s
    
    Lambdahat_0_t_s=cumsum((J_t_s/S_0_t_s)*dN_d_t_s)
    #Lambdahat_0_t_s
    
    F.T.=U_pi_phi_inf.z/sqrt(n)
    S.T.=sqrt(n)*Reduce('+',mapply('*',mapply('+',fhat_inf.z,mapply(function(x){apply(x,2,sum)},
                                                                    lapply(ghat_t.z,'*',dLambdahat_0_t),SIMPLIFY=FALSE),SIMPLIFY=FALSE),
                                   (b-beta_hat_s),SIMPLIFY=FALSE))
    T.T.=apply((S_pi_t.z*diff(c(0,Lambdahat_0_t-Lambdahat_0_t_s))),2,sum)/sqrt(n)
    
    app_path[[k]]=F.T.-S.T.-T.T.
    #app_path
  }
  
  std.boot=apply(mapply(function(x){as.vector(x)},app_path),1,sd)
  # std.boot
  
  app_std.path=lapply(app_path,function(x){x/std.boot})
  # app_std.path
  
  obs_std.path=obs_path/std.boot
  # obs_std.path
  
  #-----------------------MAXIMUM VALUE-----------------------
  max_app_path=unlist(lapply(app_path,function(x){max(abs(x))}))
  # max_app_path
  
  max_obs_path=max(abs(obs_path))
  # max_obs_path
  
  max_app_std.path=unlist(lapply(app_std.path,function(x){max(abs(x))}))
  # max_app_std.path
  
  max_obs_std.path=max(abs(obs_std.path))
  # max_obs_std.path
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_app_path>max_obs_path)*1==1))/path
  # p_value
  
  std.p_value=length(which((max_app_std.path>max_obs_std.path)*1==1))/path
  # std.p_value
  
  result=list(Time,Delta,Covari,e_i_beta,std.boot,
              app_path,app_std.path,
              obs_path,obs_std.path,
              p_value,std.p_value)
  
  names(result)=c("Time","Delta","Covari","Resid","std.boot",
                  "app_path","app_std.path",
                  "obs_path","obs_std.path",
                  "p_value","std.p_value")
  # result
  
  return(result)
}
#afttest_form()

cc=afttest_form(200,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
cc$p_value
cc$std.p_value
plotting_form(cc,"rank",50)

dd=afttest_form(200,beta_hat_ln_aft_f,std_hat_ln_aft_f,X_ln_aft_f,D_ln_aft_f,Z_ln_aft_f,given_tol)
dd$p_value
dd$std.p_value
plotting_form(dd,"rank",50)

afttest_link=function(path,b,std,Time,Delta,Covari,tol){
  # path=200;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  # path=200;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  # b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parameters
  
  Covari=matrix(Covari,nrow=n)
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  pi_i_z=list(NA)
  for(i in 1:n){
    pi_i_z[[i]]=apply(apply(Covari,2,function(x){(x<=((x[order(x)])[i]))*1}),1,prod)
  }
  pi_i_z=as.list(data.frame(t(matrix(unlist(pi_i_z),nrow=n))))
  
  N_i_t=list(NA)
  for(j in 1:n){
    N_i_t[[j]]=(e_i_beta>=e_i_beta[j])*Delta[j]
  }
  #N_i_t
  
  Y_i_t=list(NA)
  for(j in 1:n){
    Y_i_t[[j]]=(e_i_beta<=e_i_beta[j])*1
  }
  #Y_i_t
  
  N_d_t=Reduce('+',N_i_t)
  #N_d_t
  
  S_0_t=Reduce('+',Y_i_t)
  #S_0_t
  
  S_1_t=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,
                          as.list(data.frame(t(Covari))),SIMPLIFY=FALSE))
  #S_1_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t 
  
  Mhat_i_t=mapply("-",N_i_t,lapply(lapply(
    Y_i_t,'*',dLambdahat_0_t),cumsum),SIMPLIFY=FALSE)
  #Mhat_i_t
  
  Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
  #Mhat_i_inf
  
  obs_path=Reduce('+',mapply('*',pi_i_z,Mhat_i_inf,SIMPLIFY=FALSE))/sqrt(n)
  #obs_path
  
  S_pi_t.z=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,pi_i_z,SIMPLIFY=FALSE))
  #S_pi_t.z
  
  dMhat_i_t=lapply(Mhat_i_t,function(x){diff(c(0,x))})
  #dMhat_i_t
  
  #-----------------------------------------------------------
  #----------------------kernel Smoothing---------------------
  #-----------------------------------------------------------
  
  #-----------------------------g0----------------------------
  Ghat_0_t=1-exp(-Lambdahat_0_t)
  #Ghat_0_t
  
  dGhat_0_t=diff(c(0,Ghat_0_t))
  #dGhat_0_t
  
  ghat_0_t=(ksmooth(e_i_beta,dGhat_0_t,"normal",
                    bandwidth = 1.06*sd(dGhat_0_t)*n^(-0.2),x.points=e_i_beta)$y)
  #ghat_0_t
  
  ghat_t.z=list(NA)
  for(j in 1:p){
    ghat_t.z[[j]]=Reduce('+',lapply(mapply('*',pi_i_z,Covari[,j],SIMPLIFY=FALSE),
                                    function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n
  }
  #ghat_t.z
  
  #-----------------------------f0----------------------------
  Fhat_0_e=1-cumprod(1-(Delta/S_0_t))
  #Fhat_0_e
  
  dFhat_0_e=diff(c(0,Fhat_0_e))
  #dFhat_0_e
  
  Condi.Ehat=(cumsum(e_i_beta*dFhat_0_e))/(1-Fhat_0_e)
  #Condi.Ehat
  
  rhat_i=Delta*e_i_beta+(1-Delta)*Condi.Ehat
  rhat_i[is.nan(rhat_i)]=0
  #rhat_i
  
  den.f=ksmooth(rhat_i,dFhat_0_e,"normal",
                bandwidth = 1.06*sd(dFhat_0_e)*n^(-0.2),x.points=rhat_i)
  #den.f
  
  fhat_0_t=predict(loess(den.f$y~den.f$x),e_i_beta)
  #fhat_0_t
  
  fhat_t.z=list(NA)
  for(j in 1:p){
    fhat_t.z[[j]]=Reduce('+',lapply(mapply('*',pi_i_z,Delta*Covari[,j],SIMPLIFY=FALSE),
                                    function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n
  }
  #fhat_t.z
  
  fhat_inf.z=lapply(fhat_t.z,function(x){x[n,]})
  
  #-----------------------------------------------------------
  #--------Find Beta_hat_star by using optimize function------
  #-----------------------------------------------------------
  U_beta=function(beta_U){
    #beta_U=b;
    
    Time_U=Time;Delta_U=Delta;Covari_U=Covari;
    
    e_i_beta_U=as.vector(log(Time_U)+Covari_U%*%beta_U)
    
    order_resid_U=order(e_i_beta_U)
    
    Time_U=Time_U[order_resid_U]
    Covari_U=matrix(Covari_U[order_resid_U,],nrow=n)
    Delta_U=Delta_U[order_resid_U]
    e_i_beta_U=e_i_beta_U[order_resid_U]
    
    N_i_t_U=list(NA)
    for(j in 1:n){
      N_i_t_U[[j]]=(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]
    }
    #N_i_t_U
    
    dN_i_t_U=lapply(N_i_t_U,function(x){diff(c(0,x))})
    #dN_i_t_U
    
    Y_i_t_U=list(NA)
    for(j in 1:n){
      Y_i_t_U[[j]]=(e_i_beta_U<=e_i_beta_U[j])*1
    }
    #Y_i_t_U
    
    S_0_t_U=Reduce('+',Y_i_t_U)
    #S_0_t_U
    
    S_1_t_U=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t_U,as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))
    #S_1_t_U
    
    U_inf_U=apply(S_0_t_U*Reduce('+',mapply(function(x,y){x%*%t(y)},dN_i_t_U,
                                            as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))-S_1_t_U*Reduce('+',dN_i_t_U),2,sum)/n
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  # tolerance=tol+1 #initial value
  
  # while (tolerance>tol){
  
  #------------------------SAMPLE PATH-----------------------
  
  app_path=list(NA)
  
  path_check=ceiling(path/2)
  
  for (k in 1:path){
    
    if(k%%path_check==0) {
      cat("Sample Path",k,"\n")
    }
    
    phi_i=rnorm(n)
    #phi_i
    
    U_pi_phi_inf.z=apply(S_0_t*Reduce('+',mapply('*',mapply(function(x,y){x%*%t(y)},dMhat_i_t,
                                                            pi_i_z,SIMPLIFY=FALSE),phi_i,SIMPLIFY=FALSE))-S_pi_t.z*Reduce('+',mapply('*',
                                                                                                                                     dMhat_i_t,phi_i,SIMPLIFY=FALSE)),2,sum)/n
    #U_pi_phi_inf.z
    
    U_phi_inf=apply(S_0_t*Reduce('+',mapply('*',mapply(function(x,y){x%*%t(y)},dMhat_i_t,
                                                       as.list(data.frame(t(Covari))),SIMPLIFY=FALSE),phi_i,SIMPLIFY=FALSE))-S_1_t*
                      Reduce('+',mapply('*',dMhat_i_t,phi_i,SIMPLIFY=FALSE)),2,sum)/n
    #U_phi_inf
    
    if(p==1){
      beta_hat_s_list=optimize(function(BETA){sum((U_beta(BETA)-U_phi_inf)^2)},
                               c(b-2*std,b+2*std),tol = 1e-16)
      #beta_hat_s_list
      
      beta_hat_s=beta_hat_s_list$minimum
      #beta_hat_s
      
      tolerance=beta_hat_s_list$objective
      #tolerance
    }
    if(p>1){
      beta_hat_s_list=optim(b,function(BETA){sum((U_beta(BETA)-U_phi_inf)^2)})
      #beta_hat_s_list
      
      beta_hat_s=beta_hat_s_list$par
      #beta_hat_s
      
      tolerance=beta_hat_s_list$value
      #tolerance
    }
    # }
    
    e_i_beta_s=as.vector(log(Time)+Covari%*%beta_hat_s)
    
    order_resid_s=order(e_i_beta_s)
    
    Delta_s=Delta[order_resid_s]
    e_i_beta_s=e_i_beta_s[order_resid_s]
    
    N_i_t_s=list(NA)
    for(j in 1:n){
      N_i_t_s[[j]]=(e_i_beta_s>=e_i_beta_s[j])*Delta_s[j]
    }
    #N_i_t_s
    
    Y_i_t_s=list(NA)
    for(j in 1:n){
      Y_i_t_s[[j]]=(e_i_beta_s<=e_i_beta_s[j])*1
    }
    #Y_i_t_s
    
    N_d_t_s=Reduce('+',N_i_t_s)
    #N_d_s_s_t_s
    
    S_0_t_s=Reduce('+',Y_i_t_s)
    #S_0_t_s
    
    J_t_s=(S_0_t_s>0)*1
    #J_t_s
    
    dN_d_t_s=diff(c(0,N_d_t_s))
    #dN_d_t_s
    
    Lambdahat_0_t_s=cumsum((J_t_s/S_0_t_s)*dN_d_t_s)
    #Lambdahat_0_t_s
    
    F.T.=U_pi_phi_inf.z/sqrt(n)
    S.T.=sqrt(n)*Reduce('+',mapply('*',mapply('+',fhat_inf.z,mapply(function(x){apply(x,2,sum)},
      lapply(ghat_t.z,'*',dLambdahat_0_t),SIMPLIFY=FALSE),SIMPLIFY=FALSE),
      (b-beta_hat_s),SIMPLIFY=FALSE))
    T.T.=apply((S_pi_t.z*diff(c(0,Lambdahat_0_t-Lambdahat_0_t_s))),2,sum)/sqrt(n)
    
    app_path[[k]]=F.T.-S.T.-T.T.
    #app_path
  }
  
  std.boot=apply(mapply(function(x){as.vector(x)},app_path),1,sd)
  # std.boot
  
  app_std.path=lapply(app_path,function(x){x/std.boot})
  # app_std.path
  
  obs_std.path=obs_path/std.boot
  # obs_std.path
  
  #-----------------------MAXIMUM VALUE-----------------------
  max_app_path=unlist(lapply(app_path,function(x){max(abs(x))}))
  # max_app_path
  
  max_obs_path=max(abs(obs_path))
  # max_obs_path
  
  max_app_std.path=unlist(lapply(app_std.path,function(x){max(abs(x))}))
  # max_app_std.path
  
  max_obs_std.path=max(abs(obs_std.path))
  # max_obs_std.path
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_app_path>max_obs_path)*1==1))/path
  # p_value
  
  std.p_value=length(which((max_app_std.path>max_obs_std.path)*1==1))/path
  # std.p_value
  
  result=list(Time,Delta,Covari,e_i_beta,std.boot,
              app_path,app_std.path,
              obs_path,obs_std.path,
              p_value,std.p_value)
  
  names(result)=c("Time","Delta","Covari","Resid","std.boot",
                  "app_path","app_std.path",
                  "obs_path","obs_std.path",
                  "p_value","std.p_value")
  # result
  
  return(result)
}
#afttest_link()

ee=afttest_link(200,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
ee$p_value
ee$std.p_value
plotting_link(ee,"rank",50)

ff=afttest_link(200,beta_hat_ln_aft_f,std_hat_ln_aft_f,X_ln_aft_f,D_ln_aft_f,Z_ln_aft_f,given_tol)
ff$p_value
ff$std.p_value
plotting_link(ff,"rank",50)



