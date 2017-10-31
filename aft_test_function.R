#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------
W_omni=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
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
  for(i in 1:n){
    N_i_t[[i]]=(e_i_beta>=e_i_beta[i])*Delta[i]
  }
  #N_i_t
  
  Y_i_t=list(NA)
  for(i in 1:n){
    Y_i_t[[i]]=(e_i_beta<=e_i_beta[i])*1
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
  
  Mhat_i_t=mapply("-", N_i_t,lapply(lapply(
    Y_i_t,'*',dLambdahat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  obs_stat_omni=Reduce('+',mapply(function(x,y){x%*%t(y)},Mhat_i_t,pi_i_z,SIMPLIFY=FALSE))/sqrt(n)
  #obs_stat_omni
  
  result=list(e_i_beta,obs_stat_omni)
  
  names(result)=c("Resid","obs_stat_omni")

  return(result)
}
#W_omni()

W_form=function(b,Time,Delta,Covari,form=1){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
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
  for(i in 1:n){
    N_i_t[[i]]=(e_i_beta>=e_i_beta[i])*Delta[i]
  }
  #N_i_t
  
  Y_i_t=list(NA)
  for(i in 1:n){
    Y_i_t[[i]]=(e_i_beta<=e_i_beta[i])*1
  }
  #Y_i_t
  
  N_d_t=Reduce('+',N_i_t)
  #N_d_t
  
  S_0_t=Reduce('+',Y_i_t)
  #S_0_t
  
  S_1_t=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,as.list(data.frame(t(Covari))),SIMPLIFY=FALSE))
  #S_1_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t
  
  Mhat_i_t=mapply("-", N_i_t,lapply(lapply(
    Y_i_t,'*',dLambdahat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
  #Mhat_i_inf
  
  obs_stat_form=Reduce('+',mapply('*',pi_i_z,Mhat_i_inf,SIMPLIFY=FALSE))/sqrt(n)
  #obs_stat_form
  
  result=list(e_i_beta,obs_stat_form)
  
  names(result)=c("Resid","obs_stat_form")
  
  return(result)
}
#W_form()

W_link=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
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
  for(i in 1:n){
    N_i_t[[i]]=(e_i_beta>=e_i_beta[i])*Delta[i]
  }
  #N_i_t
  
  Y_i_t=list(NA)
  for(i in 1:n){
    Y_i_t[[i]]=(e_i_beta<=e_i_beta[i])*1
  }
  #Y_i_t
  
  N_d_t=Reduce('+',N_i_t)
  #N_d_t
  
  S_0_t=Reduce('+',Y_i_t)
  #S_0_t
  
  S_1_t=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,as.list(data.frame(t(Covari))),SIMPLIFY=FALSE))
  #S_1_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t
  
  Mhat_i_t=mapply("-", N_i_t,lapply(lapply(
          Y_i_t,'*',dLambdahat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
  #Mhat_i_inf
  
  obs_stat_link=(1/sqrt(n))*Reduce('+',mapply(function(x,y){x*y},
          pi_i_z,Mhat_i_inf,SIMPLIFY=FALSE))
  #obs_stat_link

  result=list(e_i_beta,obs_stat_link)
  
  names(result)=c("Resid","obs_stat_link")
  
  return(result)
}
#W_link()

#-------------------------------------------------------------
#-------------------------REALIZATION-------------------------
#-------------------------------------------------------------
What_omni=function(b,std,Time,Delta,Covari,tol){
  # b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  # b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
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
  
  tolerance=tol+1 #initial value
  
  # while (tolerance>tol){
  
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
  
  sim_stat_omni=F.T.-S.T.-T.T.
  #sim_stat_omni
  
  return(sim_stat_omni)
}
#What_omni()

What_form=function(b,std,Time,Delta,Covari,tol,form=1){
  # b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;form=1
  # b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;form=1
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
  
  fhat_inf.z=lapply(fhat_t.z,function(x){x[,n]})
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
  
  tolerance=tol+1 #initial value
  
  # while (tolerance>tol){
  
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
  
  sim_stat_form=F.T.-S.T.-T.T.
  #sim_stat_form
  
  return(sim_stat_form)
}
#What_form()

What_link=function(b,std,Time,Delta,Covari,tol){
  # b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  # b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
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
  
  fhat_inf.z=lapply(fhat_t.z,function(x){x[,n]})
  
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
  
  tolerance=tol+1 #initial value
  
  # while (tolerance>tol){
  
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
  
  sim_stat_link=F.T.-S.T.-T.T.
  #sim_stat_link
  
  return(sim_stat_link)
}
#What_link()

#-------------------------------------------------------------
#--------------------------SIMULATION-------------------------
#-------------------------------------------------------------
sample_path_omni=function(path,b,std,Time,Delta,Covari,tol){
  #path=path;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #path=path;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  #------------------------SAMPLE PATH------------------------
  
  path_check=ceiling(path/2)
  
  dataset_What=list(NA)
  for(k in 1:path){
    dataset_What[[k]]=What_omni(b,std,Time,Delta,Covari,tol)
    if(k%%path_check==0) {
      cat("Sample Path",k,"\n")
    }
  }
  
  #------------------------BOOTSTRAPPING----------------------
  std.boot=matrix(apply(mapply(function(x){as.vector(x)},dataset_What),1,sd),nrow=n)
  # std.boot
  
  dataset_std.What=lapply(dataset_What,function(x){x/std.boot})
  # dataset_std.What
  
  dataset_W_omni=W_omni(b,Time,Delta,Covari)
  # dataset_W_omni
  
  dataset_W=dataset_W_omni$obs_stat_omni
  # dataset_W
  
  dataset_std.W=dataset_W/std.boot
  # dataset_std.W
  
  #-----------------------MAXIMUM VALUE-----------------------
  max_path_What=unlist(lapply(dataset_What,function(x){max(abs(x))}))
  # max_path_What
  
  max_path_W=max(abs(dataset_W))
  # max_path_W

  max_path_std.What=unlist(lapply(dataset_std.What,function(x){max(abs(x))}))
  # max_path_std.What
  
  max_path_std.W=max(abs(dataset_std.W))
  # max_path_std.W

  #--------------------------P VALUE--------------------------
  p_value=length(which((max_path_What>max_path_W)*1==1))/path
  # p_value
  matrix(c(1,3,2,6,5,6,7,8,9),nrow=3)
  std.p_value=length(which((max_path_std.What>max_path_std.W)*1==1))/path
  # std.p_value
  
  result=list(Time,Delta,Covari,dataset_W_omni$Resid,
              dataset_W,dataset_What,
              dataset_std.W,dataset_std.What,
              std.boot,p_value,std.p_value)

  names(result)=c("Time","Delta","Covari","Resid",
                  "dataset_W","dataset_What",
                  "dataset_std.W","dataset_std.What",
                  "std.boot","p_value","std.p_value")
  # result
  
  return(result)
}
#sample_path_omni

sample_path_form=function(path,b,std,Time,Delta,Covari,tol,form=1){
  #path=path;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #path=path;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  path_check=ceiling(path/2)
  
  #------------------------SAMPLE PATH------------------------
  dataset_What=list(NA)
  for(k in 1:path){
    dataset_What[[k]]=What_form(b,std,Time,Delta,Covari,tol,form=1)
    if(k%%path_check==0) {
      cat("Sample Path",k,"\n")
    }
  }
  
  #------------------------BOOTSTRAPPING----------------------
  std.boot=matrix(apply(mapply(function(x){as.vector(x)},dataset_What),1,sd),nrow=n)
  # std.boot
  
  dataset_std.What=lapply(dataset_What,function(x){x/std.boot})
  # dataset_std.What
  
  dataset_W_form=W_form(b,Time,Delta,Covari,form=1)
  # dataset_W_form
  
  dataset_W=dataset_W_form$obs_stat_form
  # dataset_W
  
  dataset_std.W=dataset_W/std.boot
  # dataset_std.W
  
  #-----------------------MAXIMUM VALUE-----------------------
  max_path_What=unlist(lapply(dataset_What,function(x){max(abs(x))}))
  # max_path_What
  
  max_path_W=max(abs(dataset_W))
  # max_path_W
  
  max_path_std.What=unlist(lapply(dataset_std.What,function(x){max(abs(x))}))
  # max_path_std.What
  
  max_path_std.W=max(abs(dataset_std.W))
  # max_path_std.W
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_path_What>max_path_W)*1==1))/path
  # p_value
  
  std.p_value=length(which((max_path_std.What>max_path_std.W)*1==1))/path
  # std.p_value
  
  result=list(Time,Delta,Covari,dataset_W_form$Resid,
              dataset_W,dataset_What,
              dataset_std.W,dataset_std.What,
              std.boot,p_value,std.p_value)
  
  names(result)=c("Time","Delta","Covari","Resid",
                  "dataset_W","dataset_What",
                  "dataset_std.W","dataset_std.What",
                  "std.boot","p_value","std.p_value")
  # result
  
  return(result)
}
#sample_path_form

sample_path_link=function(path,b,std,Time,Delta,Covari,tol){
  #path=path;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #path=path;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  path_check=ceiling(path/2)
  
  #------------------------SAMPLE PATH------------------------
  dataset_What=list(NA)
  for(k in 1:path){
    dataset_What[[k]]=What_link(b,std,Time,Delta,Covari,tol)
    if(k%%path_check==0) {
      cat("Sample Path",k,"\n")
    }
  }
  
  #------------------------BOOTSTRAPPING----------------------
  std.boot=matrix(apply(mapply(function(x){as.vector(x)},dataset_What),1,sd),nrow=n)
  # std.boot
  
  dataset_std.What=lapply(dataset_What,function(x){x/std.boot})
  # dataset_std.What
  
  dataset_W_link=W_link(b,Time,Delta,Covari)
  # dataset_W_link
  
  dataset_W=dataset_W_link$obs_stat_link
  # dataset_W
  
  dataset_std.W=dataset_W/std.boot
  # dataset_std.W
  
  #-----------------------MAXIMUM VALUE-----------------------
  max_path_What=unlist(lapply(dataset_What,function(x){max(abs(x))}))
  # max_path_What
  
  max_path_W=max(abs(dataset_W))
  # max_path_W
  
  max_path_std.What=unlist(lapply(dataset_std.What,function(x){max(abs(x))}))
  # max_path_std.What
  
  max_path_std.W=max(abs(dataset_std.W))
  # max_path_std.W
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_path_What>max_path_W)*1==1))/path
  # p_value
  
  std.p_value=length(which((max_path_std.What>max_path_std.W)*1==1))/path
  # std.p_value
  
  result=list(Time,Delta,Covari,dataset_W_link$Resid,
               dataset_W,dataset_What,
               dataset_std.W,dataset_std.What,
               std.boot,p_value,std.p_value)
  
  names(result)=c("Time","Delta","Covari","Resid",
                  "dataset_W","dataset_What",
                  "dataset_std.W","dataset_std.What",
                  "std.boot","p_value","std.p_value")
  # result
  
  return(result)
}
#sample_path_link

#-------------------------------------------------------------
#---------------------------AFTTEST---------------------------
#-------------------------------------------------------------
afttest=function(formula,dataset,testtype="omni",path=200,tol=0.1,ftn.form){
  
  varnames=noquote(all.vars(formula))
  len.var=length(varnames)
  
  Time=dataset[,which(names(dataset)==varnames[1])]
  Delta=dataset[,which(names(dataset)==varnames[2])]
  
  n=length(Time)
  
  Covari=matrix(dataset[,which(names(dataset)==varnames[3])],nrow=n)
  
  if(len.var>=4){
    for(i in 4:len.var){
      Covari=cbind(Covari,matrix(dataset[,which(names(dataset)==varnames[i])],nrow=n))
    }
  }
  
  if(testtype=="form"){
    if(length(which(varnames==ftn.form[1]))==0){return(print("Check your code"))}
    if(is.na(ftn.form)){form=1}
    else{form=which(varnames==ftn.form[1])-2}
  }
  
  aftsrr_result=aftgee::aftsrr(formula,method="nonsm")
  b=-as.vector(aftsrr_result$beta)
  std=diag(aftsrr_result$covmat$ISMB)
  
  if(testtype=="omni"){
    return(c(list(testtype=testtype),sample_path_omni(path,b,std,Time,Delta,Covari,tol)))
  }
  if(testtype=="form"){
    return(c(list(testtype=testtype),sample_path_form(path,b,std,Time,Delta,Covari,tol,form)))
  }
  if(testtype=="link"){
    return(c(list(testtype=testtype),sample_path_link(path,b,std,Time,Delta,Covari,tol)))
  }
  # if(testtype=="aft"){
  #
  # }
  return(print("Check your code"))
}

#-------------------------------------------------------------
#---------------------------PLOTTING--------------------------
#-------------------------------------------------------------
plotting_omni=function(result,xaxix,path){
  
  result_Time=result$Time
  n=length(result_Time)
  if (xaxix=="rank"){xaxix=(1:n)[order(result_Time)]}
  else {xaxix=result_Time}
  
  med=ceiling(sqrt(length(result$std.boot))/2)
  
  dataset_What=data.frame()
  
  for (i in 1:path){
    group=i
    A=result$app_path[[i]][,med]
    AA=data.frame(group,t_i=xaxix,What=A)
    dataset_What=rbind(dataset_What,AA)
  }
  #dataset_What
  
  dataset_W=data.frame(group,t_i=xaxix,W=result$obs_path[,med])
  #dataset_W
  
  Figure1_W=
    ggplot()+
    geom_step(data=dataset_What,aes(x=t_i,y=What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
    geom_step(data=dataset_W,aes(x=t_i,y=W),colour="tomato",lwd=0.25)+
    theme_minimal()
  #Figure1_W
  
  dataset_std.What=data.frame()
  
  for (i in 1:path){
    group=i
    A=result$app_std.path[[i]][,med]
    AA=data.frame(group,t_i=xaxix,std.What=A)
    dataset_std.What=rbind(dataset_std.What,AA)
  }
  #dataset_std.What
  
  dataset_std.W=data.frame(group,t_i=xaxix,std.W=result$obs_std.path[,med])
  #dataset_std.W
  
  Figure1_std.W=
    ggplot()+
    geom_step(data=dataset_std.What,aes(x=t_i,y=std.What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
    geom_step(data=dataset_std.W,aes(x=t_i,y=std.W),colour="tomato",lwd=0.25)+
    theme_minimal()
  #Figure1_std.W
  
  return(grid.arrange(Figure1_W,Figure1_std.W,nrow=2))
}

plotting_form=function(result,xaxix,path){
  
  result_Covari=result$Covari
  n=length(result_Covari)
  if (xaxix=="rank"){xaxix=(1:n)[order(result_Covari)]}
  else {xaxix=result_Covari}
  
  dataset_What=data.frame()
  
  for (i in 1:path){
    group=i
    A=result$app_path[[i]]
    AA=data.frame(group,z_i=xaxix,What=A)
    dataset_What=rbind(dataset_What,AA)
  }
  #dataset_What
  
  dataset_W=data.frame(group,z_i=xaxix,W=result$obs_path)
  #dataset_W
  
  Figure1_W=
    ggplot()+
    geom_step(data=dataset_What,aes(x=z_i,y=What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
    geom_step(data=dataset_W,aes(x=z_i,y=W),colour="tomato",lwd=0.25)+
    theme_minimal()
  #Figure1_W
  
  dataset_std.What=data.frame()
  
  for (i in 1:path){
    group=i
    A=result$app_std.path[[i]]
    AA=data.frame(group,z_i=xaxix,std.What=A)
    dataset_std.What=rbind(dataset_std.What,AA)
  }
  #dataset_std.What
  
  dataset_std.W=data.frame(group,z_i=xaxix,std.W=result$obs_std.path)
  #dataset_std.W
  
  Figure1_std.W=
    ggplot()+
    geom_step(data=dataset_std.What,aes(x=z_i,y=std.What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
    geom_step(data=dataset_std.W,aes(x=z_i,y=std.W),colour="tomato",lwd=0.25)+
    theme_minimal()
  #Figure1_std.W
  
  return(grid.arrange(Figure1_W,Figure1_std.W,nrow=2))
}

plotting_link=function(result,xaxix,path){
  
  result_Covari=result$Covari
  n=length(result_Covari)
  if (xaxix=="rank"){xaxix=(1:n)[order(result_Covari)]}
  else {xaxix=result_Covari}
  
  dataset_What=data.frame()
  
  for (i in 1:path){
    group=i
    A=result$app_path[[i]]
    AA=data.frame(group,z_i=xaxix,What=A)
    dataset_What=rbind(dataset_What,AA)
  }
  #dataset_What
  
  dataset_W=data.frame(group,z_i=xaxix,W=result$obs_path)
  #dataset_W
  
  Figure1_W=
    ggplot()+
    geom_step(data=dataset_What,aes(x=z_i,y=What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
    geom_step(data=dataset_W,aes(x=z_i,y=W),colour="tomato",lwd=0.25)+
    theme_minimal()
  #Figure1_W
  
  dataset_std.What=data.frame()
  
  for (i in 1:path){
    group=i
    A=result$app_std.path[[i]]
    AA=data.frame(group,z_i=xaxix,std.What=A)
    dataset_std.What=rbind(dataset_std.What,AA)
  }
  #dataset_std.What
  
  dataset_std.W=data.frame(group,z_i=xaxix,std.W=result$obs_std.path)
  #dataset_std.W
  
  Figure1_std.W=
    ggplot()+
    geom_step(data=dataset_std.What,aes(x=z_i,y=std.What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
    geom_step(data=dataset_std.W,aes(x=z_i,y=std.W),colour="tomato",lwd=0.25)+
    theme_minimal()
  #Figure1_std.W
  
  return(grid.arrange(Figure1_W,Figure1_std.W,nrow=2))
}

afttestplot=function(result,xaxix,path=50){
  
  testtype=result$testtype
  
  if(testtype=="omni"){
    return(plotting_omni(result,xaxix,path))
  }
  if(testtype=="form"){
    return(plotting_form(result,xaxix,path))
  }
  if(testtype=="link"){
    return(plotting_link(result,xaxix,path))
  }
  # if(testtype=="aft"){
  #   return(plotting_aft(result,path))
  # }
}

#-------------------------------------------------------------
#---------------------------EXERCISE--------------------------
#-------------------------------------------------------------
# dataset_cox=data.frame(X_ln_cox,D_ln_cox,Z_ln_cox)
# asdf1=afttest(Surv(X_ln_cox,D_ln_cox)~Z_ln_cox,dataset_cox,"form",50,0.1,"Z_ln_cox")
# afttestplot(asdf1,"real",30)
