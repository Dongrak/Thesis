#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------
W_omni=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parametersa
  
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
  
  E_t=S_1_t/S_0_t
  #E_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t
  
  Mhat_i_t=mapply("-", N_i_t,lapply(lapply(
    Y_i_t,"*",dLambdahat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  obs_stat=Reduce('+',mapply(function(x,y){x%*%t(y)},pi_i_z,Mhat_i_t,SIMPLIFY=FALSE))/sqrt(n)
  #obs_stat
  
  result=list(Time,Delta,Covari,e_i_beta,obs_stat)
  names(result)=c("Time","Delta","Covari","Resid","obs_stat")
  
  return(result)
}
#W_omni()

W_fform=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-5*Z_wb);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parametersa
  
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
    pi_i_z[[i]]=apply(Covari,2,function(x){(x<=((x[order(x)])[i]))*1})
  }
  pi_ij_z=list(NA)
  for(j in 1:p){
    pi_ij_z[[j]]=as.list(data.frame(t(matrix(unlist(lapply(pi_i_z,function(x){x[,j]})),nrow=n))))
  }
  
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
  
  E_t=S_1_t/S_0_t
  #E_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t
  
  Mhat_i_t=mapply("-", N_i_t,lapply(lapply(
    Y_i_t,"*",dLambdahat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
  #Mhat_i_inf
  
  obs_stat=matrix(unlist(lapply(lapply(pi_ij_z,function(z){mapply(function(x,y){x*y},
          z,Mhat_i_inf,SIMPLIFY=FALSE)}),function(x){Reduce('+',x)/sqrt(n)})),nrow=n)
  #obs_stat
  
  result=list(Time,Delta,Covari,e_i_beta,obs_stat)
  names(result)=c("Time","Delta","Covari","Resid","obs_stat")
  
  return(result)
}
#W_fform()

W_linkf=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-5*Z_wb);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parametersa
  
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
  
  E_t=S_1_t/S_0_t
  #E_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t
  
  Mhat_i_t=mapply("-", N_i_t,lapply(lapply(
          Y_i_t,"*",dLambdahat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
  #Mhat_i_inf
  
  obs_stat=(1/sqrt(n))*Reduce('+',mapply(function(x,y){x*y},
          pi_i_z,Mhat_i_inf,SIMPLIFY=FALSE))
  #obs_stat
  
  result=list(Time,Delta,Covari,e_i_beta,obs_stat)
  names(result)=c("Time","Delta","Covari","Resid","obs_stat")
  
  return(result)
}
#W_linkf()

W_t.z=function(b,Time,Delta,Covari,test){
  if(test=="omni"){
    return(W_omni(b,Time,Delta,Covari))
  }
  if(test=="fform"){
    return(W_fform(b,Time,Delta,Covari))
  }
  if(test=="linkf"){
    return(W_linkf(b,Time,Delta,Covari))
  }
  if(test=="aft"){
    return(print("NOT YET..."))
    #return(W_aft(b,Time,Delta,Covari))
  }
}
#W_t.z()

#-------------------------------------------------------------
#-------------------------REALIZATION-------------------------
#-------------------------------------------------------------
What_omni=function(b,std,Time,Delta,Covari,tol){
  #b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parametersa
  
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
  
  E_t=S_1_t/S_0_t
  #E_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t 
  
  Mhat_i_t=mapply("-",N_i_t,lapply(lapply(
    Y_i_t,"*",dLambdahat_0_t),cumsum),SIMPLIFY=FALSE)
  #Mhat_i_t
  
  dN_i_t=lapply(N_i_t,function(x){diff(c(0,x))})
  #dN_i_t
  
  psi_t=S_0_t/n # Gehan's weight
  #psi_t

  U_t=Reduce('+',lapply(mapply("*",lapply(lapply(as.list(data.frame(t(Covari))),function(x)
    {t(x-t(E_t))}),"*",psi_t),dN_i_t,SIMPLIFY=FALSE),function(x){apply(x,2,cumsum)}))
  #U_t

  S_pi_t.z=Reduce('+',mapply(function(x,y){x%*%t(y)},pi_i_z,Y_i_t,SIMPLIFY=FALSE))
  #S_pi_t.z
  
  E_pi_t.z=t(t(S_pi_t.z)/S_0_t)
  #E_pi_t.z
  
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
    ghat_t.z[[j]]=Reduce('+',lapply(mapply("*",pi_i_z,Covari[,j],SIMPLIFY=FALSE),
                                    '%*%',t(ghat_0_t*Time)))/n
  }
  #ghat_t.z
  
  #-----------------------------f0----------------------------
  KM_e=cumprod(1-(Delta/S_0_t))
  #KM_e
  
  Fhat_0_e=1-KM_e
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
    fhat_t.z[[j]]=Reduce('+',lapply(mapply("*",pi_i_z,Delta*Covari[,j],SIMPLIFY=FALSE),
                                    '%*%',t(ghat_0_t*Time)))/n
  }
  #fhat_t.z
  
  #-----------------------------------------------------------
  # ghat_t.z=list(NA)
  # for(j in 1:p){
  #   ghat_t.z[[j]]=Reduce('+',lapply(lapply(pi_i_z,"*",Covari[,j]),'%*%',t(ghat_0_t*Time)))/n
  # }
  # 
  # fhat_t.z=list(NA)
  # for(j in 1:p){
  #   fhat_t.z[[j]]=Reduce('+',lapply(lapply(pi_i_z,"*",Delta*Covari[,j]),'%*%',t(ghat_0_t*Time)))/n
  # }
  #-----------------------------------------------------------
  
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
    
    E_t_U=S_1_t_U/S_0_t_U
    #E_t_U
    
    psi_t_U=S_0_t_U/n # Gehan's weight
    #psi_t_U
    
    U_t_U=Reduce('+',lapply(mapply("*",lapply(lapply(as.list(data.frame(t(Covari_U))),
          function(x){t(x-t(E_t_U))}),"*",psi_t_U),dN_i_t_U,SIMPLIFY=FALSE),
          function(x){apply(x,2,cumsum)}))
    #U_t_U
    
    U_inf_U=U_t_U[n,]
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  tolerance=tol+1 #initial value
  
  while (tolerance>tol){
    
    phi_i=rnorm(n)
    #phi_i
    
    U_pi_phi_t.z=Reduce('+',mapply("*",lapply(mapply(function(x,y){t(t(x)*y)},
                                                     lapply(pi_i_z,function(x,y){x-E_pi_t.z}),lapply(dMhat_i_t,"*",psi_t),
                                                     SIMPLIFY=FALSE),function(x){t(apply(x,1,cumsum))}),phi_i,SIMPLIFY=FALSE))
    #U_pi_phi_t.z
    
    U_phi_t=Reduce('+',mapply("*",lapply(mapply("*",lapply(as.list(data.frame(
      t(Covari))),function(x){t(x-t(E_t))}),lapply(dMhat_i_t,"*",psi_t),
      SIMPLIFY=FALSE),function(x){apply(x,2,cumsum)}),phi_i,SIMPLIFY=FALSE))
    #U_phi_t
    
    U_phi_inf=U_phi_t[n,]
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
  }
  
  e_i_beta_s=as.vector(log(Time)+Covari%*%beta_hat_s)
  
  order_resid_s=order(e_i_beta_s)
  
  Time_s=Time[order_resid_s]
  Covari_s=matrix(Covari[order_resid_s,],nrow=n)
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

  S_1_t_s=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t_s,
              as.list(data.frame(t(Covari_s))),SIMPLIFY=FALSE))
  #S_1_t
  
  E_t_s=S_1_t_s/S_0_t_s
  #E_t

  J_t_s=(S_0_t_s>0)*1
  #J_t_s
  
  dN_d_t_s=diff(c(0,N_d_t_s))
  #dN_d_t_s
  
  Lambdahat_0_t_s=cumsum((J_t_s/S_0_t_s)*dN_d_t_s)
  #Lambdahat_0_t_s
  
  dLambdahat_0_t_s=diff(c(0,Lambdahat_0_t_s))
  #dLambdahat_0_t_s

  F.T.=(1/sqrt(n))*U_pi_phi_t.z
  S.T.=sqrt(n)*Reduce("+",mapply("*",mapply("+",fhat_t.z,mapply(function(x){
    t(apply(x,2,cumsum))},lapply(ghat_t.z,function(x,y){t(x)*y},dLambdahat_0_t)
    ,SIMPLIFY=FALSE),SIMPLIFY=FALSE),(b-beta_hat_s),SIMPLIFY=FALSE))
  T.T.=(1/sqrt(n))*t(apply((t(S_pi_t.z)*diff(c(0,Lambdahat_0_t-Lambdahat_0_t_s))),2,cumsum))
  
  sim_stat=F.T.-S.T.-T.T.
  #sim_stat
  
  result=list(Time,Delta,Covari,e_i_beta,sim_stat)
  names(result)=c("Time","Delta","Covari","Resid","sim_stat")
  
  return(result)
}
#What_omni()

What_fform=function(b,std,Time,Delta,Covari,tol){
  #b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parametersa
  
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
    pi_i_z[[i]]=apply(Covari,2,function(x){(x<=((x[order(x)])[i]))*1})
  }
  pi_ij_z=list(NA)
  for(j in 1:p){
    pi_ij_z[[j]]=as.list(data.frame(t(matrix(unlist(
      lapply(pi_i_z,function(x){x[,j]})),nrow=n))))
  }
  
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
  
  E_t=S_1_t/S_0_t
  #E_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t 
  
  Mhat_i_t=mapply("-",N_i_t,lapply(lapply(
    Y_i_t,"*",dLambdahat_0_t),cumsum),SIMPLIFY=FALSE)
  #Mhat_i_t
  
  dN_i_t=lapply(N_i_t,function(x){diff(c(0,x))})
  #dN_i_t
  
  psi_t=S_0_t/n # Gehan's weight
  #psi_t
  
  U_t=Reduce('+',lapply(mapply("*",lapply(lapply(as.list(data.frame(t(Covari))),
      function(x){t(x-t(E_t))}),"*",psi_t),dN_i_t,SIMPLIFY=FALSE),function(x)
      {apply(x,2,cumsum)}))
  #U_t
  
  S_pi_t.z=lapply(lapply(pi_ij_z,function(z){mapply(function(x,y){x%*%t(y)},z,
      Y_i_t,SIMPLIFY=FALSE)}),function(x){Reduce('+',x)})
  #S_pi_t.z
  
  E_pi_t.z=lapply(S_pi_t.z,function(x){t(t(x)/S_0_t)})
  #E_pi_t.z
  
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
    ghat_t.z[[j]]=Reduce('+',lapply(mapply("*",pi_ij_z[[j]],Covari[,j],SIMPLIFY=FALSE),
                                    '%*%',t(ghat_0_t*Time)))/n
  }
  #ghat_t.z
  
  #-----------------------------f0----------------------------
  KM_e=cumprod(1-(Delta/S_0_t))
  #KM_e
  
  Fhat_0_e=1-KM_e
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
    fhat_t.z[[j]]=Reduce('+',lapply(mapply("*",pi_ij_z[[j]],Delta*Covari[,j],
                  SIMPLIFY=FALSE),'%*%',t(ghat_0_t*Time)))/n
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
    
    E_t_U=S_1_t_U/S_0_t_U
    #E_t_U
    
    psi_t_U=S_0_t_U/n # Gehan's weight
    #psi_t_U
    
    U_t_U=Reduce('+',lapply(mapply("*",lapply(lapply(as.list(data.frame(t(Covari_U))),
          function(x){t(x-t(E_t_U))}),"*",psi_t_U),dN_i_t_U,SIMPLIFY=FALSE),
                            function(x){apply(x,2,cumsum)}))
    #U_t_U
    
    U_inf_U=U_t_U[n,]
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  tolerance=tol+1 #initial value
  
  while (tolerance>tol){
    
    phi_i=rnorm(n)
    #phi_i
    
    U_pi_phi_t.z=list(NA)
    for(j in 1:p){
      U_pi_phi_t.z[[j]]=Reduce('+',mapply("*",lapply(mapply(function(x,y){t(t(x)*y)},
                                                            lapply(pi_ij_z[[j]],function(x,y){x-y},E_pi_t.z[[j]]),lapply(dMhat_i_t,"*",psi_t)
                                                            ,SIMPLIFY=FALSE),function(x){t(apply(x,1,cumsum))}),phi_i,SIMPLIFY=FALSE))
    }
    #U_pi_phi_t.z
    
    U_pi_phi_inf.z=lapply(U_pi_phi_t.z,function(x){x[,n]})
    #U_pi_phi_inf.z
    
    U_phi_t=list(NA)
    for(j in 1:p){
      U_phi_t[[j]]=Reduce('+',mapply("*",lapply(mapply("*",lapply(Covari[,j],'-',E_t[,j]),
                                                       lapply(dMhat_i_t,"*",psi_t),SIMPLIFY=FALSE),cumsum),phi_i,SIMPLIFY=FALSE))
    }
    #U_phi_t
    
    U_phi_inf=unlist(lapply(U_phi_t,function(x){x[n]}))
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
  }
  
  e_i_beta_s=as.vector(log(Time)+Covari%*%beta_hat_s)
  
  order_resid_s=order(e_i_beta_s)
  
  Time_s=Time[order_resid_s]
  Covari_s=matrix(Covari[order_resid_s,],nrow=n)
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
  
  S_1_t_s=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t_s,
                            as.list(data.frame(t(Covari_s))),SIMPLIFY=FALSE))
  #S_1_t
  
  E_t_s=S_1_t_s/S_0_t_s
  #E_t
  
  J_t_s=(S_0_t_s>0)*1
  #J_t_s
  
  dN_d_t_s=diff(c(0,N_d_t_s))
  #dN_d_t_s
  
  Lambdahat_0_t_s=cumsum((J_t_s/S_0_t_s)*dN_d_t_s)
  #Lambdahat_0_t_s
  
  dLambdahat_0_t_s=diff(c(0,Lambdahat_0_t_s))
  #dLambdahat_0_t_s
  
  F.T.=lapply(U_pi_phi_inf.z,function(x){x/sqrt(n)})
  S.T.=mapply("*",mapply("+",fhat_inf.z,mapply(function(x){apply(x,2,sum)},
        lapply(ghat_t.z,function(x,y){t(x)*y},dLambdahat_0_t),SIMPLIFY=FALSE),
        SIMPLIFY=FALSE),sqrt(n)*(b-beta_hat_s),SIMPLIFY=FALSE)
  T.T.=lapply(S_pi_t.z,function(x){apply(t(x)*diff(c(0,Lambdahat_0_t-Lambdahat_0_t_s))
        ,2,sum)/sqrt(n)})
  
  sim_stat=mapply(function(x,y,z){x-y-z},F.T.,S.T.,T.T.)
  #sim_stat
  
  result=list(Time,Delta,Covari,e_i_beta,sim_stat)
  names(result)=c("Time","Delta","Covari","Resid","sim_stat")
  
  return(result)}
#What_fform()

What_linkf=function(b,std,Time,Delta,Covari,tol){
  #b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parametersa
  
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
  
  E_t=S_1_t/S_0_t
  #E_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t 
  
  Mhat_i_t=mapply("-",N_i_t,lapply(lapply(
    Y_i_t,"*",dLambdahat_0_t),cumsum),SIMPLIFY=FALSE)
  #Mhat_i_t
  
  dN_i_t=lapply(N_i_t,function(x){diff(c(0,x))})
  #dN_i_t
  
  psi_t=S_0_t/n # Gehan's weight
  #psi_t
  
  U_t=Reduce('+',lapply(mapply("*",lapply(lapply(as.list(data.frame(t(Covari))),function(x)
        {t(x-t(E_t))}),"*",psi_t),dN_i_t,SIMPLIFY=FALSE),function(x){apply(x,2,cumsum)}))
  #U_t
  
  S_pi_t.z=Reduce('+',mapply(function(x,y){x%*%t(y)},pi_i_z,Y_i_t,SIMPLIFY=FALSE))
  #S_pi_t.z
  
  E_pi_t.z=t(t(S_pi_t.z)/S_0_t)
  #E_pi_t.z
  
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
    ghat_t.z[[j]]=Reduce('+',lapply(mapply("*",pi_i_z,Covari[,j],SIMPLIFY=FALSE),'%*%',t(ghat_0_t*Time)))/n
  }
  #ghat_t.z
  
  #-----------------------------f0----------------------------
  KM_e=cumprod(1-(Delta/S_0_t))
  #KM_e
  
  Fhat_0_e=1-KM_e
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
    fhat_t.z[[j]]=Reduce('+',lapply(mapply("*",pi_i_z,Delta*Covari[,j],SIMPLIFY=FALSE),'%*%',t(ghat_0_t*Time)))/n
  }
  #fhat_t.z
  
  fhat_inf.z=lapply(fhat_t.z,function(x){x[,n]})
  
  #-----------------------------------------------------------
  # ghat_t.z=list(NA)
  # for(j in 1:p){
  #   ghat_t.z[[j]]=Reduce('+',lapply(lapply(pi_i_z,"*",Covari[,j]),'%*%',t(ghat_0_t*Time)))/n
  # }
  # 
  # fhat_t.z=list(NA)
  # for(j in 1:p){
  #   fhat_t.z[[j]]=Reduce('+',lapply(lapply(pi_i_z,"*",Delta*Covari[,j]),'%*%',t(ghat_0_t*Time)))/n
  # }
  #-----------------------------------------------------------
  
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
    
    E_t_U=S_1_t_U/S_0_t_U
    #E_t_U
    
    psi_t_U=S_0_t_U/n # Gehan's weight
    #psi_t_U
    
    U_t_U=Reduce('+',lapply(mapply("*",lapply(lapply(as.list(data.frame(t(Covari_U))),
              function(x){t(x-t(E_t_U))}),"*",psi_t_U),dN_i_t_U,SIMPLIFY=FALSE),
              function(x){apply(x,2,cumsum)}))
    #U_t_U
    
    U_inf_U=U_t_U[n,]
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  tolerance=tol+1 #initial value
  
  while (tolerance>tol){
    
    phi_i=rnorm(n)
    #phi_i
    
    U_pi_phi_t.z=Reduce('+',mapply("*",lapply(mapply(function(x,y){t(t(x)*y)},
              lapply(pi_i_z,function(x,y){x-E_pi_t.z}),lapply(dMhat_i_t,"*",psi_t),
              SIMPLIFY=FALSE),function(x){t(apply(x,1,cumsum))}),phi_i,SIMPLIFY=FALSE))
    #U_pi_phi_t.z
    
    U_pi_phi_inf.z=U_pi_phi_t.z[,n]
    
    U_phi_t=Reduce('+',mapply("*",lapply(mapply("*",lapply(as.list(data.frame(
      t(Covari))),function(x){t(x-t(E_t))}),lapply(dMhat_i_t,"*",psi_t),
      SIMPLIFY=FALSE),function(x){apply(x,2,cumsum)}),phi_i,SIMPLIFY=FALSE))
    #U_phi_t
    
    U_phi_inf=U_phi_t[n,]
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
  }
  
  e_i_beta_s=as.vector(log(Time)+Covari%*%beta_hat_s)
  
  order_resid_s=order(e_i_beta_s)
  
  Time_s=Time[order_resid_s]
  Covari_s=matrix(Covari[order_resid_s,],nrow=n)
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
  
  S_1_t_s=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t_s,
                            as.list(data.frame(t(Covari_s))),SIMPLIFY=FALSE))
  #S_1_t
  
  E_t_s=S_1_t_s/S_0_t_s
  #E_t
  
  J_t_s=(S_0_t_s>0)*1
  #J_t_s
  
  dN_d_t_s=diff(c(0,N_d_t_s))
  #dN_d_t_s
  
  Lambdahat_0_t_s=cumsum((J_t_s/S_0_t_s)*dN_d_t_s)
  #Lambdahat_0_t_s
  
  dLambdahat_0_t_s=diff(c(0,Lambdahat_0_t_s))
  #dLambdahat_0_t_s
  
  F.T.=(1/sqrt(n))*U_pi_phi_inf.z
  S.T.=sqrt(n)*Reduce("+",mapply("*",mapply("+",fhat_inf.z,
        mapply(function(x){apply(x,2,sum)},lapply(ghat_t.z,function(x,y){t(x)*y},
        dLambdahat_0_t),SIMPLIFY=FALSE),SIMPLIFY=FALSE),(b-beta_hat_s),SIMPLIFY=FALSE))
  T.T.=(1/sqrt(n))*(apply((t(S_pi_t.z)*diff(c(0,Lambdahat_0_t-Lambdahat_0_t_s))),2,sum))
  
  sim_stat=F.T.-S.T.-T.T.
  #sim_stat
  
  result=list(Time,Delta,Covari,e_i_beta,sim_stat)
  names(result)=c("Time","Delta","Covari","Resid","sim_stat")
  
  return(result)
}
#What_linkf()

What_t.z=function(b,std,Time,Delta,Covari,test,tol){
  if(test=="omni"){
    return(What_omni(b,std,Time,Delta,Covari,tol))
  }
  if(test=="fform"){
    return(What_fform(b,std,Time,Delta,Covari,tol))
  }
  if(test=="linkf"){
    return(What_linkf(b,std,Time,Delta,Covari,tol))
  }
  if(test=="aft"){
    return(print("NOT YET..."))
    #return(What_aft(b,std,Time,Delta,Covari,tol))
  }
}
#What_t.z()

##############################################################
# omnibus test
# a=W_t.z(beta_hat_wb,X_wb,D_wb,Z_wb,"omni")$obs_stat[100,]
# for(i in 1:30){
#   plot(What_t.z(beta_hat_wb,std_hat_wb,X_wb,D_wb,Z_wb,"omni",
#   given_tol)$sim_stat[100,],type="s",col="grey",ylim=c(-2,2));par(new=TRUE)
# }
# plot(a,type="s",col="red",ylim=c(-2,2))

##############################################################
# functional form
# a=W_t.z(beta_hat_wb,X_wb,D_wb,Z_wb,"fform")$obs_stat[,1]
# for(i in 1:30){
#   plot(What_t.z(beta_hat_wb,std_hat_wb,X_wb,D_wb,Z_wb,"fform",
#   given_tol)$sim_stat[,1],type="s",col="grey",ylim=c(-2,2));par(new=TRUE)
# }
# plot(a,type="s",col="red",ylim=c(-2,2))

##############################################################
# link function
# a=W_t.z(beta_hat_wb,X_wb,D_wb,Z_wb,"linkf")$obs_stat
# for(i in 1:30){
#   plot(What_t.z(beta_hat_wb,std_hat_wb,X_wb,D_wb,Z_wb,"linkf",
#   given_tol)$sim_stat,type="s",col="grey",ylim=c(-2,2));par(new=TRUE)
# }
# plot(a,type="s",col="red",ylim=c(-2,2))

##############################################################

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
    dataset_What[[k]]=What_omni(b,std,Time,Delta,Covari,tol)$sim_stat
    if(k%%path_check==0) {
      cat("Sample Path",k,"\n")
    }
  }
  
  #------------------------BOOTSTRAPPING----------------------
  std.boot=matrix(apply(mapply(function(x){as.vector(x)},dataset_What),1,sd),nrow=n)
  # std.boot
  
  dataset_std.What=lapply(dataset_What,function(x){x/std.boot})
  # dataset_std.What
  
  dataset_W=W_omni(b,Time,Delta,Covari)$obs_stat
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
  
  #-----------------------------------------------------------
  #  p  : the ratio of (What>=W)*1
  # H_0 : the data follow the assumption of the aft model.
  #
  # if p>0.05, cannot reject the null hypothesis. i.e. accept it 
  # if p=<0.05, reject the null hypothesis.
  #
  # absolute/maximum 기준으로 What이 큰것의 비율(p)이 
  # 0.96이면 당연히 accetp
  # 0.04이면 당연히 reject
  # 0.45이면 accetp
  # p_alpha는 acceptance rate을 구하는 것이다! 
  #-----------------------------------------------------------
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_path_What>max_path_W)*1==1))/path
  # p_value
  
  std.p_value=length(which((max_path_std.What>max_path_std.W)*1==1))/path
  # std.p_value
  
  result=list(dataset_What=dataset_What,dataset_std.What=dataset_std.What,
              dataset_W=dataset_W,dataset_std.W=dataset_std.W,
              std.boot=std.boot,p_value=p_value,std.p_value=std.p_value)
  # result
  
  return(result)
}
#sample_path_omni

sample_path_fform=function(path,b,std,Time,Delta,Covari,tol){
  #path=path;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #path=path;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  path_check=ceiling(path/2)
  
  #------------------------SAMPLE PATH------------------------
  dataset_What=list(NA)
  for(k in 1:path){
    dataset_What[[k]]=What_fform(b,std,Time,Delta,Covari,tol)$sim_stat
    if(k%%path_check==0) {
      cat("Sample Path",k,"\n")
    }
  }
  
  #------------------------BOOTSTRAPPING----------------------
  std.boot=matrix(apply(mapply(function(x){as.vector(x)},dataset_What),1,sd),nrow=n)
  # std.boot
  
  dataset_std.What=lapply(dataset_What,function(x){x/std.boot})
  # dataset_std.What
  
  dataset_W=W_fform(b,Time,Delta,Covari)$obs_stat
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
  
  #-----------------------------------------------------------
  #  p  : the ratio of (What>=W)*1
  # H_0 : the data follow the assumption of the aft model.
  #
  # if p>0.05, cannot reject the null hypothesis. i.e. accept it 
  # if p=<0.05, reject the null hypothesis.
  #
  # absolute/maximum 기준으로 What이 큰것의 비율(p)이 
  # 0.96이면 당연히 accetp
  # 0.04이면 당연히 reject
  # 0.45이면 accetp
  # p_alpha는 acceptance rate을 구하는 것이다! 
  #-----------------------------------------------------------
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_path_What>max_path_W)*1==1))/path
  # p_value
  
  std.p_value=length(which((max_path_std.What>max_path_std.W)*1==1))/path
  # std.p_value
  
  result=list(dataset_What=dataset_What,dataset_std.What=dataset_std.What,
              dataset_W=dataset_W,dataset_std.W=dataset_std.W,
              std.boot=std.boot,p_value=p_value,std.p_value=std.p_value)
  # result
  
  return(result)
}
#sample_path_fform

sample_path_linkf=function(path,b,std,Time,Delta,Covari,tol){
  #path=path;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #path=path;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  path_check=ceiling(path/2)
  
  #------------------------SAMPLE PATH------------------------
  dataset_What=list(NA)
  for(k in 1:path){
    dataset_What[[k]]=What_linkf(b,std,Time,Delta,Covari,tol)$sim_stat
    if(k%%path_check==0) {
      cat("Sample Path",k,"\n")
    }
  }
  
  #------------------------BOOTSTRAPPING----------------------
  std.boot=matrix(apply(mapply(function(x){as.vector(x)},dataset_What),1,sd),nrow=n)
  # std.boot
  
  dataset_std.What=lapply(dataset_What,function(x){x/std.boot})
  # dataset_std.What
  
  dataset_W=W_linkf(b,Time,Delta,Covari)$obs_stat
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
  
  #-----------------------------------------------------------
  #  p  : the ratio of (What>=W)*1
  # H_0 : the data follow the assumption of the aft model.
  #
  # if p>0.05, cannot reject the null hypothesis. i.e. accept it 
  # if p=<0.05, reject the null hypothesis.
  #
  # absolute/maximum 기준으로 What이 큰것의 비율(p)이 
  # 0.96이면 당연히 accetp
  # 0.04이면 당연히 reject
  # 0.45이면 accetp
  # p_alpha는 acceptance rate을 구하는 것이다! 
  #-----------------------------------------------------------
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_path_What>max_path_W)*1==1))/path
  # p_value
  
  std.p_value=length(which((max_path_std.What>max_path_std.W)*1==1))/path
  # std.p_value
  
  result=list(dataset_What=dataset_What,dataset_std.What=dataset_std.What,
              dataset_W=dataset_W,dataset_std.W=dataset_std.W,
              std.boot=std.boot,p_value=p_value,std.p_value=std.p_value)
  # result
  
  return(result)
}
#sample_path_linkf

sample_path=function(path,b,std,Time,Delta,Covari,test,tol){
  if(test=="omni"){
    return(sample_path_omni(path,b,std,Time,Delta,Covari,tol))
  }
  if(test=="fform"){
    return(sample_path_fform(path,b,std,Time,Delta,Covari,tol))
  }
  if(test=="linkf"){
    return(sample_path_linkf(path,b,std,Time,Delta,Covari,tol))
  }
  if(test=="aft"){
    return(print("NOT YET..."))
    #return(sample_path_aft(path,b,std,Time,Delta,Covari,tol))
  }
}
#sample_path

##############################################################
# omnibus test
#####aft
# path1=30
# result_omni_aft=sample_path_omni(path1,beta_hat_ln_aft,std_hat_ln_aft,
# X_ln_aft,D_ln_aft,Z_ln_aft,0.1)
# for(i in 1:path1){
#   plot(result_omni_aft$dataset_std.What[[i]][(n/2),],ylim=c(-3,3),type="s",
# col="grey");par(new=TRUE)
# }
# plot(result_omni_aft$dataset_std.W[(n/2),],ylim=c(-3,3),type="s",col="red")

#####cox
# path1=30
# result_omni_cox=sample_path_omni(path1,beta_hat_ln_cox,std_hat_ln_cox,
# X_ln_cox,D_ln_cox,Z_ln_cox,0.1)
# for(i in 1:path1){
#   plot(result_omni_cox$dataset_std.What[[i]][(n/2),],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_omni_cox$dataset_std.W[(n/2),],ylim=c(-3,3),type="s",col="red")

##############################################################
# functional form
#####aft
# path1=30
# result_fform_aft=sample_path_fform(path1,beta_hat_ln_aft,std_hat_ln_aft,
# X_ln_aft,D_ln_aft,Z_ln_aft,0.1)
# for(i in 1:path1){
#   plot(result_fform_aft$dataset_std.What[[i]][,1],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_fform_aft$dataset_std.W[,1],ylim=c(-3,3),type="s",col="red")
#####cox
# path1=30
# result_fform_cox=sample_path_fform(path1,beta_hat_ln_cox,std_hat_ln_cox,
# X_ln_cox,D_ln_cox,Z_ln_cox,0.1)
# for(i in 1:path1){
#   plot(result_fform_cox$dataset_std.What[[i]][,1],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_fform_cox$dataset_std.W[,1],ylim=c(-3,3),type="s",col="red")

##############################################################
# link function
#####aft
# path1=30
# result_linkf_aft=sample_path_linkf(path1,beta_hat_ln_aft,std_hat_ln_aft,
# X_ln_aft,D_ln_aft,Z_ln_aft,0.1)
# for(i in 1:path1){
#   plot(result_linkf_aft$dataset_std.What[[i]],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_linkf_aft$dataset_std.W,ylim=c(-3,3),type="s",col="red")
#####cox
# path1=30
# result_linkf_cox=sample_path_linkf(path1,beta_hat_ln_cox,std_hat_ln_cox,
# X_ln_cox,D_ln_cox,Z_ln_cox,0.1)
# for(i in 1:path1){
#   plot(result_linkf_cox$dataset_std.What[[i]],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_linkf_cox$dataset_std.W,ylim=c(-3,3),type="s",col="red")

#-------------------------------------------------------------
#---------------------------PLOTTING--------------------------
#-------------------------------------------------------------
plotting_omni=function(result,path){}

plotting_std.omni=function(result,path){}

plotting_fform=function(result,path){}

plotting_std.fform=function(result,path){}

plotting_linkf=function(result,path){}

plotting_std.linkf=function(result,path){}


plotting_asdfasdfasdf=function(result,standardization="standardized",path=50){
  
  if (standardization==0) {
    dataset_What=data.frame()
    for (i in 1:n.path){
      group=i
      A=result$dataset_What[,i]
      AA=data.frame(group,t_i=1:n,What=A)
      dataset_What=rbind(dataset_What,AA)
    }
    #dataset_What
    
    dataset_W=data.frame(group,t_i=1:n,W=result$dataset_W)
    #dataset_W
    
    Figure1_W=
      ggplot()+
      geom_step(data=dataset_What,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
      geom_step(data=dataset_W,aes(x=t_i,y=W),colour="tomato")
    #Figure1_W
    
    return(Figure1_W)
  }
  if (standardization==1) {
    dataset_std.What=data.frame()
    for (i in 1:n.path){
      group=i
      A=result$dataset_std.What[,i]
      AA=data.frame(group,t_i=1:n,std.What=A)
      dataset_std.What=rbind(dataset_std.What,AA)
    }
    #dataset_std.What
    
    dataset_std.W=data.frame(group,t_i=1:n,std.W=result$dataset_std.W)
    #dataset_std.W
    
    Figure1_std.W=
      ggplot()+
      geom_step(data=dataset_std.What,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
      geom_step(data=dataset_std.W,aes(x=t_i,y=std.W),colour="tomato")
    return(Figure1_std.W)
  }
}

plotting=function(result,standardization="standardized",path=50){
  if(test=="omni"){
    if(standardization=="standardized"){
      return(plotting_std.omni(result,path=50))
    }
    if(standardization=="unstandardized"){
      return(plotting_omni(result,path=50))
    }
  }
  if(test=="fform"){
    if(standardization=="standardized"){
      return(plotting_std.fform(result,path=50))
    }
    if(standardization=="unstandardized"){
      return(plotting_fform(result,path=50))
    }
  }
  if(test=="linkf"){
    if(standardization=="standardized"){
      return(plotting_std.linkf(result,path=50))
    }
    if(standardization=="unstandardized"){
      return(plotting_linkf(result,path=50))
    }
  }
  # if(test=="aft"){
  #   if(standardization=="standardized"){
  #     return(plotting_std.aft(result,path=50))
  #   }
  #   if(standardization=="unstandardized"){
  #     return(plotting_aft(result,path=50))
  #   }
  # }
}
#plotting

##############################################################
##############################################################
##############################################################
############################미완성############################
##############################################################
##############################################################
##############################################################
W_aft=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  pi_ij_z=list(NA)
  for(j in 1:p){
    pi_ij_z[[j]]=Covari[,j]
  }
  
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
  
  E_t=S_1_t/S_0_t
  #E_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t
  
  Mhat_i_t=mapply("-", N_i_t,lapply(lapply(
    Y_i_t,"*",dLambdahat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  pi_ij_z.Mhat_i_t=list(NA)
  for(j in 1:p){
    pi_ij_z.Mhat_i_t[[j]]=lapply(Mhat_i_t,function(x){pi_ij_z[[j]]%*%t(x)})
  }
  #pi_ij_z.Mhat_i_t
  # z by t matrix
  
  obs_stat=list(NA)
  for(j in 1:p){
    obs_stat[[j]]=Reduce('+',pi_ij_z.Mhat_i_t[[j]])/sqrt(n)
  }
  #obs_stat
  
  result=list(Time,Delta,Covari,e_i_beta,obs_stat)
  names(result)=c("Time","Delta","Covari","Resid","obs_stat")
  
  return(result)
}
#W_aft()

What_aft=function(b,std,Time,Delta,Covari,tol){}
#What_aft()

sample_path_aft=function(path,b,std,Time,Delta,Covari,tol){}
# result_aft=sample_path_aft(1,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,0.1)

plotting_aft=function(result,path){}

plotting_std.aft=function(result,path){}

