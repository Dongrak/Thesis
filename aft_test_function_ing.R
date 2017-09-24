#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------
W_omni=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time)
  p=length(b)
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  pi_ij_z=list(NA)
  pi_i_z=list(NA)
  for(j in 1:p){
    Covari_order_j=as.vector(Covari[,j])[order(Covari[,j])]
    for(i in 1:n){
      pi_i_z[[i]]=((Covari_order_j>=Covari_order_j[i]))[order_resid]
    }
    pi_ij_z[[j]]=pi_i_z
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
  
  S_1_t=Reduce('+',mapply(
    "*", Y_i_t, Covari, SIMPLIFY = FALSE))
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
    pi_ij_z.Mhat_i_t[[j]]=mlapply(list(pi_ij_z[[j]],Mhat_i_t),function(x,y){x%*%t(y)})
  }
  #pi_ij_z.Mhat_i_t
  # z by t matrix
  
  obs.stat=list(NA)
  for(j in 1:p){
    obs.stat[[j]]=Reduce('+',pi_ij_z.Mhat_i_t[[j]])/sqrt(n)
  }
  #obs.stat
  
  result=list(Time,Delta,Covari,e_i_beta,obs.stat)
  names(result)=c("Time","Delta","Covari","Resid","obs.stat")
  
  return(result)
}
#W_omni()

W_ftnform=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time)
  p=length(b)
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  pi_ij_z=list(NA)
  pi_i_z=list(NA)
  for(j in 1:p){
    Covari_order_j=as.vector(Covari[,j])[order(Covari[,j])]
    for(i in 1:n){
      pi_i_z[[i]]=((Covari_order_j>=Covari_order_j[i])*1)[order_resid]
    }
    pi_ij_z[[j]]=pi_i_z
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
  
  S_1_t=Reduce('+',mapply(
    "*", Y_i_t, Covari, SIMPLIFY = FALSE))
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
  
  pi_ij_z.Mhat_i_inf=list(NA)
  for(j in 1:p){
    pi_ij_z.Mhat_i_inf[[j]]=mapply("*",pi_ij_z[[j]],Mhat_i_inf,SIMPLIFY = FALSE)
  }
  #pi_ij_z.Mhat_i_inf
  # z by 1 matrix
  
  obs.stat=lapply(pi_ij_z.Mhat_i_inf,function(x){Reduce('+',x)/sqrt(n)})
  #obs.stat
  
  result=list(Time,Delta,Covari,e_i_beta,obs.stat)
  names(result)=c("Time","Delta","Covari","Resid","obs.stat")
  
  return(result)
}
#W_ftnform()

W_linkftn=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time)
  p=length(b)
  
  if(p==1){return(print("ERROR MESSAGE : the number  needs to be greater than one"))}
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  pi_ij_z=list(NA)
  pi_i_z=list(NA)
  for(j in 1:p){
    Covari_order_j=as.vector(Covari[,j])[order(Covari[,j])]
    for(i in 1:n){
      pi_i_z[[i]]=((Covari_order_j>=Covari_order_j[i])*1)[order_resid]
    }
    pi_ij_z[[j]]=pi_i_z
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
  
  S_1_t=Reduce('+',mapply(
    "*", Y_i_t, Covari, SIMPLIFY = FALSE))
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
  
  pi_ij_z.Mhat_i_inf=list(NA)
  for(j in 1:p){
    pi_ij_z.Mhat_i_inf[[j]]=mapply("*",pi_ij_z[[j]],Mhat_i_inf,SIMPLIFY = FALSE)
  }
  #pi_ij_z.Mhat_i_inf
  # z by 1 matrix
  
  obs.stat=lapply(pi_ij_z.Mhat_i_inf,function(x){Reduce('+',x)/sqrt(n)})
  #obs.stat
  
  result=list(Time,Delta,Covari,e_i_beta,obs.stat)
  names(result)=c("Time","Delta","Covari","Resid","obs.stat")
  
  return(result)
}
#W_linkftn()

################################################ mi wan sung
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
  
  S_1_t=Reduce('+',mapply(
    "*", Y_i_t, Covari, SIMPLIFY = FALSE))
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
  
  obs.stat=list(NA)
  for(j in 1:p){
    obs.stat[[j]]=Reduce('+',pi_ij_z.Mhat_i_t[[j]])/sqrt(n)
  }
  #obs.stat
  
  result=list(Time,Delta,Covari,e_i_beta,obs.stat)
  names(result)=c("Time","Delta","Covari","Resid","obs.stat")
  
  return(result)
}
#W_aft()

W_t.z=function(b,Time,Delta,Covari,test){
  if(test=="omni"){
    return(W_omni(b,Time,Delta,Covari))
  }
  if(test=="ftnform"){
    return(W_ftnform(b,Time,Delta,Covari))
  }
  if(test=="linkftn"){
    return(W_linkftn(b,Time,Delta,Covari))
  }
  if(test=="aft"){
    return(print("NOT YET..."))
    #return(W_aft(b,Time,Delta,Covari))
  }
}
#W_t.z()
W_t.z(beta_hat_wb,X_wb,D_wb,Z_wb,"omni")

#-------------------------------------------------------------
#-------------------------REALIZATION-------------------------
#-------------------------------------------------------------
What_omni=function(b,std,Time,Delta,Covari,tol){
  #b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  b=c(beta_hat_wb,beta_hat_wb^2);Covari=c(Z_wb,Z_wb^2);
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time) # the number of individuals
  p=length(b) # the number of parametersa
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  pi_ij_z=list(NA)
  pi_i_z=list(NA)
  for(j in 1:p){
    Covari_order_j=as.vector(Covari[,j])[order(Covari[,j])]
    for(i in 1:n){
      pi_i_z[[i]]=((Covari_order_j>=Covari_order_j[i])*1)[order_resid]
    }
    pi_ij_z[[j]]=pi_i_z
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
  
  S_j_1_t=list(NA)
  for(j in 1:p){
    S_j_1_t[[j]]=Reduce('+',mapply(
      "*", Y_i_t, Covari[,j], SIMPLIFY = FALSE))
  }
  #S_j_1_t
  
  E_j_t=list(NA)
  for(j in 1:p){
    E_j_t[[j]]=S_j_1_t[[j]]/S_0_t
  }
  #E_j_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t 
  
  Mhat_i_t=mapply("-",N_i_t,lapply(lapply(
    Y_i_t,"*",dLambdahat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  dN_i_t=lapply(N_i_t,function(x){diff(c(0,x))})
  #dN_i_t
  
  psi_t=S_0_t/n # Gehan's weight
  #psi_t
  
  U_j_t=list(NA)
  for(j in 1:p){
    U_j_t[[j]]=Reduce('+',lapply(mapply("*",dN_i_t,lapply(
      lapply(Covari[,j],'-',E_j_t[[j]]),"*",psi_t), SIMPLIFY = FALSE),cumsum))
  }
  #U_j_t
  
  pi_ij_z.Y_i_t=list(NA)
  for(j in 1:p){
    pi_ij_z.Y_i_t[[j]]=mlapply(list(pi_ij_z[[j]],Y_i_t),function(x,y){x%*%t(y)})
  }
  #pi_ij_z.Y_i_t
  
  S_j_pi_t.z=list(NA)
  for(j in 1:p){
    S_j_pi_t.z[[j]]=Reduce('+',pi_ij_z.Y_i_t[[j]])
  }
  #S_j_pi_t.z
  
  E_j_pi_t.z=lapply(S_j_pi_t.z,function(x,y){t(t(x)/y)},S_0_t)
  #E_j_pi_t.z
  
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
  
  ghat_j_t.z=list(NA)
  for(j in 1:p){
    ghat_j_t.z[[j]]=Reduce('+',lapply(lapply(pi_ij_z[[j]],"*",Covari[,j]),'%*%',t(ghat_0_t*Time)))/n
  }
  #ghat_j_t.z
  
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
  
  fhat_j_t.z=list(NA)
  for(j in 1:p){
    fhat_j_t.z[[j]]=Reduce('+',lapply(lapply(pi_ij_z[[j]],"*",Delta*Covari[,j]),'%*%',t(ghat_0_t*Time)))/n
  }
  #fhat_j_t.z
  
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
    
    Y_i_t_U=list(NA)
    for(j in 1:n){
      Y_i_t_U[[j]]=(e_i_beta_U<=e_i_beta_U[j])*1
    }
    #Y_i_t_U
    
    N_d_t_U=Reduce('+',N_i_t_U)
    #N_d_t_U
    
    S_0_t_U=Reduce('+',Y_i_t_U)
    #S_0_t_U
    
    S_j_1_t_U=list(NA)
    for(j in 1:p){
      S_j_1_t_U[[j]]=Reduce('+',mapply(
        "*", Y_i_t_U, Covari_U[,j], SIMPLIFY = FALSE))
    }
    #S_j_1_t_U
    
    E_j_t_U=list(NA)
    for(j in 1:p){
      E_j_t_U[[j]]=S_j_1_t_U[[j]]/S_0_t_U
    }
    #E_j_t_U
    
    J_t_U=(S_0_t_U>0)*1
    #J_t_U
    
    dN_d_t_U=diff(c(0,N_d_t_U))
    #dN_d_t_U
    
    Lambdahat_0_t_U=cumsum((J_t_U/S_0_t_U)*dN_d_t_U)
    #Lambdahat_0_t_U
    
    dLambdahat_0_t_U=diff(c(0,Lambdahat_0_t_U))
    #dLambdahat_0_t _U
    
    Mhat_i_t=mapply("-", N_i_t_U,lapply(
      lapply(Y_i_t_U,"*",dLambdahat_0_t_U),cumsum), SIMPLIFY = FALSE)
    #Mhat_i_t _U
    
    dN_i_t_U=lapply(N_i_t_U,function(x){diff(c(0,x))})
    #dN_i_t_U
    
    psi_t_U=S_0_t_U/n # Gehan's weight
    #psi_t_U
    
    U_j_t_U=list(NA)
    for(j in 1:p){
      U_j_t_U[[j]]=Reduce('+',lapply(mapply("*",dN_i_t_U,lapply(
        lapply(Covari_U[,j],'-',E_j_t_U[[j]]),"*",psi_t_U), SIMPLIFY = FALSE)
        ,cumsum))
    }
    #U_j_t_U
    
    U_j_inf_U=unlist(lapply(U_j_t_U,function(x){x[n]}))
    #U_j_inf_U
    
    return(U_j_inf_U)
  }
  #U_beta()
  
  tolerance=tol+1 #initial value
  
  #while (tolerance>tol){
  
  phi_i=rnorm(n)
  #phi_i
  
  U_j_pi_phi_t.z=list(NA)
  for(j in 1:p){
    U_j_pi_phi_t.z[[j]]=Reduce('+',mapply(function(x){t(apply(x,2,cumsum))}
                                          ,mapply("*",mapply("*",dMhat_i_t,lapply(lapply(
                                            pi_ij_z[[j]],function(x,y){x-y},E_j_pi_t.z[[j]]),
                                            function(x,y){t(t(x)*y)},psi_t),SIMPLIFY=FALSE)
                                            ,phi_i,SIMPLIFY=FALSE),SIMPLIFY=FALSE))
  }
  
  U_j_phi_t=list(NA)
  for(j in 1:p){
    U_j_phi_t[[j]]=Reduce('+',lapply(mapply("*",mapply(
      "*",dMhat_i_t,lapply(lapply(Covari[,j],'-',E_j_t[[j]]),"*",psi_t)
      ,SIMPLIFY = FALSE),phi_i,SIMPLIFY = FALSE),cumsum))
  } 
  #U_j_phi_t
  
  U_j_phi_inf=unlist(lapply(U_j_phi_t,max))
  #U_j_phi_inf
  if(p==1){
    beta_hat_s_list=optimize(function(BETA){sum((U_beta(BETA)-U_j_phi_inf)^2)},
                             c(b-2*std,b+2*std),tol = 1e-16)
    #beta_hat_s_list
    
    beta_hat_s=beta_hat_s_list$minimum;beta_hat_s
    #beta_hat_s
    
    tolerance=beta_hat_s_list$objective
    #tolerance
    
  }
  if(p>1){
    beta_hat_s_list=optim(b,function(BETA){sum((U_beta(BETA)-U_j_phi_inf)^2)})
    #beta_hat_s_list
    
    beta_hat_s=beta_hat_s_list$par;beta_hat_s
    #beta_hat_s
    
    tolerance=beta_hat_s_list$objective
    #tolerance
  }
  #beta_hat_s_list
  #}
  
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
  
  S_1_t_s=Reduce('+',mapply("*", Y_i_t_s, Covari_s,
                            SIMPLIFY = FALSE))
  #S_1_t_s
  
  E_t_s=S_1_t_s/S_0_t_s
  #E_t_s
  
  J_t_s=(S_0_t_s>0)*1
  #J_t_s
  
  dN_d_t_s=diff(c(0,N_d_t_s))
  #dN_d_t_s
  
  Lambdahat_0_t_s=cumsum((J_t_s/S_0_t_s)*dN_d_t_s)
  #Lambdahat_0_t_s
  
  dLambdahat_0_t_s=diff(c(0,Lambdahat_0_t_s))
  #dLambdahat_0_t_s
  
  order_Time=order(Time)
  #order_Time
  
  F.T.=lapply(U_j_pi_phi_t.z,"/",sqrt(n))
  S.T.=mapply("*",mapply("+",fhat_j_t.z,mapply(function(x){t(apply(x,2,cumsum))}
                                               ,lapply(ghat_j_t.z,function(x,y){t(t(x)*y)},dLambdahat_0_t),
                                               SIMPLIFY=FALSE),SIMPLIFY=FALSE),(b-beta_hat_s)*sqrt(n),SIMPLIFY=FALSE)
  T.T.=lapply(mapply(function(x){t(apply(x,2,cumsum))},lapply(S_j_pi_t.z,function(x,y)
  {t(t(x)*y)},diff(c(0,Lambdahat_0_t-Lambdahat_0_t_s))),SIMPLIFY=FALSE),"/",sqrt(n))
  
  sim_stat=mapply("-",mapply("-",F.T.,S.T.,SIMPLIFY=FALSE),T.T.,SIMPLIFY=FALSE)
  #sim_stat
  
  result=list(Time,Delta,Covari,e_i_beta,sim_stat)
  names(result)=c("Time","Delta","Covari","Resid","sim_stat")
  
  return(result)
}
#What_omni()

What_ftnform=function(b,std,Time,Delta,Covari,tol){}
#What_ftnform()

What_linkftn=function(b,std,Time,Delta,Covari,tol){}
#What_linkftn()

What_aft=function(b,std,Time,Delta,Covari,tol){}
#What_aft()

What_t.z=function(b,std,Time,Delta,Covari,test,tol){
  if(test=="omni"){
    return(What_omni(b,std,Time,Delta,Covari,tol))
  }
  if(test=="ftnform"){
    return(What_ftnform(b,std,Time,Delta,Covari,tol))
  }
  if(test=="linkftn"){
    return(What_linkftn(b,std,Time,Delta,Covari,tol))
  }
  if(test=="aft"){
    return(print("NOT YET..."))
    #return(What_aft(b,std,Time,Delta,Covari,tol))
  }
}
#What_t.z()

#-------------------------------------------------------------
#-------------------------SAMPLE PATH-------------------------
#-------------------------------------------------------------
sample_path_omni=function(path,b,std,Time,Delta,Covari,tol){
  #path=path;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #path=path;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  
  p=length(b)
  
  path=200
  
  #------------------------SAMPLE PATH------------------------
  dataset_What=rep(list(list(NA)),p)
  for(k in 1:path){
    stat=What_omni(b,std,Time,Delta,Covari,tol)$sim_stat
    for(j in 1:p){
      dataset_What[[j]][[k]]=stat[[j]]
    }
    if(k%%100==0) {
      cat("Sample Path",k,"\n")
    }
  }
  
  #------------------------BOOTSTRAPPING----------------------
  std.boot=list(NA)
  for(j in 1:p){
    std.boot[[j]]=matrix((apply(mapply(function(x){as.vector(x)},dataset_What[[j]]),1,sd)),nrow=n)
  }
  # std.boot
  
  dataset_std.What=list(NA)
  for(j in 1:p){
    dataset_std.What[[j]]=lapply(dataset_What[[j]],function(x){x/std.boot[[j]]})
  }
  # dataset_std.What
  
  dataset_W=W_omni(b,Time,Delta,Covari)$obs.stat
  # dataset_W
  
  dataset_std.W=list(NA)
  for(j in 1:p){
    dataset_std.W[[j]]=dataset_W[[j]]/std.boot[[j]]
  }
  # dataset_std.W
  
  #-----------------------MAXIMUM VALUE-----------------------
  max_path_What=list(NA)
  for(j in 1:p){
    max_path_What[[j]]=Reduce("pmax",lapply(dataset_What[[j]],function(x){abs(x)}))
  }
  # max_path_What
  
  max_path_W=list(NA)
  for(j in 1:p){
    max_path_W[[j]]=max(abs(dataset_W[[j]]))
  }
  # max_path_W
  
  max_path_std.What=list(NA)
  for(j in 1:p){
    max_path_std.What[[j]]=Reduce("pmax",lapply(dataset_std.What[[j]],function(x){abs(x)}))
  }
  # max_path_std.What
  
  max_path_std.W=list(NA)
  for(j in 1:p){
    max_path_std.W[[j]]=max(abs(dataset_std.W[[j]]))
  }
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
  p_value=list(NA)
  for(j in 1:p){
    p_value[[j]]=length(which((max_path_What[[j]]>max_path_W[[j]])*1==1))/(n^2)
  }
  # p_value
  
  std.p_value=list(NA)
  for(j in 1:p){
    std.p_value[[j]]=length(which((max_path_std.What[[j]]>max_path_std.W[[j]])*1==1))/(n^2)
  }
  # std.p_value
  
  result=list(dataset_What=dataset_What,dataset_std.What=dataset_std.What,
              dataset_W=dataset_W,dataset_std.W=dataset_std.W,
              std.boot=std.boot,p_value=p_value,std.p_value=std.p_value)
  # result
  
  return(result)
}
# sample_path_omni(1,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,100000)

sample_path_ftnform=function(path,b,std,Time,Delta,Covari,tol){}
# sample_path_ftnform(1,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,100000)

sample_path_linkftn=function(path,b,std,Time,Delta,Covari,tol){}
# sample_path_linkftn(1,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,100000)

sample_path_aft=function(path,b,std,Time,Delta,Covari,tol){}
# sample_path_aft(1,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,100000)

sample_path=function(path,b,std,Time,Delta,Covari,test,tol){
  if(test=="omni"){
    return(sample_path_omni(path,b,std,Time,Delta,Covari,tol))
  }
  if(test=="ftnform"){
    return(sample_path_ftnform(path,b,std,Time,Delta,Covari,tol))
  }
  if(test=="linkftn"){
    return(sample_path_linkftn(path,b,std,Time,Delta,Covari,tol))
  }
  if(test=="aft"){
    return(print("NOT YET..."))
    #return(sample_path_aft(path,b,std,Time,Delta,Covari,tol))
  }
}
# sample_path(1,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,"omni",100000)

