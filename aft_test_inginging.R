What_j_t.z_omni=function(b,std,Time,Delta,Covari,tol){
  #b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(beta_hat_wb,beta_hat_wb^2);Covari=c(Z_wb,Z_wb^2);
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time) # the number of individuals
  p=length(b) # the number of parameters
  
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
  optim_U_beta=function(x){
    beta_U=x
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

    U_j_pi_phi_t.z1=list(NA)
    for(j in 1:p){
      U_j_pi_phi_t.z1[[j]]=Reduce('+',mapply(function(x){t(apply(x,2,cumsum))}
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
  
  What_j_t.z=mapply("-",mapply("-",F.T.,S.T.,SIMPLIFY=FALSE),T.T.,SIMPLIFY=FALSE)
  #What_j_t.z
  
  result=list(Time,Delta,Covari,e_i_beta,What_j_t.z)
  names(result)=c("Time","Delta","Covari","Resid","What_j_t.z")
  
  return(result)
}
#What_j_t.z_omni()
What_j_t.z_omni(c(beta_hat_wb,beta_hat_wb^2),std_hat_wb,X_wb,D_wb,c(Z_wb,Z_wb^2),given_tol)$What_j_t.z

