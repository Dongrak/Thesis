W_t_aft=function(b,Time,Delta,Covari,weight="f1"){
  #b=beta_hat_cox;Time=T_cox;Delta=D_cox;Covari=Z_cox;weight="f1"
  #b=beta_hat_aft;Time=T_aft;Delta=D_aft;Covari=Z_aft;weight="f1"
  
  # Covari is n by J matrix consited of the covariates
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time) # the number of individuals
  p=length(b) # the number of parameters
  
  e_i_beta=as.vector(log(Time)-Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  if(p==1){Covari=Covari[order_resid]}
  if(p>1){Covari=Covari[order_resid,]}
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  if(weight=="11"){
    w_ij_z=1
  }
  
  if(p==1){
    if(weight=="fi"){
      w_ij_z=list(NA)
      w_i_z=list(NA)
      Covari_order_j=as.vector(Covari)[order(Covari)]
      for(i in 1:n){
        w_i_z[[i]]=((Covari_order_j>=Covari_order_j[i])*Covari_order_j[i])[order_resid]
      }
      w_ij_z[[p]]=w_i_z
    }
    if(weight=="f1"){
      w_ij_z=list(NA)
      w_ij_z[[p]]=Covari
    }
    if(weight=="1i"){
      w_ij_z=list(NA)
      w_i_z=list(NA)
      Covari_order_j=as.vector(Covari)[order(Covari)]
      for(i in 1:n){
        w_i_z[[i]]=((Covari_order_j>=Covari_order_j[i])*1)[order_resid]
      }
      w_ij_z[[p]]=w_i_z
    }
  }
  
  if(p>1){
    if(weight=="fi"){
      w_ij_z=list(NA)
      w_i_z=list(NA)
      for(j in 1:p){
        Covari_order_j=as.vector(Covari[,j])[order(Covari[,j])]
        for(i in 1:n){
          w_i_z[[i]]=((Covari_order_j>=Covari_order_j[i])*Covari_order_j[i])[order_resid]
        }
        w_ij_z[[j]]=w_i_z
      }
      #w_ij_z
    }
    if(weight=="f1"){
      w_ij_z=list(NA)
      for(j in 1:p){
        w_ij_z[[j]]=Covari[,j]
      }
    }
    if(weight=="1i"){
      w_ij_z=list(NA)
      w_i_z=list(NA)
      for(j in 1:p){
        Covari_order_j=as.vector(Covari[,j])[order(Covari[,j])]
        for(i in 1:n){
          w_i_z[[i]]=((Covari_order_j>=Covari_order_j[i])*1)[order_resid]
        }
        w_ij_z[[j]]=w_i_z
      }
    }
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
  
  Ahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Ahat_0_t
  
  dAhat_0_t=diff(c(0,Ahat_0_t))
  #dAhat_0_t
  
  Mhat_i_t=mapply("-", N_i_t,lapply(lapply(
    Y_i_t,"*",dAhat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  w_ij_z.Mhat_i_t=list(NA)
  for(j in 1:p){
    w_ij_z.Mhat_i_t[[j]]=mlapply(list(w_ij_z[[j]],Mhat_i_t),function(x,y){x%*%t(y)})
  }
  #w_ij_z.Mhat_i_t
  # z by t matrix
  
  W_j_t.z=list(NA)
  for(j in 1:p){
    W_j_t.z[[j]]=Reduce('+',w_ij_z.Mhat_i_t[[j]])
  }
  #W_j_t.z
  
  result=list(Time,Delta,Covari,e_i_beta,W_j_t.z)
  names(result)=c("Time","Delta","Covari","Resid","W_j_t.z")
  
  return(result)
}
#W_t.z_omni()
W_t_aft(beta_hat_aft,T_aft,D_aft,Z_aft,"f1")
