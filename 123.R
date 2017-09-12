W_t.z_wb=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  
  # Covari is n by J matrix consited of the covariates
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time) # the number of individuals
  p=length(b) # the number of parameters
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  if(p==1){Covari=Covari[order_resid]}
  if(p>1){Covari=Covari[order_resid,]}
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  if(p==1){
    w_ij_z=list(NA)
    w_ij_z[[p]]=Covari
  }
  if(p>1){
    w_ij_z=list(NA)
    for(j in 1:p){
      w_ij_z[[j]]=Covari[,j]
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
    w_ij_z.Mhat_i_t[[j]]=lapply(Mhat_i_t,function(x){w_ij_z[[j]]%*%t(x)})
  }
  #w_ij_z.Mhat_i_t
  # z by t matrix

  W_j_t.z=list(NA)
  for(j in 1:p){
    W_j_t.z[[j]]=Reduce('+',w_ij_z.Mhat_i_t[[j]])/sqrt(n)
  }
  #W_j_t.z
  
  ##############U_t
  
  result=list(Time,Delta,Covari,e_i_beta,W_j_t.z)
  names(result)=c("Time","Delta","Covari","Resid","W_j_t.z")
  
  return(result)
}
#W_t.z_wb()
W_t.z_wb(beta_hat_wb,T_wb,D_wb,Z_wb)
