#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------
W_j_t.z_omni=function(b,Time,Delta,Covari){
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
  
  W_j_t.z=list(NA)
  for(j in 1:p){
    W_j_t.z[[j]]=Reduce('+',pi_ij_z.Mhat_i_t[[j]])/sqrt(n)
  }
  #W_j_t.z
  
  result=list(Time,Delta,Covari,e_i_beta,W_j_t.z)
  names(result)=c("Time","Delta","Covari","Resid","W_j_t.z")
  
  return(result)
}
#W_j_t.z_omni()

W_j_z_ftnform=function(b,Time,Delta,Covari){
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
  
  W_j_inf.z=lapply(pi_ij_z.Mhat_i_inf,function(x){Reduce('+',x)/sqrt(n)})
  #W_j_inf.z
  
  result=list(Time,Delta,Covari,e_i_beta,W_j_inf.z)
  names(result)=c("Time","Delta","Covari","Resid","W_j_inf.z")
  
  return(result)
}
#W_j_z_ftnform()

W_j_z_linkftn=function(b,Time,Delta,Covari){
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
  
  W_j_inf.z=lapply(pi_ij_z.Mhat_i_inf,function(x){Reduce('+',x)/sqrt(n)})
  #W_j_inf.z
  
  result=list(Time,Delta,Covari,e_i_beta,W_j_inf.z)
  names(result)=c("Time","Delta","Covari","Resid","W_j_inf.z")
  
  return(result)
}
#W_j_z_linkftn()

################################################ mi wan sung
W_j_t.z_aft=function(b,Time,Delta,Covari){
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
  
  W_j_t.z=list(NA)
  for(j in 1:p){
    W_j_t.z[[j]]=Reduce('+',pi_ij_z.Mhat_i_t[[j]])/sqrt(n)
  }
  #W_j_t.z
  
  ##############U_t
  
  result=list(Time,Delta,Covari,e_i_beta,W_j_t.z)
  names(result)=c("Time","Delta","Covari","Resid","W_j_t.z")
  
  return(result)
}
#W_j_t.z_aft()

W_j_t.z=function(b,Time,Delta,Covari,test){
  if(test=="omni"){
    return(W_j_t.z_omni(b,Time,Delta,Covari))
  }
  if(test=="ftnform"){
    return(W_j_z_ftnform(b,Time,Delta,Covari))
  }
  if(test=="linkftn"){
    return(W_j_z_linkftn(b,Time,Delta,Covari))
  }
  if(test=="aft"){
    return(print("NOT YET..."))
    #return(W_t.z_aft(b,Time,Delta,Covari))
  }
}
#W_j_t.z()
W_j_t.z(beta_hat_wb,X_wb,D_wb,Z_wb,"omni")

#-------------------------------------------------------------
#-------------------------SAMPLE PATH-------------------------
#-------------------------------------------------------------
What_t.z_omni=function(b,std,Time,Delta,Covari,tol){}
#What_t.z_omni()

What_z_ftnform=function(b,std,Time,Delta,Covari,tol){}
#What_z_ftnform()

What_z_linkftn=function(b,std,Time,Delta,Covari,tol){}
#What_z_linkftn()

################################################ mi wan sung
What_t.z_aft=function(b,std,Time,Delta,Covari,tol){}
#What_t.z_aft()

What_t.z=function(b,std,Time,Delta,Covari,test,tol){
  if(test=="omni"){
    return(What_t.z_omni(b,std,Time,Delta,Covari,tol))
  }
  if(test=="ftnform"){
    return(What_z_ftnform(b,std,Time,Delta,Covari,tol))
  }
  if(test=="linkftn"){
    return(What_z_linkftn(b,std,Time,Delta,Covari,tol))
  }
  if(test=="aft"){
    return(print("NOT YET..."))
    #return(What_t.z_aft(b,std,Time,Delta,Covari,tol))
  }
}
#What_t.z()


