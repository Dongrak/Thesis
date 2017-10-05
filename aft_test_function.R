#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------
W_t=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;
  
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time)
  p=length(b)
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]

  w_i=1*(Covari<=Covari[order(Covari)][100])
  
  N_i_s_t.beta=list(NA)
  for(j in 1:n){
    N_i_s_t.beta[[j]]=(e_i_beta>=e_i_beta[j])*Delta[j]
  }
  #N_i_s_t.beta
  
  Y_i_s_t.beta=list(NA)
  for(j in 1:n){
    Y_i_s_t.beta[[j]]=(e_i_beta<=e_i_beta[j])*1
  }
  #Y_i_s_t.beta
  
  N_d_s_t.beta=Reduce('+',N_i_s_t.beta)
  #N_d_s_t.beta
  
  S_0_s_t.beta=Reduce('+',Y_i_s_t.beta)
  #S_0_s_t.beta
  
  S_1_s_t.beta=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_s_t.beta,
                                 as.list(data.frame(t(Covari))),SIMPLIFY=FALSE))
  #S_1_s_t.beta
  
  E_s_t.beta=S_1_s_t.beta/S_0_s_t.beta
  #E_s_t.beta
  
  J_t.beta=(S_0_s_t.beta>0)*1
  #J_t.beta
  
  dN_d_s_t.beta=diff(c(0,N_d_s_t.beta))
  #dN_d_s_t.beta
  
  Ahat_0_t.beta=cumsum((J_t.beta/S_0_s_t.beta)*dN_d_s_t.beta)
  #Ahat_0_t.beta
  
  dAhat_0_t.beta=diff(c(0,Ahat_0_t.beta))
  #dAhat_0_t.beta
  
  Mhat_i_s_t.beta=mapply("-", N_i_s_t.beta,lapply(lapply(
    Y_i_s_t.beta,"*",dAhat_0_t.beta),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_s_t.beta
  
  W_t_test=(Reduce('+',mapply('*',Mhat_i_s_t.beta,w_i,SIMPLIFY = FALSE)))/sqrt(n)
  #W_t_test
  
  return(W_t_test)

  plot(W_t_test)
}
#W_t()
ee=W_t(beta_hat_wb,X_wb,D_wb,Z_wb)

#-------------------------------------------------------------
#-------------------------SAMPLE PATH-------------------------
#-------------------------------------------------------------
What_t=function(b,std,Time,Delta,Covari,weight,test,tol){
  #b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;weight=given_weight;test=given_test;tol=given_tol;
  #b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;weight=given_weight;test=given_test;tol=given_tol;
  
  n=length(Time)
  
  e_i_beta=log(Time)+Covari*b
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=Covari[order_resid]
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  if (weight=="a"){w_i=Covari*(Covari<=median(Covari))}
  if (weight=="b"){w_i=Covari}
  if (weight=="c"){w_i=1*(Covari<=median(Covari))}
  if (weight=="d"){w_i=1}
  
  N_i_s_t.beta=list(NA)
  for(j in 1:n){
    N_i_s_t.beta[[j]]=(e_i_beta>=e_i_beta[j])*Delta[j]
  }
  #N_i_s_t.beta
  
  Y_i_s_t.beta=list(NA)
  for(j in 1:n){
    Y_i_s_t.beta[[j]]=(e_i_beta<=e_i_beta[j])*1
  }
  #Y_i_s_t.beta
  
  N_d_s_t.beta=Reduce('+',N_i_s_t.beta)
  #N_d_s_t.beta
  
  S_0_s_t.beta=Reduce('+',Y_i_s_t.beta)
  #S_0_s_t.beta
  
  S_1_s_t.beta=Reduce('+',mapply(
    "*", Y_i_s_t.beta, Covari, SIMPLIFY = FALSE))
  #S_1_s_t.beta
  
  E_s_t.beta=S_1_s_t.beta/S_0_s_t.beta
  #E_s_t.beta
  
  J_t.beta=(S_0_s_t.beta>0)*1
  #J_t.beta
  
  dN_d_s_t.beta=diff(c(0,N_d_s_t.beta))
  #dN_d_s_t.beta
  
  Ahat_0_t.beta=cumsum((J_t.beta/S_0_s_t.beta)*dN_d_s_t.beta)
  #Ahat_0_t.beta
  
  dAhat_0_t.beta=diff(c(0,Ahat_0_t.beta))
  #dAhat_0_t.beta 
  
  Mhat_i_s_t.beta=mapply("-",N_i_s_t.beta,lapply(lapply(
    Y_i_s_t.beta,"*",dAhat_0_t.beta),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_s_t.beta
  
  dN_i_s_t.beta=lapply(N_i_s_t.beta,function(x){diff(c(0,x))})
  #dN_i_s_t.beta
  
  Q_t=S_0_s_t.beta/n # Gehan's weight
  #Q_t Q_t=1 #
  
  U_t.beta=Reduce('+',lapply(mapply("*",dN_i_s_t.beta,lapply(
    lapply(Covari,'-',E_s_t.beta),"*",Q_t), SIMPLIFY = FALSE),cumsum))
  #U_t.beta

  S_w_s_t.beta=Reduce('+',mapply(
    "*", Y_i_s_t.beta,w_i, SIMPLIFY = FALSE))
  #S_w_s_t.beta
  
  E_w_s_t.beta=S_w_s_t.beta/S_0_s_t.beta
  #E_w_s_t.beta
  
  dMhat_i_s_t.beta=lapply(Mhat_i_s_t.beta,function(x){diff(c(0,x))})
  #dMhat_i_s_t.beta
  
  #-----------------------------------------------------------
  #----------------------kernel Smoothing---------------------
  #-----------------------------------------------------------
  
  #-----------------------------g0----------------------------
  Ghat_0_t=1-exp(-Ahat_0_t.beta)
  #Ghat_0_t
  
  dGhat_0_t=diff(c(0,Ghat_0_t))
  #dGhat_0_t
  
  ghat_0_t=(ksmooth(e_i_beta,dGhat_0_t,"normal",
                    bandwidth = 1.06*sd(dGhat_0_t)*n^(-0.2),x.points=e_i_beta)$y)
  #ghat_0_t
  
  fhat_Y_t=Reduce('+',lapply(w_i*Covari,'*',ghat_0_t*Time))/n
  #fhat_Y_t
  
  #-----------------------------f0----------------------------
  KM_e=cumprod(1-(Delta/S_0_s_t.beta))
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

  fhat_N_t=Reduce('+',lapply(Delta*w_i*Covari,'*',fhat_0_t*Time))/n
  #fhat_N_t
  
  #-----------------------------------------------------------
  #--------Find Beta_hat_star by using optimize function------
  #-----------------------------------------------------------
  U_beta=function(beta_U=b,Time_U=Time,Delta_U=Delta,Covari_U=Covari){
    #beta_U=b;Time_U=Time;Delta_U=Delta;Covari_U=Covari;
    
    e_i_beta_U=log(Time_U)+Covari_U*beta_U
    
    order_resid_U=order(e_i_beta_U)
    
    Time_U=Time_U[order_resid_U]
    Covari_U=Covari_U[order_resid_U]
    Delta_U=Delta_U[order_resid_U]
    e_i_beta_U=e_i_beta_U[order_resid_U]
    
    N_i_s_t.beta_U=list(NA)
    for(j in 1:n){
      N_i_s_t.beta_U[[j]]=(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]
    }
    #N_i_s_t.beta_U
    
    Y_i_s_t.beta_U=list(NA)
    for(j in 1:n){
      Y_i_s_t.beta_U[[j]]=(e_i_beta_U<=e_i_beta_U[j])*1
    }
    #Y_i_s_t.beta_U
    
    N_d_s_t.beta_U=Reduce('+',N_i_s_t.beta_U)
    #N_d_s_t.beta_U
    
    S_0_s_t.beta_U=Reduce('+',Y_i_s_t.beta_U)
    #S_0_s_t.beta_U
    
    S_1_s_t.beta_U=Reduce('+',mapply(
      "*",Y_i_s_t.beta_U, Covari_U, SIMPLIFY = FALSE))
    #S_1_s_t.beta_U
    
    E_s_t.beta_U=S_1_s_t.beta_U/S_0_s_t.beta_U
    #E_s_t.beta
    
    J_t.beta_U=(S_0_s_t.beta_U>0)*1
    #J_t.beta_U
    
    dN_d_s_t.beta_U=diff(c(0,N_d_s_t.beta_U))
    #dN_d_s_t.beta_U
    
    Ahat_0_t.beta_U=cumsum((J_t.beta_U/S_0_s_t.beta_U)*dN_d_s_t.beta_U)
    #Ahat_0_t.beta_U
    
    dAhat_0_t.beta_U=diff(c(0,Ahat_0_t.beta_U))
    #dAhat_0_t.beta _U
    
    Mhat_i_s_t.beta=mapply("-", N_i_s_t.beta_U,lapply(
      lapply(Y_i_s_t.beta_U,"*",dAhat_0_t.beta_U),cumsum), SIMPLIFY = FALSE)
    #Mhat_i_s_t.beta _U
    
    dN_i_s_t.beta_U=lapply(N_i_s_t.beta_U,function(x){diff(c(0,x))})
    #dN_i_s_t.beta_U

    Q_t_U=S_0_s_t.beta_U/n # Gehan's weight
    #Q_t_U Q_t_U=1 #
    
    U_t.beta_U=Reduce('+',lapply(mapply("*",dN_i_s_t.beta_U,lapply(
      lapply(Covari_U,'-',E_s_t.beta_U),"*",Q_t_U), SIMPLIFY = FALSE),cumsum))
    #U_t.beta_U
    
    U_inf.beta_U=U_t.beta_U[n]
    #U_inf.beta_U
   
    return(U_inf.beta_U)
  }
  #U_beta()
  
  tolerance=tol+1 #initial value
  
  while (tolerance>tol){
    
    G_i=rnorm(n)
    #G_i
    
    U_w_G_t.beta=Reduce('+',lapply(mapply("*",mapply(
      "*",dMhat_i_s_t.beta,lapply(lapply(w_i,'-',E_w_s_t.beta),"*",Q_t)
      , SIMPLIFY = FALSE),G_i,SIMPLIFY = FALSE),cumsum))
    #U_w_G_t.beta
    
    U_G_t.beta=Reduce('+',lapply(mapply("*",mapply(
      "*",dMhat_i_s_t.beta,lapply(lapply(Covari,'-',E_s_t.beta),"*",Q_t)
      , SIMPLIFY = FALSE),G_i,SIMPLIFY = FALSE),cumsum))
    #U_G_t.beta

    U_G_inf.beta=U_G_t.beta[n]
    #U_G_t.inf.beta
    
    beta_hat_s_list=optimize(function(beta){abs(U_beta(beta_U=beta)-U_G_inf.beta)},
                             c(b-5*std,b+5*std),
                             tol = 1e-16);beta_hat_s_list
    #beta_hat_s_list
    
    tolerance=beta_hat_s_list$objective;tolerance
    #tolerance
  }
  
  beta_hat_s=beta_hat_s_list$minimum;beta_hat_s
  #beta_hat_s
  
  e_i_beta_s=log(Time)+Covari*beta_hat_s
  
  order_resid_s=order(e_i_beta_s)
  
  Time_s=Time[order_resid_s]
  Covari_s=Covari[order_resid_s]
  Delta_s=Delta[order_resid_s]
  e_i_beta_s=e_i_beta_s[order_resid_s]
  
  N_i_s_t.beta_s=list(NA)
  for(j in 1:n){
    N_i_s_t.beta_s[[j]]=(e_i_beta_s>=e_i_beta_s[j])*Delta_s[j]
  }
  #N_i_s_t.beta_s
  
  Y_i_s_t.beta_s=list(NA)
  for(j in 1:n){
    Y_i_s_t.beta_s[[j]]=(e_i_beta_s<=e_i_beta_s[j])*1
  }
  #Y_i_s_t.beta_s
  
  N_d_s_t.beta_s=Reduce('+',N_i_s_t.beta_s)
  #N_d_s_s_t.beta_s
  
  S_0_s_t.beta_s=Reduce('+',Y_i_s_t.beta_s)
  #S_0_s_t.beta_s
  
  S_1_s_t.beta_s=Reduce('+',mapply("*", Y_i_s_t.beta_s, Covari_s,
                                   SIMPLIFY = FALSE))
  #S_1_s_t.beta_s
  
  E_s_t.beta_s=S_1_s_t.beta_s/S_0_s_t.beta_s
  #E_s_t.beta_s
  
  J_t.beta_s=(S_0_s_t.beta_s>0)*1
  #J_t.beta_s
  
  dN_d_s_t.beta_s=diff(c(0,N_d_s_t.beta_s))
  #dN_d_s_t.beta_s
  
  Ahat_0_t.beta_s=cumsum((J_t.beta_s/S_0_s_t.beta_s)*dN_d_s_t.beta_s)
  #Ahat_0_t.beta_s
  
  dAhat_0_t.beta_s=diff(c(0,Ahat_0_t.beta_s))
  #dAhat_0_t.beta_s

  order_Time=order(Time)
  #order_Time
  
  order_Time_s=order(Time_s)
  #order_Time_s
  
  #order_Covari=order(Covari)
  
  F.T.=(1/sqrt(n))*U_w_G_t.beta
  S.T.=sqrt(n)*(fhat_N_t+cumsum(fhat_Y_t*dAhat_0_t.beta))*(b-beta_hat_s)
  T.T.=(1/sqrt(n))*cumsum(S_w_s_t.beta*diff(c(0,Ahat_0_t.beta-Ahat_0_t.beta_s)))
  
  What_t=F.T.-S.T.-T.T.
  #What_t
  
  #return(G_i)
  #return(beta_hat_s_list)
  #return(AA)
  #return(BB)
  #return(CC)
  #return(What_t)
  #return(list(What_t=What_t,beta_hat_s=beta_hat_s))
  
  if (test=="omni"){return(What_t[order_Time])}
  if (test=="ftn.form"){return(What_t)}
  
}
#What_t()

sample_path=function(path,b,std,Time,Delta,Covari,weight,test,tol){
  #path=path;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;weight=given_weight;test=given_test;tol=given_tol;
  #path=path;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;weight=given_weight;test=given_test;tol=given_tol;
 
  #------------------------SAMPLE PATH------------------------
  dataset_What=matrix(What_t(b,std,Time,Delta,Covari,weight,test,tol))
  for(k in 2:path){
    dataset_What=cbind(dataset_What,What_t(b,std,Time,Delta,Covari,weight,test,tol))
    if(k%%100==0) {
      cat("Sample Path",k,"\n")
    }
  }
  colnames(dataset_What)=paste0(rep("path"),1:path)
  rownames(dataset_What)=c(1:n)
  # dataset_What

  #------------------------BOOTSTRAPPING----------------------
  std.boot=as.vector(apply(dataset_What,1,sd))
  # std.boot
  
  dataset_std.What=dataset_What/std.boot
  # dataset_std.What
  
  dataset_W=W_t(b,Time,Delta,Covari,weight,test)
  # dataset_W
  
  dataset_std.W=dataset_W/std.boot
  # dataset_std.W
  
  #-----------------------MAXIMUM VALUE-----------------------
  max_path_What=as.vector(apply(abs(dataset_What),2,max))
  # max_path_What
  
  max_path_W=max(abs(dataset_W))
  # max_path_W
  
  max_path_std.What=as.vector(apply(abs(dataset_std.What),2,max))
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
  
  std_p_value=length(which((max_path_std.What>max_path_std.W)*1==1))/path
  # std_p_value
  
  result=list(dataset_What=dataset_What,dataset_std.What=dataset_std.What,
              dataset_W=dataset_W,dataset_std.W=dataset_std.W,
              std.boot=std.boot,p_value=p_value,std_p_value=std_p_value)
  # result
  
  return(result)
}
#sample_path

plotting=function(result,standardization,n.path){
  
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
#plotting

#-------------------------------------------------------------
#-------------------------OMNIBUS TEST------------------------
#-------------------------------------------------------------
#system.time(sample_path(path,beta_hat_wb,std_hat_wb,T_wb,D_wb,Z_wb,given_weight,given_test,given_tol))
#--------------------------CENSORING--------------------------
result_wb=sample_path(path,beta_hat_wb,std_hat_wb,X_wb,D_wb,Z_wb,given_weight,given_test,given_tol)
# Probability_wb
result_wb$p_value
result_wb$std_p_value
# PLOT : W_wb vs What_wb
Figure1_W_wb=plotting(result_wb,0,50);Figure1_W_wb
# PLOT : std.W_wb vs std.What_wb
Figure1_std.W_wb=plotting(result_wb,1,50);Figure1_std.W_wb

result_gg=sample_path(path,beta_hat_gg,std_hat_gg,X_gg,D_gg,Z_gg,given_weight,given_test,given_tol)
# Probability_gg
result_gg$p_value
result_gg$std_p_value
# PLOT : W_gg vs What_gg
Figure1_W_gg=plotting(result_gg,0,50);Figure1_W_gg
# PLOT : W_gg vs What_gg
Figure1_std.W_gg=plotting(result_gg,1,50);Figure1_std.W_gg

