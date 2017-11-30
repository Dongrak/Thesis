#-------------------------------------------------------------
#---------------------------AFTTEST---------------------------
#-------------------------------------------------------------
afttest_omni=function(path,b,std,Time,Delta,Covari,tol){
  # path=200;b=beta_hat_ln_aft;std=std_hat_ln_aft;Time=X_ln_aft;Delta=D_ln_aft;Covari=Z_ln_aft;tol=given_tol;
  # path=200;b=beta_hat_ln_cox;std=std_hat_ln_cox;Time=X_ln_cox;Delta=D_ln_cox;Covari=Z_ln_cox;tol=given_tol;
  # path=200;b=c(1.3,1.1);Covari=c(Z_ln_aft,Z_ln_aft^2-4*Z_ln_aft);
  
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
  order_Covari=apply(Covari,2,function(x){order(x)}) 
  pi_i_z=sapply(1:n,function(j){apply(sapply(1:p,function(i){c(rep(0,(which(
    order_Covari[,i]==j)-1)),rep(1,(n+1-which(order_Covari[,i]==j))))}),1,prod)},simplify=F)
  
  N_i_t=sapply(1:n,function(j){(e_i_beta>=e_i_beta[j])*Delta[j]},simplify=F)
  #N_i_t
  
  Y_i_t=sapply(1:n,function(j){(e_i_beta<=e_i_beta[j])*1},simplify=F)
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
  
  ghat_t.z=sapply(1:p,function(j){Reduce('+',lapply(mapply('*',pi_i_z,Covari[,j],
                                                           SIMPLIFY=FALSE),function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n},simplify=F)
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
  
  fhat_t.z=sapply(1:p,function(j){Reduce('+',lapply(mapply('*',pi_i_z,Delta*Covari[,j],
                                                           SIMPLIFY=FALSE),function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n},simplify=F)
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
    
    N_i_t_U=sapply(1:n,function(j){(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]},simplify=F)
    #N_i_t_U
    
    Y_i_t_U=sapply(1:n,function(j){(e_i_beta_U<=e_i_beta_U[j])*1},simplify=F)
    #Y_i_t_U
    
    S_0_t_U=Reduce('+',Y_i_t_U)
    #S_0_t_U
    
    S_1_t_U=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t_U,as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))
    #S_1_t_U
    
    dN_i_t_U=lapply(N_i_t_U,function(x){diff(c(0,x))})
    #dN_i_t_U
    
    U_inf_U=apply(S_0_t_U*Reduce('+',mapply(function(x,y){x%*%t(y)},dN_i_t_U,
                                            as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))-S_1_t_U*Reduce('+',dN_i_t_U),2,sum)/n
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  #------------------------SAMPLE PATH-----------------------
  
  app_path=list(NA)
  
  co=detectCores(logical=FALSE)-1 # number of core if logical is False else it means thread
  registerDoParallel(co)
  cl=makeCluster(co)
  app_path=foreach(k=1:path,.inorder=FALSE) %dopar% {
    
    # path_check=ceiling(path/2)
    # for (k inC 1:path){
    # if(k%%path_check==0) {
    #   cat("Sample Path",k,"\n")
    # }
    
    # tolerance=tol+1 #initial value
    
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
    
    N_i_t_s=sapply(1:n,function(j){(e_i_beta_s>=e_i_beta_s[j])*Delta_s[j]},simplify=F)
    #N_i_t_s
    
    Y_i_t_s=sapply(1:n,function(j){(e_i_beta_s<=e_i_beta_s[j])*1},simplify=F)
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
    
    app_path=F.T.-S.T.-T.T.
    # app_path[[k]]=F.T.-S.T.-T.T.
    #app_path
  }
  stopCluster(cl)
  closeAllConnections()
  
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
  
  # result=list(Time,Delta,Covari,e_i_beta,std.boot,
  #             app_path,app_std.path,
  #             obs_path,obs_std.path,
  #             p_value,std.p_value)
  # 
  # names(result)=c("Time","Delta","Covari","Resid","std.boot",
  #                 "app_path","app_std.path",
  #                 "obs_path","obs_std.path",
  #                 "p_value","std.p_value")
  
  result=list(p_value,std.p_value);names(result)=c("p_value","std.p_value");
  result
  
  rm(list=(ls()[ls()!="result"]));gc();
  return(result)
}
#afttest_omni()

afttest_form=function(path,b,std,Time,Delta,Covari,tol,form=1){
  # path=200;b=beta_hat_ln_aft;std=std_hat_ln_aft;Time=X_ln_aft;Delta=D_ln_aft;Covari=Z_ln_aft;tol=given_tol;
  # path=200;b=beta_hat_ln_cox;std=std_hat_ln_cox;Time=X_ln_cox;Delta=D_ln_cox;Covari=Z_ln_cox;tol=given_tol;
  # path=200;b=c(1.3,1.1);Covari=c(Z_ln_aft,Z_ln_aft^2-4*Z_ln_aft);
  
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
  order_Covari_form=order(Covari[,form])
  pi_i_z=sapply(1:n,function(j){c(rep(0,(which(order_Covari_form==j)-1)),rep(1,(n+1-which(order_Covari_form==j))))},simplify=F)
  
  N_i_t=sapply(1:n,function(j){(e_i_beta>=e_i_beta[j])*Delta[j]},simplify=F)
  #N_i_t
  
  Y_i_t=sapply(1:n,function(j){(e_i_beta<=e_i_beta[j])*1},simplify=F)
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
  
  ghat_t.z=sapply(1:p,function(j){Reduce('+',lapply(mapply('*',pi_i_z,Covari[,j],
                                                           SIMPLIFY=FALSE),function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n},simplify=F)
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
  
  fhat_t.z=sapply(1:p,function(j){Reduce('+',lapply(mapply('*',pi_i_z,Delta*Covari[,j],
                                                           SIMPLIFY=FALSE),function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n},simplify=F)
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
    
    N_i_t_U=sapply(1:n,function(j){(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]},simplify=F)
    #N_i_t_U
    
    Y_i_t_U=sapply(1:n,function(j){(e_i_beta_U<=e_i_beta_U[j])*1},simplify=F)
    #Y_i_t_U
    
    S_0_t_U=Reduce('+',Y_i_t_U)
    #S_0_t_U
    
    S_1_t_U=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t_U,as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))
    #S_1_t_U
    
    dN_i_t_U=lapply(N_i_t_U,function(x){diff(c(0,x))})
    #dN_i_t_U
    
    U_inf_U=apply(S_0_t_U*Reduce('+',mapply(function(x,y){x%*%t(y)},dN_i_t_U,
                                            as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))-S_1_t_U*Reduce('+',dN_i_t_U),2,sum)/n
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  #------------------------SAMPLE PATH-----------------------
  
  app_path=list(NA)
  
  co=detectCores(logical=FALSE)-1 # number of core if logical is False else it means thread
  registerDoParallel(co)
  cl=makeCluster(co)
  app_path=foreach(k=1:path,.inorder=FALSE) %dopar% {
    
    # path_check=ceiling(path/2)
    # for (k inC 1:path){
    # if(k%%path_check==0) {
    #   cat("Sample Path",k,"\n")
    # }
    
    # tolerance=tol+1 #initial value
    
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
    
    N_i_t_s=sapply(1:n,function(j){(e_i_beta_s>=e_i_beta_s[j])*Delta_s[j]},simplify=F)
    #N_i_t_s
    
    Y_i_t_s=sapply(1:n,function(j){(e_i_beta_s<=e_i_beta_s[j])*1},simplify=F)
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
    
    app_path=F.T.-S.T.-T.T.
    # app_path[[k]]=F.T.-S.T.-T.T.
    #app_path
  }
  stopCluster(cl)
  closeAllConnections()
  
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
  
  # result=list(Time,Delta,Covari,e_i_beta,std.boot,
  #             app_path,app_std.path,
  #             obs_path,obs_std.path,
  #             p_value,std.p_value)
  # 
  # names(result)=c("Time","Delta","Covari","Resid","std.boot",
  #                 "app_path","app_std.path",
  #                 "obs_path","obs_std.path",
  #                 "p_value","std.p_value")
  
  result=list(p_value,std.p_value);names(result)=c("p_value","std.p_value");
  # result
  
  rm(list=(ls()[ls()!="result"]));gc();
  return(result)
}
#afttest_form()

afttest_link=function(path,b,std,Time,Delta,Covari,tol){
  # path=200;b=beta_hat_ln_aft;std=std_hat_ln_aft;Time=X_ln_aft;Delta=D_ln_aft;Covari=Z_ln_aft;tol=given_tol;
  # path=200;b=beta_hat_ln_cox;std=std_hat_ln_cox;Time=X_ln_cox;Delta=D_ln_cox;Covari=Z_ln_cox;tol=given_tol;
  # path=200;b=c(1.3,1.1);Covari=c(Z_ln_aft,Z_ln_aft^2-4*Z_ln_aft);
  
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
  order_Covari=apply(Covari,2,function(x){order(x)}) 
  pi_i_z=sapply(1:n,function(j){apply(sapply(1:p,function(i){c(rep(0,(which(
    order_Covari[,i]==j)-1)),rep(1,(n+1-which(order_Covari[,i]==j))))}),1,prod)},simplify=F)
  
  N_i_t=sapply(1:n,function(j){(e_i_beta>=e_i_beta[j])*Delta[j]},simplify=F)
  #N_i_t
  
  Y_i_t=sapply(1:n,function(j){(e_i_beta<=e_i_beta[j])*1},simplify=F)
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
  
  ghat_t.z=sapply(1:p,function(j){Reduce('+',lapply(mapply('*',pi_i_z,Covari[,j],
                                                           SIMPLIFY=FALSE),function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n},simplify=F)
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
  
  fhat_t.z=sapply(1:p,function(j){Reduce('+',lapply(mapply('*',pi_i_z,Delta*Covari[,j],
                                                           SIMPLIFY=FALSE),function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n},simplify=F)
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
    
    N_i_t_U=sapply(1:n,function(j){(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]},simplify=F)
    #N_i_t_U
    
    Y_i_t_U=sapply(1:n,function(j){(e_i_beta_U<=e_i_beta_U[j])*1},simplify=F)
    #Y_i_t_U
    
    dN_i_t_U=lapply(N_i_t_U,function(x){diff(c(0,x))})
    #dN_i_t_U
    
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
  
  #------------------------SAMPLE PATH-----------------------
  
  app_path=list(NA)
  
  co=detectCores(logical=FALSE)-1 # number of core if logical is False else it means thread
  registerDoParallel(co)
  cl=makeCluster(co)
  app_path=foreach(k=1:path,.inorder=FALSE) %dopar% {
    
    # path_check=ceiling(path/2)
    # for (k inC 1:path){
    # if(k%%path_check==0) {
    #   cat("Sample Path",k,"\n")
    # }
    
    # tolerance=tol+1 #initial value
    
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
    
    N_i_t_s=sapply(1:n,function(j){(e_i_beta_s>=e_i_beta_s[j])*Delta_s[j]},simplify=F)
    #N_i_t_s
    
    Y_i_t_s=sapply(1:n,function(j){(e_i_beta_s<=e_i_beta_s[j])*1},simplify=F)
    #Y_i_t_s
    
    N_d_t_s=Reduce('+',N_i_t_s)
    
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
    
    app_path=F.T.-S.T.-T.T.
    # app_path[[k]]=F.T.-S.T.-T.T.
    #app_path
  }
  stopCluster(cl)
  closeAllConnections()
  
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
  
  # result=list(Time,Delta,Covari,e_i_beta,std.boot,
  #             app_path,app_std.path,
  #             obs_path,obs_std.path,
  #             p_value,std.p_value)
  # 
  # names(result)=c("Time","Delta","Covari","Resid","std.boot",
  #                 "app_path","app_std.path",
  #                 "obs_path","obs_std.path",
  #                 "p_value","std.p_value")
  
  result=list(p_value,std.p_value);names(result)=c("p_value","std.p_value");
  # result
  
  rm(list=(ls()[ls()!="result"]));gc();
  return(result)
}
#afttest_link()

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
    return(c(list(testtype=testtype),afttest_omni(path,b,std,Time,Delta,Covari,tol)))
  }
  if(testtype=="form"){
    return(c(list(testtype=testtype),afttest_form(path,b,std,Time,Delta,Covari,tol,form)))
  }
  if(testtype=="link"){
    return(c(list(testtype=testtype),afttest_link(path,b,std,Time,Delta,Covari,tol)))
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
  result_log_Time=log(result_Time)
  n=length(result_Time)
  if (xaxix=="rank"){xaxix=(1:n)[order(result_Time)]}
  else if (xaxix=="real"){xaxix=result_Time}
  else (xaxix=result_Time)
  
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
