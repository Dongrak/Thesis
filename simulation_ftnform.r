#-------------------------------------------------------------
#---------------------------SETTING---------------------------
#-------------------------------------------------------------
#rm(list=ls());gc();

memory.limit(16*2^20)

options(max.print=999999)
options(error=NULL)

#install.packages("ggplot2")
#install.packages("survival")
#install.packages("aftgee")
#install.packages("ENmisc")
#install.packages("plotly")

library(ggplot2)
library(survival)
library(aftgee)
library(ENmisc)
library(plotly)

simulation=200
n=200
path=200
alpha=0.05

given_tol=0.1

#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------
W_ftnform=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-5*Z_wb);
  
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
#W_ftnform()

#-------------------------------------------------------------
#-------------------------REALIZATION-------------------------
#-------------------------------------------------------------
What_ftnform=function(b,std,Time,Delta,Covari,tol){
  #b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parametersa
  
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
#What_ftnform()

#-------------------------------------------------------------
#-------------------------SAMPLE PATH-------------------------
#-------------------------------------------------------------
sample_path_ftnform=function(path,b,std,Time,Delta,Covari,tol){
  #path=path;b=beta_hat_gg;std=std_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg;tol=given_tol;
  #path=path;b=beta_hat_wb;std=std_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb;tol=given_tol;
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
  #------------------------SAMPLE PATH------------------------
  dataset_What=list(NA)
  for(k in 1:path){
    dataset_What[[k]]=What_ftnform(b,std,Time,Delta,Covari,tol)$sim_stat
    if(k%%100==0) {
      cat("Sample Path",k,"\n")
    }
  }
  
  #------------------------BOOTSTRAPPING----------------------
  std.boot=matrix(apply(mapply(function(x){as.vector(x)},dataset_What),1,sd),nrow=n)
  # std.boot
  
  dataset_std.What=lapply(dataset_What,function(x){x/std.boot})
  # dataset_std.What
  
  dataset_W=W_ftnform(b,Time,Delta,Covari)$obs_stat
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
#sample_path_ftnform

#-------------------------------------------------------------
#-----------------------FUNCTIONAL FORM-----------------------
#-------------------------------------------------------------
simulation_ftnform=function(simulation,n,path,alpha,tol){
  #simulation=simulation;n=n;path=path;alpha=alpha;tol=given_tol;
  
  result=list(NA)
  
  for(k in 1:simulation){
    if(k%%1==0) {
      cat("simulation",k,"\n")
    }
    # -------------------------------------------------------------
    # ------------------------DATA GENERATE------------------------
    # -------------------------------------------------------------
    # n=200
    beta_0=1
    gamma_0=0.1
    Z=matrix(rnorm(n,3,1),nrow=n)
    
    #-------------------LOG NORMAL DISTRIBUTION-------------------
    T_ln_aft_f=as.vector(exp(-beta_0*Z-gamma_0*(Z^2))*qlnorm(runif(n),5,1))
    C_ln_aft_f=as.vector(exp(-beta_0*Z-gamma_0*(Z^2))*qlnorm(runif(n),6.5,1))
    X_ln_aft_f=C_ln_aft_f*(T_ln_aft_f>C_ln_aft_f)+T_ln_aft_f*(T_ln_aft_f<=C_ln_aft_f)
    D_ln_aft_f=0*(T_ln_aft_f>C_ln_aft_f)+1*(T_ln_aft_f<=C_ln_aft_f)
    Z_ln_aft_f=Z
    
    #------------Estimate Beta_hat_ln_aft_f by using Aftgee-----------
    aftsrr_beta_ln_aft_f=aftsrr(Surv(X_ln_aft_f,D_ln_aft_f)~Z_ln_aft_f,method="nonsm")
    beta_hat_ln_aft_f=-as.vector(aftsrr_beta_ln_aft_f$beta);beta_hat_ln_aft_f
    std_hat_ln_aft_f=diag(aftsrr_beta_ln_aft_f$covmat$ISMB);std_hat_ln_aft_f
    
    # result_ln_aft_f
    result_ln_aft_f=sample_path_ftnform(path,beta_hat_ln_aft_f,std_hat_ln_aft_f,
                                        X_ln_aft_f,D_ln_aft_f,Z_ln_aft_f,given_tol)
    
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
    
    p_mean=rbind(c(result_ln_aft_f$p_value,result_ln_aft_f$std.p_value))
    colnames(p_mean)=c("W","std.W")
    rownames(p_mean)=c("p_ln_aft_mean")
    #p_mean
    
    p_alpha=(p_mean>=alpha)*1
    colnames(p_alpha)=c("W","std.W")
    rownames(p_alpha)=c("p_ln_aft_alpha")
    #p_alpha
    
    p_value=list(p_mean,p_alpha)
    #p_value
    
    result[[k]]=list(result_ln_aft_f,p_value)
  }
  return(result)
}
#simulation_ftnform

prob.table_ftnform=function(simul_result){
  simul=length(simul_result)
  
  p_mean_set=list(NA)
  p_alpha_set=list(NA)
  
  for(k in 1:simul){
    p_mean_set[[k]]=simul_result[[k]][[2]][[1]]
    p_alpha_set[[k]]=simul_result[[k]][[2]][[2]]
  }
  
  p_mean=Reduce("+",p_mean_set)/simul
  p_alpha=Reduce("+",p_alpha_set)/simul
  
  return(list(p_mean,p_alpha))
}
#prob.table_ftnform

date()
simulation_result_ftnform1=simulation_ftnform(simulation,n,path,alpha,given_tol)
prob.table_ftnform(simulation_result_ftnform1)
date()
simulation_result_ftnform2=simulation_ftnform(simulation,n,path,alpha,given_tol)
prob.table_ftnform(simulation_result_ftnform2)
date()
simulation_result_ftnform3=simulation_ftnform(simulation,n,path,alpha,given_tol)
prob.table_ftnform(simulation_result_ftnform3)
date()
simulation_result_ftnform4=simulation_ftnform(simulation,n,path,alpha,given_tol)
prob.table_ftnform(simulation_result_ftnform4)
date()
simulation_result_ftnform5=simulation_ftnform(simulation,n,path,alpha,given_tol)
prob.table_ftnform(simulation_result_ftnform5)
date()
simulation_result_ftnform=c(simulation_result_ftnform1,simulation_result_ftnform2,
                            simulation_result_ftnform3,simulation_result_ftnform4,
                            simulation_result_ftnform5)
prob.table_ftnform(simulation_result_ftnform)
date()


