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

simulation=100
n=250
path=150
alpha=0.05

given_tol=0.1

#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------
W_omni=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  #b=c(1.3,1.1);Covari=c(Z_wb,Z_wb^2-Z_wb);
  
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

#-------------------------------------------------------------
#-------------------------REALIZATION-------------------------
#-------------------------------------------------------------
What_omni=function(b,std,Time,Delta,Covari,tol){
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

#-------------------------------------------------------------
#-------------------------SAMPLE PATH-------------------------
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


#-------------------------------------------------------------
#--------------------------SIMULATION-------------------------
#-------------------------------------------------------------
simulation_omni=function(simulation,n,path,alpha,tol){
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
    Z=matrix(rnorm(n,3,1),nrow=n)
    
    #-------------------LOG NORMAL DISTRIBUTION-------------------
    T_ln_aft=as.vector(exp(-beta_0*Z)*qlnorm(runif(n),5,1))
    C_ln_aft=as.vector(exp(-beta_0*Z)*qlnorm(runif(n),6.5,1))
    X_ln_aft=C_ln_aft*(T_ln_aft>C_ln_aft)+T_ln_aft*(T_ln_aft<=C_ln_aft)
    D_ln_aft=0*(T_ln_aft>C_ln_aft)+1*(T_ln_aft<=C_ln_aft)
    Z_ln_aft=Z
    
    T_ln_cox=as.vector(qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5,1,lower.tail=FALSE))
    C_ln_cox=as.vector(qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5.7,1,lower.tail=FALSE))
    X_ln_cox=C_ln_cox*(T_ln_cox>C_ln_cox)+T_ln_cox*(T_ln_cox<=C_ln_cox)
    D_ln_cox=0*(T_ln_cox>C_ln_cox)+1*(T_ln_cox<=C_ln_cox)
    Z_ln_cox=Z
    
    #------------Estimate Beta_hat_ln_aft_f by using Aftgee-----------
    aftsrr_beta_ln_aft=aftsrr(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,method="nonsm")
    beta_hat_ln_aft=-as.vector(aftsrr_beta_ln_aft$beta);beta_hat_ln_aft
    std_hat_ln_aft=diag(aftsrr_beta_ln_aft$covmat$ISMB);std_hat_ln_aft
    
    aftsrr_beta_ln_cox=aftsrr(Surv(X_ln_cox,D_ln_cox)~Z_ln_cox,method="nonsm")
    beta_hat_ln_cox=-as.vector(aftsrr_beta_ln_cox$beta);beta_hat_ln_cox
    std_hat_ln_cox=diag(aftsrr_beta_ln_cox$covmat$ISMB);std_hat_ln_cox
    
    # result_ln_aft
    result_ln_aft=sample_path_omni(path,beta_hat_ln_aft,std_hat_ln_aft,
                                   X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
    
    # result_ln_cox
    result_ln_cox=sample_path_omni(path,beta_hat_ln_cox,std_hat_ln_cox,
                                   X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
    
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
    
    p_mean=rbind(c(result_ln_aft$p_value,result_ln_aft$std.p_value),
                 c(result_ln_cox$p_value,result_ln_cox$std.p_value))
    colnames(p_mean)=c("W","std.W")
    rownames(p_mean)=c("p_ln_aft_mean","p_ln_cox_mean")
    #p_mean
    
    p_alpha=(p_mean>=alpha)*1
    colnames(p_alpha)=c("W","std.W")
    rownames(p_alpha)=c("p_ln_aft_alpha","p_ln_cox_alpha")
    #p_alpha
    
    p_value=list(p_mean,p_alpha)
    #p_value
    
    #result[[k]]=list(result_ln_aft,result_ln_cox,p_value)
    result[[k]]=list(p_value)
  }
  return(result)
}
#simulation_omni

prob.table_omni=function(simul_result){
  simul=length(simul_result)
  
  p_mean_set=list(NA)
  p_alpha_set=list(NA)
  
  for(k in 1:simul){
    # p_mean_set[[k]]=simul_result[[k]][[3]][[1]]
    # p_alpha_set[[k]]=simul_result[[k]][[3]][[2]]
    p_mean_set[[k]]=simul_result[[k]][[1]][[1]]
    p_alpha_set[[k]]=simul_result[[k]][[1]][[2]]
  }
  
  p_mean=Reduce("+",p_mean_set)/simul
  p_alpha=Reduce("+",p_alpha_set)/simul
  
  return(list(p_mean,p_alpha))
}
#prob.table_omni

date()
simulation_result_omni1=simulation_omni(simulation,n,path,alpha,given_tol)
prob.table_omni(simulation_result_omni1)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_omni1_n500p150sim100")
date()
simulation_result_omni2=simulation_omni(simulation,n,path,alpha,given_tol)
prob.table_omni(simulation_result_omni2)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_omni2_n500p150sim100")
date()
simulation_result_omni3=simulation_omni(simulation,n,path,alpha,given_tol)
prob.table_omni(simulation_result_omni3)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_omni3_n500p150sim100")
date()
simulation_result_omni4=simulation_omni(simulation,n,path,alpha,given_tol)
prob.table_omni(simulation_result_omni4)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_omni4_n500p150sim100")
date()
simulation_result_omni5=simulation_omni(simulation,n,path,alpha,given_tol)
prob.table_omni(simulation_result_omni5)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_omni5_n500p150sim100")
date()
simulation_result_omni=c(simulation_result_omni1,simulation_result_omni2,
                         simulation_result_omni3,simulation_result_omni4,
                         simulation_result_omni5)
prob.table_omni(simulation_result_omni)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_omni_n500p150sim500")
date()