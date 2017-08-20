#-------------------------------------------------------------
#---------------------------SETTING---------------------------
#-------------------------------------------------------------
rm(list=ls());gc();

memory.limit(16*2^20)

options(max.print=999999)
options(error=NULL)

#install.packages("ggplot2")
#install.packages("survival")
#install.packages("aftgee")
#install.packages("BB")

library(ggplot2)
library(survival)
library(aftgee)
library(BB)

#-------------------------------------------------------------
#-----------------------COX DATA ANALYSIS---------------------
#-------------------------------------------------------------
W_t=function(b,Time,Delta,Covari,weight,test){
  #b=beta_hat_cox;weight=w_i;Time=T_cox;Delta=D_cox;Covari=Z_cox;
  
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
  
  Mhat_i_s_t.beta=mapply("-", N_i_s_t.beta,lapply(lapply(
    Y_i_s_t.beta,"*",dAhat_0_t.beta),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_s_t.beta
  
  W_t=n^(-0.5)*(Reduce('+',mapply('*',Mhat_i_s_t.beta,w_i,SIMPLIFY = FALSE)))
  #W_t
  
  if (test=="omni"){return(W_t[order(Time)])}
  if (test=="ftn.form"){return(W_t[order(Covari)])}
  
}
#W_t()

#-------------------------------------------------------------
#-----------------------AFT DATA ANALYSIS---------------------
#-------------------------------------------------------------
What_t=function(b,std,Time,Delta,Covari,weight,test,tol){
  #b=beta_hat_aft;weight=given_weight;Time=T_aft;Delta=D_aft;Covari=Z_aft;test=given_test;tol=given_tol;

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
  
  Q_t=1 #Q_t=S_0_s_t.beta/n # Gehan's weight
  #Q_t
  
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
  KM_e=cumprod(1-Delta/S_0_s_t.beta)
  #KM_e
  
  Fhat_0_e=1-KM_e
  #Fhat_0_e
  
  dFhat_0_e=diff(c(0,Fhat_0_e))
  #dFhat_0_e
  
  Condi.Ehat=(cumsum(e_i_beta*dFhat_0_e))/(1-Fhat_0_e)
  #Condi.Ehat
  
  rhat_i=Delta*e_i_beta+(1-Delta)*Condi.Ehat
  #rhat_i
  
  fhat_0_t=ksmooth(e_i_beta,dFhat_0_e,"normal",
                   bandwidth = 1.06*sd(dFhat_0_e)*n^(-0.2),x.points=e_i_beta)$y
  #fhat_0_t

  fhat_N_t=Reduce('+',lapply(Delta*w_i*Covari,'*',fhat_0_t*Time))/n
  #fhat_N_t
  
  #-----------------------------------------------------------
  #--------Find Beta_hat_star by using optimize function------
  #-----------------------------------------------------------
  U_beta=function(beta_U=b,Time_U=Time,Delta_U=Delta,Covari_U=Covari){
    #beta_U=b;Time_U=Time;Delta_U=Delta;Covari_U=Covari;Q_t_U=Q_t;
    
    e_i_beta_U=log(Time_U)+Covari_U*beta_U
    Time_U=Time_U[order(e_i_beta_U)]
    Covari_U=Covari_U[order(e_i_beta_U)]
    Delta_U=Delta_U[order(e_i_beta_U)]
    e_i_beta_U=e_i_beta_U[order(e_i_beta_U)]
    
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
    
    S_1_s_t.beta_U=Reduce('+',mapply("*",
                                     Y_i_s_t.beta_U, Covari_U, SIMPLIFY = FALSE))
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
    #Q_t_U
    
    U_t.beta_U=Reduce('+',lapply(mapply("*",dN_i_s_t.beta_U,lapply(
      lapply(Covari_U,'-',E_s_t.beta_U),"*",Q_t_U), SIMPLIFY = FALSE),cumsum))
    #U_t.beta_U
    
    U_t.beta_U_order=U_t.beta_U[order(Time_U)]
    #U_t.beta_U_order
    
    U_inf.beta_U=U_t.beta_U_order[n]
    #U_inf.beta_U
   
    return(U_inf.beta_U)
  }
  #U_beta()
  
  tolerance=tol+1 #initial value
  
  while (tolerance>tol){
    
    G_i=rnorm(n)
    #G_i
    
    U_w_G_t.beta=Reduce('+',lapply(lapply(mapply(
      "*",dMhat_i_s_t.beta,lapply(lapply(w_i,'-',E_w_s_t.beta),"*",Q_t)
      , SIMPLIFY = FALSE),"*",G_i),cumsum))
    #U_w_G_t.beta
    
    U_G_t.beta=Reduce('+',lapply(lapply(mapply(
      "*",dMhat_i_s_t.beta,lapply(lapply(Covari,'-',E_s_t.beta),"*",Q_t)
      , SIMPLIFY = FALSE),"*",G_i),cumsum))
    #U_G_t.beta
    
    U_G_t.beta_order=U_G_t.beta[order(Time)]
    #U_G_t.beta_order
    
    U_G_t.inf.beta=U_G_t.beta_order[n]
    #U_G_t.inf.beta
    
    beta_hat_s_list=optimize(function(beta){abs(U_beta(beta)-U_G_t.inf.beta)},
                             c(b-std,b+std),
                             tol = 1e-16)
    #beta_hat_s_list
    
    tolerance=beta_hat_s_list$objective
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
  
  AA=(1/sqrt(n))*U_w_G_t.beta;AA
  BB=sqrt(n)*(fhat_N_t+cumsum(fhat_Y_t*dAhat_0_t.beta))*(b-beta_hat_s);BB
  CC=(1/sqrt(n))*cumsum(S_w_s_t.beta*diff(c(0,Ahat_0_t.beta-Ahat_0_t.beta_s)));CC
  
  What_t=AA-BB-CC
  #What_t
  
  #return(G_i)
  #return(beta_hat_s_list)
  #return(AA)
  #return(BB)
  #return(CC)
  #return(What_t)
  #return(list(What_t=What_t,beta_hat_s=beta_hat_s))
  
  if (test=="omni"){return(What_t[order(Time)])}
  if (test=="ftn.form"){return(What_t[order(Covari)])}
  
}
#What_t()

sample_path=function(path,b,std,Time,Delta,Covari,weight,test,tol){
  # path=path;b=beta_hat_aft;weight=given_weight;Time=T_aft;Delta=D_aft;Covari=Z_aft;test=given_test;tol=given_tol;
  
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
  
  #------------------------SAMPLE PATH------------------------
  #------------------------SAMPLE PATH------------------------
  
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
  
  std.max_path_What=as.vector(apply(abs(dataset_std.What),2,max))
  # std.max_path_What
  
  std.max_path_W=max(abs(dataset_std.W))
  # std.max_path_W
  
  #-----------------------------------------------------------
  #  p  : the ratio of (What>=W)*1
  # H_0 : the data follow the assumption of the aft model.
  #
  # if p>0.95, cannot reject the null hypothesis. i.e. accept it 
  # if p<0.95, reject the null hypothesis.
  #
  # absolute/maximum 기준으로 What이 큰것의 비율(p)이 
  # 0.96이면 당연히 accetp
  # 0.04이면 당연히 reject
  # 0.45이면 reject <= W가 55번이나 튀어나간거다!
  #-----------------------------------------------------------
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_path_What>max_path_W)*1==1))/path
  # p_value
  
  std.p_value=length(which((std.max_path_What>std.max_path_W)*1==1))/path
  # std.p_value
  
  result=list(dataset_What=dataset_What,dataset_std.What=dataset_std.What,
              dataset_W=dataset_W,dataset_std.W=dataset_std.W,
              std.boot=std.boot,p_value=p_value,std.p_value=std.p_value)
  # result
  
  return(result)
}
#sample_path

plotting=function(result,standardization){
  
  if (standardization==0) {
    dataset_What=data.frame()
    for (i in 1:path){
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
      geom_line(data=dataset_What,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
      geom_line(data=dataset_W,aes(x=t_i,y=W),colour="tomato")
    #Figure1_W
    
    return(Figure1_W)
  }
  if (standardization==1) {
    dataset_std.What=data.frame()
    for (i in 1:path){
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
      geom_line(data=dataset_std.What,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
      geom_line(data=dataset_std.W,aes(x=t_i,y=std.W),colour="tomato")
    return(Figure1_std.W)
  }
}
#plotting

#-------------------------------------------------------------
#------------------------DATA GENERATE------------------------
#-------------------------------------------------------------
path=200

n=200
id=c(1:n) # identification
beta_0=1 # beta_0
gamma_0=0.1 # gamma_0
Z=rnorm(n,3,1) # covariate

#-----------------------AFT DATA GENERATE---------------------
T_s_aft=exp(-beta_0*Z)*qlnorm(runif(n),5,1) # lognormal baseline hazard
C_aft=exp(-beta_0*Z)*qlnorm(runif(n),6.5,1) # censoring time for aft model
T_aft=C_aft*(T_s_aft>C_aft)+T_s_aft*(T_s_aft<=C_aft) # observed time for aft model
D_aft=0*(T_s_aft>C_aft)+1*(T_s_aft<=C_aft) # delta
D_aft=D_aft[order(T_aft)]
Z_aft=Z[order(T_aft)]
T_aft=T_aft[order(T_aft)]

#-------------AFT DATA GENERATE (FUNCTIONAL FORM)-------------
T_s_aft_f=exp(-beta_0*Z-gamma_0*Z^2)*qlnorm(runif(n),5,1) # lognormal baseline hazard
C_aft_f=exp(-beta_0*Z-gamma_0*Z^2)*qlnorm(runif(n),6.5,1) # censoring time for aft model
T_aft_f=C_aft_f*(T_s_aft_f>C_aft_f)+T_s_aft_f*(T_s_aft_f<=C_aft_f) # observed time for aft model
D_aft_f=0*(T_s_aft_f>C_aft_f)+1*(T_s_aft_f<=C_aft_f) # delta
D_aft_f=D_aft_f[order(T_aft_f)]
Z_aft_f=Z[order(T_aft_f)]
T_aft_f=T_aft_f[order(T_aft_f)]

#-----------------------COX DATA GENERATE---------------------
T_s_cox=qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5,1,lower.tail = FALSE) # lognormal baseline hazard
C_cox=qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5.7,1,lower.tail = FALSE) # censoring time for cox model
T_cox=C_cox*(T_s_cox > C_cox)+T_s_cox*(T_s_cox<=C_cox) # observed time for cox model
D_cox=0*(T_s_cox > C_cox)+1*(T_s_cox<=C_cox) # delta
D_cox=D_cox[order(T_cox)]
Z_cox=Z[order(T_cox)]
T_cox=T_cox[order(T_cox)]

#-------------------------------------------------------------
#------------Estimate Beta_hat_aft by using Aftgee------------
#-------------------------------------------------------------
aftsrr_beta_aft=aftsrr(Surv(T_aft,D_aft)~Z_aft,method="nonsm")
beta_hat_aft=-unlist(summary(aftsrr_beta_aft))$coefficients1;beta_hat_aft
std_hat_aft=unlist(summary(aftsrr_beta_aft))$coefficients2;std_hat_aft

#-------------------------------------------------------------
#-----------Estimate Beta_hat_aft_f by using Aftgee-----------
#-------------------------------------------------------------
aftsrr_beta_aft_f=aftsrr(Surv(T_aft_f,D_aft_f)~Z_aft_f,method="nonsm")
beta_hat_aft_f=-unlist(summary(aftsrr_beta_aft_f))$coefficients1;beta_hat_aft_f
std_hat_aft_f=unlist(summary(aftsrr_beta_aft_f))$coefficients2;std_hat_aft_f

#-------------------------------------------------------------
#------------Estimate Beta_hat_cox by using Aftgee------------
#-------------------------------------------------------------
aftsrr_beta_cox=aftsrr(Surv(T_cox,D_cox)~Z_cox,method="nonsm")
beta_hat_cox=-unlist(summary(aftsrr_beta_cox))$coefficients1;beta_hat_cox
std_hat_cox=unlist(summary(aftsrr_beta_cox))$coefficients2;std_hat_cox

#-------------------------------------------------------------
#------------------------WEIGHT&TOLERANCE---------------------
#-------------------------------------------------------------
given_tol=0.1

given_weight="c"
#(weight=="a"){w_i=Covari*(Covari<=median(Covari))}
#(weight=="b"){w_i=Covari}  
#(weight=="c"){w_i=1*(Covari<=median(Covari))}
#(weight=="d"){w_i=1}

given_test="omni"
# test="link.ftn"
# test="ftn.form"
# test="....???"

#-------------------------------------------------------------
#-------------------------OMNIBUS TEST------------------------
#-------------------------------------------------------------
#system.time(sample_path(path,beta_hat_aft,std_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,given_tol))
#--------------------------CENSORING--------------------------
result_aft=sample_path(path,beta_hat_aft,std_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,given_tol)
result_cox=sample_path(path,beta_hat_cox,std_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,given_tol)

result_aft$p_value
result_aft$std.p_value

result_cox$p_value
result_cox$std.p_value

# PLOT : W_aft vs What_aft
Figure1_W_aft=plotting(result_aft,0);Figure1_W_aft

# PLOT : std.W_aft vs std.What_aft
Figure1_std.W_aft=plotting(result_aft,1);Figure1_std.W_aft

# PLOT : W_cox vs What_cox
Figure1_W_cox=plotting(result_cox,0);Figure1_W_cox

# PLOT : W_cox vs What_cox
Figure1_std.W_cox=plotting(result_cox,1);Figure1_std.W_cox



############################고쳐야함 







#-------------------------NONCENSORING------------------------
dataset_What_aft_NC=sample_path_What(path,beta_hat_aft,T_s_aft,rep(1,n),Z_aft,given_weight,given_test,given_tol)
dataset_W_aft_NC=sample_path_W(beta_hat_aft,T_s_aft,rep(1,n),Z_aft,given_weight,given_test,dataset_What_aft_NC)

kol_typ_test_aft_NC=kolmogorov(dataset_W_aft_NC,dataset_What_aft_NC);kol_typ_test_aft_NC

p_aft_NC=kol_typ_test_aft_NC[3,];p_aft_NC

# PLOT : W_aft_NC vs What_aft_NC
Figure1_W_aft_NC=
  ggplot()+
  geom_line(data=dataset_What_aft_NC,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
  geom_line(data=dataset_W_aft_NC,aes(x=t_i,y=W),colour="tomato")
Figure1_W_aft_NC

# PLOT : std.W_aft_NC vs std.What_aft_NC
Figure1_std.W_aft_NC=
  ggplot()+
  geom_line(data=dataset_What_aft_NC,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
  geom_line(data=dataset_W_aft_NC,aes(x=t_i,y=std.W),colour="tomato")
Figure1_std.W_aft_NC

#-------------------------------------------------------------
#----------------------FUNCTION FORM TEST---------------------
#-------------------------------------------------------------
# dataset_W()
dataset_What_aft_f=sample_path_What(path,beta_hat_aft_f,T_aft_f,D_aft_f,Z_aft_f,given_weight,given_test,given_tol)

# dataset_What()
dataset_W_aft_f=sample_path_W(beta_hat_aft_f,T_aft_f,D_aft_f,Z_aft_f,given_weight,given_test,dataset_What_aft_f)

dataset_W_aft_f_order=dataset_W_aft_f$W[order(Z_aft_f)]
Z_aft_f_order=Z_aft_f[order(Z_aft_f)]

dataset_W_aft_order=dataset_W_aft$W[order(Z_aft)]
Z_aft_order=Z_aft[order(Z_aft)]

plot(Z_aft_f_order,dataset_W_aft_f_order,col="red",ylim = c(-0.5,0.5),type="l")
par(new=TRUE)
plot(Z_aft_order,dataset_W_aft_order,col="black",ylim = c(-0.5,0.5),type="l")


dataset_std.W_aft_f_order=dataset_W_aft_f$std.W[order(Z_aft_f)]
Z_aft_f_order=Z_aft_f[order(Z_aft_f)]

dataset_std.W_aft_order=dataset_W_aft$std.W[order(Z_aft)]
Z_aft_order=Z_aft[order(Z_aft)]

plot(Z_aft_f_order,dataset_std.W_aft_f_order,col="red",ylim=c(-1.5,2.5),type="l")
par(new=TRUE)
plot(Z_aft_order,dataset_std.W_aft_order,col="black",ylim=c(-1.5,2.5),type="l")
