#-------------------------------------------------------------
#---------------------------SETTING---------------------------
#-------------------------------------------------------------
rm(list=ls())

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
  Time=Time[order(e_i_beta)]
  Covari=Covari[order(e_i_beta)]
  Delta=Delta[order(e_i_beta)]
  e_i_beta=e_i_beta[order(e_i_beta)]
  
  if (weight=="a"){w_i=Covari*(Covari<=median(Covari))}
  if (weight=="b"){w_i=Covari}
  if (weight=="c"){w_i=1*(Covari<=median(Covari))}
  if (weight=="d"){w_i=1}
  
  N_i_s_t.beta=list()
  for(j in 1:n){
    N_i_s_t.beta[[j]]=(e_i_beta>=e_i_beta[j])*Delta[j]
  }
  #N_i_s_t.beta
  
  Y_i_s_t.beta=list()
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
  
  if (test=="omni"){return(W_t)}
  if (test=="ftn.form"){return(W_t[order(Covari)])}
  
}
#W_t()

#-------------------------------------------------------------
#-----------------------AFT DATA ANALYSIS---------------------
#-------------------------------------------------------------
What_t=function(b,Time,Delta,Covari,weight,tol,test){
  #b=beta_hat_aft;weight=given_weight;Time=T_aft;Delta=D_aft;Covari=Z_aft;tol=given_tol;test=given_test;
  
  n=length(Time)
  
  e_i_beta=log(Time)+Covari*b
  Time=Time[order(e_i_beta)]
  Covari=Covari[order(e_i_beta)]
  Delta=Delta[order(e_i_beta)]
  e_i_beta=e_i_beta[order(e_i_beta)]
  
  if (weight=="a"){w_i=Covari*(Covari<=median(Covari))}
  if (weight=="b"){w_i=Covari}
  if (weight=="c"){w_i=1*(Covari<=median(Covari))}
  if (weight=="d"){w_i=1}
  
  N_i_s_t.beta=list()
  for(j in 1:n){
    N_i_s_t.beta[[j]]=(e_i_beta>=e_i_beta[j])*Delta[j]
  }
  #N_i_s_t.beta
  
  Y_i_s_t.beta=list()
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
  
  #Q_i_t=S_0_s_t.beta/n # Gehan's weight
  Q_i_t=1
  
  U_t.beta=Reduce('+',lapply(mapply('*',dN_i_s_t.beta,
                                    Q_i_t*(Covari-E_s_t.beta), SIMPLIFY = FALSE),cumsum))
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
  b_g0=b
  
  e_i_beta_g0=log(Time)+Covari*b_g0
  Time_g0=Time[order(e_i_beta_g0)]
  Covari_g0=Covari[order(e_i_beta_g0)]
  Delta_g0=Delta[order(e_i_beta_g0)]
  e_i_beta_g0=e_i_beta_g0[order(e_i_beta_g0)]
  
  N_i_s_t.beta_g0=list()
  for(j in 1:n){
    N_i_s_t.beta_g0[[j]]=(e_i_beta_g0>=e_i_beta_g0[j])*Delta_g0[j]
  }
  #N_i_s_t.beta_g0
  
  Y_i_s_t.beta_g0=list()
  for(j in 1:n){
    Y_i_s_t.beta_g0[[j]]=(e_i_beta_g0<=e_i_beta_g0[j])*1
  }
  #Y_i_s_t.beta_g0
  
  N_d_s_t.beta_g0=Reduce('+',N_i_s_t.beta_g0)
  #N_d_s_t.beta_g0
  
  S_0_s_t.beta_g0=Reduce('+',Y_i_s_t.beta_g0)
  #S_0_s_t.beta_g0
  
  S_1_s_t.beta_g0=Reduce('+',mapply(
    "*", Y_i_s_t.beta_g0, Covari_g0, SIMPLIFY = FALSE))
  #S_1_s_t.beta_g0
  
  E_s_t.beta_g0=S_1_s_t.beta_g0/S_0_s_t.beta_g0
  #E_s_t.beta_g0
  
  J_t.beta_g0=(S_0_s_t.beta_g0>0)*1
  #J_t.beta_g0
  
  dN_d_s_t.beta_g0=diff(c(0,N_d_s_t.beta_g0))
  #dN_d_s_t.beta_g0
  
  Ahat_0_t.beta_g0=cumsum((J_t.beta_g0/S_0_s_t.beta_g0)*dN_d_s_t.beta_g0)
  #Ahat_0_t.beta_g0
  
  Ghat_0_t=1-exp(-Ahat_0_t.beta_g0)
  #Ghat_0_t
  
  dGhat_0_t=diff(c(0,Ghat_0_t))
  #dGhat_0_t
  
  ghat_0_t=(ksmooth(Time_g0,dGhat_0_t,"normal",
                    bandwidth = 1.06*sd(dGhat_0_t)*n^(-0.2),x.points=Time_g0)$y)
  #ghat_0_t
  
  fhat_Y_t=ghat_0_t*Time*(sum(w_i*Covari)/n)
  #fhat_Y_t
  
  b_f0=b
  
  e_i_beta_f0=log(Time)+Covari*b_f0
  Time_f0=Time[order(e_i_beta_f0)]
  Covari_f0=Covari[order(e_i_beta_f0)]
  Delta_f0=Delta[order(e_i_beta_f0)]
  e_i_beta_f0=e_i_beta_f0[order(e_i_beta_f0)]
  
  N_i_e_f0=list()
  for(j in 1:n){
    N_i_e_f0[[j]]=(e_i_beta_f0>=e_i_beta_f0[j])*Delta_f0[j]
  }
  #N_i_e_f0
  
  Y_i_e_f0=list()
  for(j in 1:n){
    Y_i_e_f0[[j]]=(e_i_beta_f0<=e_i_beta_f0[j])*1
  }
  #Y_i_e_f0
  
  N_d_e_f0=Reduce('+',N_i_e_f0)
  #N_d_e_f0
  
  Y_d_e_f0=Reduce('+',Y_i_e_f0)
  #Y_d_e_f0
  
  KM_e_f0=cumprod(1-Delta_f0/Y_d_e_f0)
  #KM_e_f0
  
  Fhat_0_e_f0=1-KM_e_f0
  #Fhat_0_e_f0
  
  dFhat_0_e_f0=diff(c(0,Fhat_0_e_f0))
  #dFhat_0_e_f0
  
  Condi.Ehat_f0=(cumsum(e_i_beta_f0*dFhat_0_e_f0))/(1-Fhat_0_e_f0)
  #Condi.Ehat_f0
  
  rhat_i_f0=Delta_f0*e_i_beta_f0+(1-Delta_f0)*Condi.Ehat_f0
  #rhat_i_f0
  
  fhat_0_t=ksmooth(Time_f0,dFhat_0_e_f0,"normal",
                   bandwidth = 1.06*sd(dFhat_0_e_f0)*n^(-0.2),x.points=Time_f0)$y
  #fhat_0_t
  
  fhat_N_t=fhat_0_t*Time*(sum(Delta*w_i*Covari)/n)
  #fhat_N_t
  
  #-----------------------------------------------------------
  #--------Find Beta_hat_star by using optimize function------
  #-----------------------------------------------------------
  U_beta=function(beta_U=b,Time_U=Time,Delta_U=Delta,Covari_U=Covari){
    #beta_U=b;Time_U=Time;Delta_U=Delta;Covari_U=Covari;Q_i_t_U=Q_i_t;
    
    e_i_beta_U=log(Time_U)+Covari_U*beta_U
    Time_U=Time_U[order(e_i_beta_U)]
    Covari_U=Covari_U[order(e_i_beta_U)]
    Delta_U=Delta_U[order(e_i_beta_U)]
    e_i_beta_U=e_i_beta_U[order(e_i_beta_U)]
    
    N_i_s_t.beta_U=list()
    for(j in 1:n){
      N_i_s_t.beta_U[[j]]=(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]
    }
    #N_i_s_t.beta_U
    
    Y_i_s_t.beta_U=list()
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
    
    Q_i_t_U=S_0_s_t.beta_U/n # Gehan's weight
    #Q_i_t=1
    
    U_t.beta_U=Reduce('+',lapply(mapply('*',dN_i_s_t.beta_U,
                                        Q_i_t_U*(Covari_U-E_s_t.beta_U),SIMPLIFY = FALSE),cumsum))
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
    
    U_w_G_t.beta=Reduce('+',mapply('*',lapply(mapply('*',dMhat_i_s_t.beta,
                                                     Q_i_t*(w_i-E_w_s_t.beta),SIMPLIFY = FALSE),cumsum),G_i,SIMPLIFY = FALSE))
    #U_w_G_t.beta
    
    U_G_t.beta=Reduce('+',mapply('*',lapply(mapply('*',dMhat_i_s_t.beta,
                                                   Q_i_t*(Covari-E_s_t.beta),SIMPLIFY = FALSE),cumsum),G_i,SIMPLIFY = FALSE))
    #U_G_t.beta
    
    U_G_t.inf.beta=U_G_t.beta[n]
    #U_G_t.inf.beta
    
    beta_hat_s_list=optimize(function(beta){abs(U_beta(beta)-U_G_t.inf.beta)},
                             c(b-std_hat_aft,b+std_hat_aft),
                             tol = 1e-16)
    #beta_hat_s_list
    
    tolerance=beta_hat_s_list$objective
    #tolerance
  }
  
  beta_hat_s=beta_hat_s_list$minimum;beta_hat_s
  #beta_hat_s
  
  e_i_beta_s=log(Time)+Covari*beta_hat_s
  Time_s=Time[order(e_i_beta_s)]
  Covari_s=Covari[order(e_i_beta_s)]
  Delta_s=Delta[order(e_i_beta_s)]
  e_i_beta_s=e_i_beta_s[order(e_i_beta_s)]
  
  N_i_s_t.beta_s=list()
  for(j in 1:n){
    N_i_s_t.beta_s[[j]]=(e_i_beta_s>=e_i_beta_s[j])*Delta_s[j]
  }
  #N_i_s_t.beta_s
  
  Y_i_s_t.beta_s=list()
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
  BB=sqrt(n)*(fhat_N_t+cumsum(fhat_N_t*dAhat_0_t.beta))*(b-beta_hat_s);BB
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
  
  if (test=="omni"){return(What_t)}
  if (test=="ftn.form"){return(What_t[order(Covari)])}
  
}
#What_t()

simulation_What=function(sim,b,Time,Delta,Covari,weight,test,tol){
  
  dataset_What=data.frame()
  for(k in 1:sim){
    group=k
    A=What_t(b,Time,Delta,Covari,weight,tol,test)
    AA=data.frame(group,t_i=1:n,What=A)
    dataset_What=rbind(dataset_What,AA)
    if(k%%100==0) {
      cat("Simulation",k,"\n")
    }
  }
  
  bootstrap=function(dataset_What){
    at.time.t=list()
    std.at.time.t=vector()
    for(k in 1:n){
      at.time.t[[k]]=data.frame(as.matrix(dataset_What)[which(dataset_What$t_i==k),])
      std.at.time.t[k]=sd(sample(at.time.t[[k]]$What,n,replace=TRUE))
    }
    return(std.at.time.t)
  }
  
  mean.std.boot=bootstrap(dataset_What)
  
  AAA=data.frame((dataset_What$What/mean.std.boot),mean.std.boot)
  dataset_What=cbind(dataset_What,AAA)
  colnames(dataset_What)=c("group","t_i","What","std.What","std.boot")
  
  return(dataset_What)
}

simulation_W=function(b,Time,Delta,Covari,weight,test,dataset_What){
  mean.std.boot=dataset_What$std.boot[1:n]
  B=W_t(b,Time,Delta,Covari,weight,test)
  dataset_W=data.frame(0,1:n,B)
  AAAA=data.frame((dataset_W[3]/mean.std.boot),mean.std.boot)
  dataset_W=cbind(dataset_W,AAAA)
  colnames(dataset_W)=c("group","t_i","W","std.W","std.boot")
  return(dataset_W)
}

# KOLMOGOROV TYPE SUPREMUM TEST
kolmogorov=function(dataset_W,dataset_What){
  
  sim=length(dataset_What[,1])/length(dataset_W[,1])
  
  # kolmogorov type supremum test for W
  max_W=max(abs(dataset_W$W))
  max_What=c()
  for(k in 1: sim){
    max_What[k]=max(abs(dataset_What[which(dataset_What$group==k),]$What))
  }
  n_W=length(which((max_What>max_W)==1))
  p_W=n_W/sim
  
  # kolmogorov type supremum test for std.W
  max_std_W=max(abs(dataset_W$std.W))
  max_std_What=c()
  for(k in 1: sim){
    max_std_What[k]=max(abs(dataset_What[which(dataset_What$group==k),]$std.What))
  }
  n_std.W=length(which((max_std_What>max_std_W)==1))
  p_std.W=n_std.W/sim
  
  # kolmogorov type supremum test result
  test_result=rbind(c(sim,sim),c(n_W,n_std.W),c(p_W,p_std.W))
  colnames(test_result)=c("W","std.W")
  rownames(test_result)=c("sim","n","p")
  return(test_result)
}

#-------------------------------------------------------------
#------------------------DATA GENERATE------------------------
#-------------------------------------------------------------
n=20
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

sim=200

#-------------------------------------------------------------
#-------------------------OMNIBUS TEST------------------------
#-------------------------------------------------------------

#-------------------------NONCENSORING------------------------
dataset_What_aft_NC=simulation_What(sim,beta_hat_aft,T_s_aft,rep(1,n),Z_aft,given_weight,given_test,given_tol)
dataset_W_aft_NC=simulation_W(beta_hat_aft,T_s_aft,rep(1,n),Z_aft,given_weight,given_test,dataset_What_aft_NC)

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

#--------------------------CENSORING--------------------------
# dataset_What()
dataset_What_aft=simulation_What(sim,beta_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,given_tol)
dataset_What_cox=simulation_What(sim,beta_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,given_tol)

# dataset_W()
dataset_W_aft=simulation_W(beta_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,dataset_What_aft)
dataset_W_cox=simulation_W(beta_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,dataset_What_cox)

kol_typ_test_aft=kolmogorov(dataset_W_aft,dataset_What_aft);kol_typ_test_aft
kol_typ_test_cox=kolmogorov(dataset_W_cox,dataset_What_cox);kol_typ_test_cox

p_aft=kol_typ_test_aft[3,];p_aft
p_cox=kol_typ_test_cox[3,];p_cox

# dataset_What50_aft=dataset_What_aft[1:(n*50),]
# dataset_What50_cox=dataset_What_cox[1:(n*50),]

# PLOT : W_aft vs What_aft
Figure1_W_aft=
  ggplot()+
  geom_line(data=dataset_What_aft,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
  geom_line(data=dataset_W_aft,aes(x=t_i,y=W),colour="tomato")
Figure1_W_aft

# PLOT : W_cox vs What_cox
Figure1_W_cox=
  ggplot()+
  geom_line(data=dataset_What_cox,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
  geom_line(data=dataset_W_cox,aes(x=t_i,y=W),colour="tomato")
Figure1_W_cox

# PLOT : std.W_aft vs std.What_aft
Figure1_std.W_aft=
  ggplot()+
  geom_line(data=dataset_What_aft,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
  geom_line(data=dataset_W_aft,aes(x=t_i,y=std.W),colour="tomato")
Figure1_std.W_aft.aft

# PLOT : W_cox vs What_cox
Figure1_std.W_cox=
  ggplot()+
  geom_line(data=dataset_What_cox,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
  geom_line(data=dataset_W_cox,aes(x=t_i,y=std.W),colour="tomato")
Figure1_std.W_cox

#-------------------------------------------------------------
#----------------------FUNCTION FORM TEST---------------------
#-------------------------------------------------------------
# dataset_W()
dataset_What_aft_f=simulation_What(sim,beta_hat_aft_f,T_aft_f,D_aft_f,Z_aft_f,given_weight,given_test,given_tol)

# dataset_What()
dataset_W_aft_f=simulation_W(beta_hat_aft_f,T_aft_f,D_aft_f,Z_aft_f,given_weight,given_test,dataset_What_aft_f)

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


#######################simulation&bootstarpping code need to be chainged
dataset_What=data.frame(What_t(b,Time,Delta,Covari,weight,tol,test))
columnnames[1]=paste("What",1)
for(k in 2:sim){
  
  dataset_What=cbind(dataset_What,What_t(b,Time,Delta,Covari,weight,tol,test))
  columnnames[k]=paste("What",k)
  colnames(dataset_What)=columnnames
  
  if(k%%100==0) {
    cat("Simulation",k,"\n")
  }
}

sd(sample(dataset_What[1,],sim,replace=TRUE))

