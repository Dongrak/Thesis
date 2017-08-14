<<<<<<< HEAD
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
#------------------------WEIGHT&TOLERANCE---------------------
#-------------------------------------------------------------
given_tol=0.1

given_weight="c"
#(weight=="a"){w_i=Covari*(Covari<=median(Covari))}
#(weight=="b"){w_i=Covari}  
#(weight=="c"){w_i=1*(Covari<=median(Covari))}
#(weight=="d"){w_i=1}

#-------------------------------------------------------------
#-----------------------COX DATA ANALYSIS---------------------
#-------------------------------------------------------------
W_t=function(b=beta_hat_cox,weight=given_weight,Time=T_cox,Delta=D_cox,Covari=Z_cox){
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
  
  Y_d_s_t.beta=Reduce('+',Y_i_s_t.beta)
  #Y_d_s_t.beta
  
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
  
  return(W_t)
}
#W_t()

#-------------------------------------------------------------
#-----------------------AFT DATA ANALYSIS---------------------
#-------------------------------------------------------------
What_t=function(b=beta_hat_aft,weight=given_weight,Time=T_aft,Delta=D_aft,Covari=Z_aft,tol=given_tol){
  #b=beta_hat_aft;weight=given_weight;Time=T_aft;Delta=D_aft;Covari=Z_aft;tol=given_tolerance;
  
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
  
  Y_d_s_t.beta=Reduce('+',Y_i_s_t.beta)
  #Y_d_s_t.beta
  
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
  
  Y_d_s_t.beta_g0=Reduce('+',Y_i_s_t.beta_g0)
  #Y_d_s_t.beta_g0
  
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
    
    Y_d_s_t.beta_U=Reduce('+',Y_i_s_t.beta_U)
    #Y_d_s_t.beta_U
    
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
  return(What_t)
  #return(list(What_t=What_t,beta_hat_s=beta_hat_s))
}
#What_t()

simulation_What=function(sim=50,b=beta_hat_aft,weight=given_weight,Time=T_aft,Delta=D_aft,Covari=Z_aft,tol=given_tolerance){
  
  dataset_What=data.frame()
  for(k in 1:sim){
    group=k
    A=What_t(b,weight,Time,Delta,Covari)
    AA=data.frame(group,t_i=1:n,What=A)
    dataset_What=rbind(dataset_What,AA)
    if(k%%100==0) {
      cat("Simulation",k,"\n")
    }
  }
  
  bootstrap=function(){
    at.time.t=list()
    std.at.time.t=vector()
    for(k in 1:n){
      at.time.t[[k]]=data.frame(as.matrix(dataset_What)[which(dataset_What$t_i==k),])
      std.at.time.t[k]=sd(at.time.t[[k]]$What)
    }
    return(std.at.time.t)
  }
  
  mean.std.boot=bootstrap()
  
  AAA=data.frame((dataset_What$What/mean.std.boot),mean.std.boot)
  dataset_What=cbind(dataset_What,AAA)
  colnames(dataset_What)=c("group","t_i","What","std.What","std.boot")
  
  return(dataset_What)
}

simulation_W=function(b=beta_hat_aft,weight=given_weight,Time=T_aft,Delta=D_aft,Covari=Z_aft,dataset_What){
  mean.std.boot=dataset_What$std.boot[1:n]
  B=W_t(b,weight,Time,Delta,Covari)
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
    max_What[k]=max(abs(as.matrix(dataset_What)[which(dataset_What$group==k),3]))
  }
  n_W=length(which((max_What>max_W)==1))
  p_W=n_W/sim
  
  # kolmogorov type supremum test for std.W
  max_std_W=max(abs(dataset_W$std.W))
  max_std_What=c()
  for(k in 1: sim){
    max_std_What[k]=max(abs(as.matrix(dataset_What)[which(dataset_What$group==k),4]))
  }
  n_std.W=length(which((max_std_What>max_std_W)==1))
  p_std.W=n_std.W/sim
  
  # kolmogorov type supremum test result
  test_result=rbind(c(sim,sim),c(n_W,n_std.W),c(p_W,p_std.W))
  colnames(test_result)=c("W","std.W")
  rownames(test_result)=c("sim","n","p")
  return(test_result)
}

n_vector=c(100,200,500,1000,2000)
sim_vector=c(200,500,1000,2000)

#n_vector=c(10,20,30,40,50)
#sim_vector=c(15,25,35,45,55)

result_aft.aft=list()
result_aft.cox=list()
result_cox.aft=list()
result_cox.cox=list()

for(k in 1:length(n_vector)){
  cat("n_vector : ",k,"\n")
  
  n=n_vector[k]
  sim=max(sim_vector)
  
  #-------------------------------------------------------------
  #------------------------DATA GENERATE------------------------
  #-------------------------------------------------------------
  id=c(1:n) # identification
  beta_0=1 # beta_0
  Z=rnorm(n,3,1) # covariate
  
  #-----------------------AFT DATA GENERATE---------------------
  T_s_aft=exp(-beta_0*Z)*qlnorm(runif(n),5,1) # lognormal baseline hazard
  C_aft=exp(-beta_0*Z)*qlnorm(runif(n),6.5,1) # censoring time for aft model
  T_aft=C_aft*(T_s_aft>C_aft)+T_s_aft*(T_s_aft<=C_aft) # observed time for aft model
  D_aft=0*(T_s_aft>C_aft)+1*(T_s_aft<=C_aft) # delta
  D_aft=D_aft[order(T_aft)]
  Z_aft=Z[order(T_aft)]
  T_aft=T_aft[order(T_aft)]
  
  #-----------------------COX DATA GENERATE---------------------
  T_s_cox=qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5,1,lower.tail = FALSE) # lognormal baseline hazard
  C_cox=qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5.7,1,lower.tail = FALSE) # censoring time for cox model
  T_cox=C_cox*(T_s_cox > C_cox)+T_s_cox*(T_s_cox<=C_cox) # observed time for cox model
  D_cox=0*(T_s_cox > C_cox)+1*(T_s_cox <= C_cox) # delta
  D_cox=D_cox[order(T_cox)]
  Z_cox=Z[order(T_cox)]
  T_cox=T_cox[order(T_cox)]
  
  #-------------------------------------------------------------
  #--------------Estimate Beta_hat_aft by using Aftgee--------------
  #-------------------------------------------------------------
  aftsrr_beta_aft=aftsrr(Surv(T_aft,D_aft)~Z_aft,method="nonsm")
  beta_hat_aft=-unlist(summary(aftsrr_beta_aft))$coefficients1;beta_hat_aft
  std_hat_aft=unlist(summary(aftsrr_beta_aft))$coefficients2;std_hat_aft
  
  #-------------------------------------------------------------
  #--------------Estimate Beta_hat_cox by using Aftgee--------------
  #-------------------------------------------------------------
  aftsrr_beta_cox=aftsrr(Surv(T_cox,D_cox)~Z_cox,method="nonsm")
  beta_hat_cox=-unlist(summary(aftsrr_beta_cox))$coefficients1;beta_hat_cox
  std_hat_cox=unlist(summary(aftsrr_beta_cox))$coefficients2;std_hat_cox
  
  #dataset_What()
  dataset_What_aft=simulation_What(sim,beta_hat_aft,given_weight,T_aft,D_aft,Z_aft)
  dataset_What_cox=simulation_What(sim,beta_hat_cox,given_weight,T_cox,D_aft,Z_cox)
  
  #dataset_W()
  dataset_W_aft=simulation_W(beta_hat_aft,given_weight,T_aft,D_aft,Z_aft,dataset_What_aft)
  dataset_W_cox=simulation_W(beta_hat_cox,given_weight,T_cox,D_aft,Z_cox,dataset_What_cox)
  
  dataset_What50_aft=dataset_What_aft[1:(n*50),]
  dataset_What50_cox=dataset_What_cox[1:(n*50),]
  
  # PLOT : W_aft vs What_aft
  Figure1_W_aft.aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_aft.aft
  
  # PLOT : W_cox vs What_aft
  Figure1_W_cox.aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_cox.aft
  
  # PLOT : W_aft vs What_cox
  Figure1_W_aft.cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_aft.cox
  
  # PLOT : W_cox vs What_cox
  Figure1_W_cox.cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_cox.cox
  
  # PLOT : std.W_aft vs std.What_aft
  Figure1_std.W_aft.aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_aft.aft
  
  # PLOT : std.W_cox vs std.What_aft
  Figure1_std.W_cox.aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_cox.aft
  
  # PLOT : std.W_aft vs std.What_cox
  Figure1_std.W_aft.cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_aft.cox
  
  # PLOT : W_cox vs What_cox
  Figure1_std.W_cox.cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_cox.cox

  dataset_What200_aft=dataset_What_aft[1:(n*sim_vector[1]),]
  dataset_What200_cox=dataset_What_cox[1:(n*sim_vector[1]),]
  
  dataset_What500_aft=dataset_What_aft[1:(n*sim_vector[2]),]
  dataset_What500_cox=dataset_What_cox[1:(n*sim_vector[2]),]
  
  dataset_What1000_aft=dataset_What_aft[1:(n*sim_vector[3]),]
  dataset_What1000_cox=dataset_What_cox[1:(n*sim_vector[3]),]
  
  dataset_What2000_aft=dataset_What_aft[1:(n*sim_vector[4]),]
  dataset_What2000_cox=dataset_What_cox[1:(n*sim_vector[4]),]
  
  result_aft.aft[[k]]=list(Figure1_W_aft.aft,Figure1_std.W_aft.aft,
                           dataset_W_aft,dataset_What_aft,
                           kolmogorov(dataset_W_aft,dataset_What200_aft),
                           kolmogorov(dataset_W_aft,dataset_What500_aft),
                           kolmogorov(dataset_W_aft,dataset_What1000_aft),
                           kolmogorov(dataset_W_aft,dataset_What2000_aft))
  
  result_aft.cox[[k]]=list(Figure1_W_aft.cox,Figure1_std.W_aft.cox,
                           dataset_W_aft,dataset_What_cox,
                           kolmogorov(dataset_W_aft,dataset_What200_cox),
                           kolmogorov(dataset_W_aft,dataset_What500_cox),
                           kolmogorov(dataset_W_aft,dataset_What1000_cox),
                           kolmogorov(dataset_W_aft,dataset_What2000_cox))
  
  result_cox.aft[[k]]=list(Figure1_W_cox.aft,Figure1_std.W_cox.aft,
                           dataset_W_cox,dataset_What_aft,
                           kolmogorov(dataset_W_cox,dataset_What200_aft),
                           kolmogorov(dataset_W_cox,dataset_What500_aft),
                           kolmogorov(dataset_W_cox,dataset_What1000_aft),
                           kolmogorov(dataset_W_cox,dataset_What2000_aft))
  
  result_cox.cox[[k]]=list(Figure1_W_cox.cox,Figure1_std.W_cox.cox,
                           dataset_W_cox,dataset_What_cox,
                           kolmogorov(dataset_W_cox,dataset_What200_cox),
                           kolmogorov(dataset_W_cox,dataset_What500_cox),
                           kolmogorov(dataset_W_cox,dataset_What1000_cox),
                           kolmogorov(dataset_W_cox,dataset_What2000_cox))
}

#n_vector=c(100,200,500,1000,2000)
#sim_vector=c(200,500,1000,2000)

table_fuction=function(result){
hh=matrix(nrow=length(n_vector),ncol=length(sim_vector))
rownames(hh)=c("n=100","n=200","n=500","n=1000","n=2000")
colnames(hh)=c("sim=200","sim=500","sim=1000","sim=2000")
for(k in 1:length(n_vector)){
  for(j in 1:length(sim_vector)){
    n=n_vector[k]
    sim=sim_vector[j]
    
    aa=result[[k]]
    i=4+j
    bb=aa[[i]] #########################위 코드 변경으로 나중에는 1- => 없애기 
    cc=sprintf("%1.3f",bb[3,1])
    dd=sprintf("%1.3f",bb[3,2])
    ee=paste(cc,dd)
    hh[k,j]=ee
    
  }
}
return(hh)
}

table_aft.aft=table_fuction(result_aft.aft);table_aft.aft
#table_aft.cox=table_fuction(result_aft.cox);table_aft.cox
#table_cox.aft=table_fuction(result_cox.aft);table_cox.aft
table_cox.cox=table_fuction(result_cox.cox);table_cox.cox

result_aft.aft[[5]][2]
#result_cox.aft[[1]][1]
#result_aft.cox[[1]][1]
result_cox.cox[[5]][2]
=======
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
#------------------------WEIGHT&TOLERANCE---------------------
#-------------------------------------------------------------
given_tol=0.1

given_weight="c"
#(weight=="a"){w_i=Covari*(Covari<=median(Covari))}
#(weight=="b"){w_i=Covari}  
#(weight=="c"){w_i=1*(Covari<=median(Covari))}
#(weight=="d"){w_i=1}

#-------------------------------------------------------------
#-----------------------COX DATA ANALYSIS---------------------
#-------------------------------------------------------------
W_t=function(b=beta_hat_cox,weight=given_weight,Time=T_cox,Delta=D_cox,Covari=Z_cox){
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
  
  Y_d_s_t.beta=Reduce('+',Y_i_s_t.beta)
  #Y_d_s_t.beta
  
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
  
  return(W_t)
}
#W_t()

#-------------------------------------------------------------
#-----------------------AFT DATA ANALYSIS---------------------
#-------------------------------------------------------------
What_t=function(b=beta_hat_aft,weight=given_weight,Time=T_aft,Delta=D_aft,Covari=Z_aft,tol=given_tol){
  #b=beta_hat_aft;weight=given_weight;Time=T_aft;Delta=D_aft;Covari=Z_aft;tol=given_tolerance;
  
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
  
  Y_d_s_t.beta=Reduce('+',Y_i_s_t.beta)
  #Y_d_s_t.beta
  
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
  
  Y_d_s_t.beta_g0=Reduce('+',Y_i_s_t.beta_g0)
  #Y_d_s_t.beta_g0
  
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
    
    Y_d_s_t.beta_U=Reduce('+',Y_i_s_t.beta_U)
    #Y_d_s_t.beta_U
    
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
  return(What_t)
  #return(list(What_t=What_t,beta_hat_s=beta_hat_s))
}
#What_t()

simulation_What=function(sim=50,b=beta_hat_aft,weight=given_weight,Time=T_aft,Delta=D_aft,Covari=Z_aft,tol=given_tolerance){
  
  dataset_What=data.frame()
  for(k in 1:sim){
    group=k
    A=What_t(b,weight,Time,Delta,Covari)
    AA=data.frame(group,t_i=1:n,What=A)
    dataset_What=rbind(dataset_What,AA)
    if(k%%100==0) {
      cat("Simulation",k,"\n")
    }
  }
  
  bootstrap=function(){
    at.time.t=list()
    std.at.time.t=vector()
    for(k in 1:n){
      at.time.t[[k]]=data.frame(as.matrix(dataset_What)[which(dataset_What$t_i==k),])
      std.at.time.t[k]=sd(at.time.t[[k]]$What)
    }
    return(std.at.time.t)
  }
  
  mean.std.boot=bootstrap()
  
  AAA=data.frame((dataset_What$What/mean.std.boot),mean.std.boot)
  dataset_What=cbind(dataset_What,AAA)
  colnames(dataset_What)=c("group","t_i","What","std.What","std.boot")
  
  return(dataset_What)
}

simulation_W=function(b=beta_hat_aft,weight=given_weight,Time=T_aft,Delta=D_aft,Covari=Z_aft,dataset_What){
  mean.std.boot=dataset_What$std.boot[1:n]
  B=W_t(b,weight,Time,Delta,Covari)
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
    max_What[k]=max(abs(as.matrix(dataset_What)[which(dataset_What$group==k),3]))
  }
  n_W=length(which((max_What>max_W)==1))
  p_W=n_W/sim
  
  # kolmogorov type supremum test for std.W
  max_std_W=max(abs(dataset_W$std.W))
  max_std_What=c()
  for(k in 1: sim){
    max_std_What[k]=max(abs(as.matrix(dataset_What)[which(dataset_What$group==k),4]))
  }
  n_std.W=length(which((max_std_What>max_std_W)==1))
  p_std.W=n_std.W/sim
  
  # kolmogorov type supremum test result
  test_result=rbind(c(sim,sim),c(n_W,n_std.W),c(p_W,p_std.W))
  colnames(test_result)=c("W","std.W")
  rownames(test_result)=c("sim","n","p")
  return(test_result)
}

n_vector=c(100,200,500,1000,2000)
sim_vector=c(200,500,1000,2000)

#n_vector=c(10,20,30,40,50)
#sim_vector=c(15,25,35,45,55)

result_aft.aft=list()
result_aft.cox=list()
result_cox.aft=list()
result_cox.cox=list()

for(k in 1:length(n_vector)){
  cat("n_vector : ",k,"\n")
  
  n=n_vector[k]
  sim=max(sim_vector)
  
  #-------------------------------------------------------------
  #------------------------DATA GENERATE------------------------
  #-------------------------------------------------------------
  id=c(1:n) # identification
  beta_0=1 # beta_0
  Z=rnorm(n,3,1) # covariate
  
  #-----------------------AFT DATA GENERATE---------------------
  T_s_aft=exp(-beta_0*Z)*qlnorm(runif(n),5,1) # lognormal baseline hazard
  C_aft=exp(-beta_0*Z)*qlnorm(runif(n),6.5,1) # censoring time for aft model
  T_aft=C_aft*(T_s_aft>C_aft)+T_s_aft*(T_s_aft<=C_aft) # observed time for aft model
  D_aft=0*(T_s_aft>C_aft)+1*(T_s_aft<=C_aft) # delta
  D_aft=D_aft[order(T_aft)]
  Z_aft=Z[order(T_aft)]
  T_aft=T_aft[order(T_aft)]
  
  #-----------------------COX DATA GENERATE---------------------
  T_s_cox=qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5,1,lower.tail = FALSE) # lognormal baseline hazard
  C_cox=qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5.7,1,lower.tail = FALSE) # censoring time for cox model
  T_cox=C_cox*(T_s_cox > C_cox)+T_s_cox*(T_s_cox<=C_cox) # observed time for cox model
  D_cox=0*(T_s_cox > C_cox)+1*(T_s_cox <= C_cox) # delta
  D_cox=D_cox[order(T_cox)]
  Z_cox=Z[order(T_cox)]
  T_cox=T_cox[order(T_cox)]
  
  #-------------------------------------------------------------
  #--------------Estimate Beta_hat_aft by using Aftgee--------------
  #-------------------------------------------------------------
  aftsrr_beta_aft=aftsrr(Surv(T_aft,D_aft)~Z_aft,method="nonsm")
  beta_hat_aft=-unlist(summary(aftsrr_beta_aft))$coefficients1;beta_hat_aft
  std_hat_aft=unlist(summary(aftsrr_beta_aft))$coefficients2;std_hat_aft
  
  #-------------------------------------------------------------
  #--------------Estimate Beta_hat_cox by using Aftgee--------------
  #-------------------------------------------------------------
  aftsrr_beta_cox=aftsrr(Surv(T_cox,D_cox)~Z_cox,method="nonsm")
  beta_hat_cox=-unlist(summary(aftsrr_beta_cox))$coefficients1;beta_hat_cox
  std_hat_cox=unlist(summary(aftsrr_beta_cox))$coefficients2;std_hat_cox
  
  #dataset_What()
  dataset_What_aft=simulation_What(sim,beta_hat_aft,given_weight,T_aft,D_aft,Z_aft)
  dataset_What_cox=simulation_What(sim,beta_hat_cox,given_weight,T_cox,D_aft,Z_cox)
  
  #dataset_W()
  dataset_W_aft=simulation_W(beta_hat_aft,given_weight,T_aft,D_aft,Z_aft,dataset_What_aft)
  dataset_W_cox=simulation_W(beta_hat_cox,given_weight,T_cox,D_aft,Z_cox,dataset_What_cox)
  
  dataset_What50_aft=dataset_What_aft[1:(n*50),]
  dataset_What50_cox=dataset_What_cox[1:(n*50),]
  
  # PLOT : W_aft vs What_aft
  Figure1_W_aft.aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_aft.aft
  
  # PLOT : W_cox vs What_aft
  Figure1_W_cox.aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_cox.aft
  
  # PLOT : W_aft vs What_cox
  Figure1_W_aft.cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_aft.cox
  
  # PLOT : W_cox vs What_cox
  Figure1_W_cox.cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_cox.cox
  
  # PLOT : std.W_aft vs std.What_aft
  Figure1_std.W_aft.aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_aft.aft
  
  # PLOT : std.W_cox vs std.What_aft
  Figure1_std.W_cox.aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_cox.aft
  
  # PLOT : std.W_aft vs std.What_cox
  Figure1_std.W_aft.cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_aft.cox
  
  # PLOT : W_cox vs What_cox
  Figure1_std.W_cox.cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_cox.cox

  dataset_What200_aft=dataset_What_aft[1:(n*sim_vector[1]),]
  dataset_What200_cox=dataset_What_cox[1:(n*sim_vector[1]),]
  
  dataset_What500_aft=dataset_What_aft[1:(n*sim_vector[2]),]
  dataset_What500_cox=dataset_What_cox[1:(n*sim_vector[2]),]
  
  dataset_What1000_aft=dataset_What_aft[1:(n*sim_vector[3]),]
  dataset_What1000_cox=dataset_What_cox[1:(n*sim_vector[3]),]
  
  dataset_What2000_aft=dataset_What_aft[1:(n*sim_vector[4]),]
  dataset_What2000_cox=dataset_What_cox[1:(n*sim_vector[4]),]
  
  result_aft.aft[[k]]=list(Figure1_W_aft.aft,Figure1_std.W_aft.aft,
                           dataset_W_aft,dataset_What_aft,
                           kolmogorov(dataset_W_aft,dataset_What200_aft),
                           kolmogorov(dataset_W_aft,dataset_What500_aft),
                           kolmogorov(dataset_W_aft,dataset_What1000_aft),
                           kolmogorov(dataset_W_aft,dataset_What2000_aft))
  
  result_aft.cox[[k]]=list(Figure1_W_aft.cox,Figure1_std.W_aft.cox,
                           dataset_W_aft,dataset_What_cox,
                           kolmogorov(dataset_W_aft,dataset_What200_cox),
                           kolmogorov(dataset_W_aft,dataset_What500_cox),
                           kolmogorov(dataset_W_aft,dataset_What1000_cox),
                           kolmogorov(dataset_W_aft,dataset_What2000_cox))
  
  result_cox.aft[[k]]=list(Figure1_W_cox.aft,Figure1_std.W_cox.aft,
                           dataset_W_cox,dataset_What_aft,
                           kolmogorov(dataset_W_cox,dataset_What200_aft),
                           kolmogorov(dataset_W_cox,dataset_What500_aft),
                           kolmogorov(dataset_W_cox,dataset_What1000_aft),
                           kolmogorov(dataset_W_cox,dataset_What2000_aft))
  
  result_cox.cox[[k]]=list(Figure1_W_cox.cox,Figure1_std.W_cox.cox,
                           dataset_W_cox,dataset_What_cox,
                           kolmogorov(dataset_W_cox,dataset_What200_cox),
                           kolmogorov(dataset_W_cox,dataset_What500_cox),
                           kolmogorov(dataset_W_cox,dataset_What1000_cox),
                           kolmogorov(dataset_W_cox,dataset_What2000_cox))
}

#n_vector=c(100,200,500,1000,2000)
#sim_vector=c(200,500,1000,2000)

table_fuction=function(result){
hh=matrix(nrow=length(n_vector),ncol=length(sim_vector))
rownames(hh)=c("n=100","n=200","n=500","n=1000","n=2000")
colnames(hh)=c("sim=200","sim=500","sim=1000","sim=2000")
for(k in 1:length(n_vector)){
  for(j in 1:length(sim_vector)){
    n=n_vector[k]
    sim=sim_vector[j]
    
    aa=result[[k]]
    i=4+j
    bb=aa[[i]] #########################위 코드 변경으로 나중에는 1- => 없애기 
    cc=sprintf("%1.3f",bb[3,1])
    dd=sprintf("%1.3f",bb[3,2])
    ee=paste(cc,dd)
    hh[k,j]=ee
    
  }
}
return(hh)
}

table_aft.aft=table_fuction(result_aft.aft);table_aft.aft
#table_aft.cox=table_fuction(result_aft.cox);table_aft.cox
#table_cox.aft=table_fuction(result_cox.aft);table_cox.aft
table_cox.cox=table_fuction(result_cox.cox);table_cox.cox

result_aft.aft[[5]][2]
#result_cox.aft[[1]][1]
#result_aft.cox[[1]][1]
result_cox.cox[[5]][2]
>>>>>>> 47f262838fd808749c38e967afd5232cc54fea73
