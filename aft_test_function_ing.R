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
#install.packages("evd")
#install.packages("ENmisc")

library(ggplot2)
library(survival)
library(aftgee)
library(evd)
library(ENmisc)

#-------------------------------------------------------------
#------------------------DATA GENERATE------------------------
#-------------------------------------------------------------
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
path=200

given_tol=0.1

# given_weight="fi"
# (weight=="11"){w_ij.z=1*1} 
# (weight=="fi"){w_ij.z=f()*I()} # omnibus
# (weight=="f1"){w_ij.z=f()*1} # ph & aft
# (weight=="1i"){w_ij.z=1*I()} # ftnform & linkftn

# given_test="omni"
# given_test="ftnform"
# given_test="linkftn"
# given_test="aft"

#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------
W_t.z_omni=function(b,Time,Delta,Covari,weight="fi"){
  #b=beta_hat_cox;Time=T_cox;Delta=D_cox;Covari=Z_cox;weight="fi"
  #b=beta_hat_aft;Time=T_aft;Delta=D_aft;Covari=Z_aft;weight="fi"
  
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
W_t.z_omni(beta_hat_aft,T_aft,D_aft,Z_aft,"fi")

W_z_ftnform=function(b,Time,Delta,Covari,weight="1i"){
  #b=beta_hat_cox;Time=T_cox;Delta=D_cox;Covari=Z_cox;weight="1i"
  #b=beta_hat_aft;Time=T_aft;Delta=D_aft;Covari=Z_aft;weight="1i"
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
  
  Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
  #Mhat_i_inf
  
  w_ij_z.Mhat_i_inf=list(NA)
  for(j in 1:p){
    w_ij_z.Mhat_i_inf[[j]]=mapply("*",w_ij_z[[j]],Mhat_i_inf,SIMPLIFY = FALSE)
  }
  #w_ij_z.Mhat_i_inf
  # z by 1 matrix
  
  W_j_inf.z=lapply(w_ij_z.Mhat_i_inf,function(x){Reduce('+',x)})
  #W_j_inf.z
  
  result=list(Time,Delta,Covari,e_i_beta,W_j_inf.z)
  names(result)=c("Time","Delta","Covari","Resid","W_j_inf.z")
  
  return(result)
}
#W_z_ftnform()
W_z_ftnform(beta_hat_aft,T_aft,D_aft,Z_aft,"1i")

W_z_linkftn=function(b,Time,Delta,Covari,weight="1i"){
  #b=beta_hat_cox;Time=T_cox;Delta=D_cox;Covari=Z_cox;weight="1i"
  #b=beta_hat_aft;Time=T_aft;Delta=D_aft;Covari=Z_aft;weight="1i"
  
  # Covari is n by J matrix consited of the covariates
  Covari=matrix(Covari,nrow=n)
  
  n=length(Time) # the number of individuals
  p=length(b) # the number of parameters
  
  if(p==1){return(print("ERROR MESSAGE : the number  needs to be greater than one"))}
  
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
  
  Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
  #Mhat_i_inf
  
  w_ij_z.Mhat_i_inf=list(NA)
  for(j in 1:p){
    w_ij_z.Mhat_i_inf[[j]]=mapply("*",w_ij_z[[j]],Mhat_i_inf,SIMPLIFY = FALSE)
  }
  #w_ij_z.Mhat_i_inf
  # z by 1 matrix
  
  W_j_inf.z=lapply(w_ij_z.Mhat_i_inf,function(x){Reduce('+',x)})
  #W_j_inf.z
  
  result=list(Time,Delta,Covari,e_i_beta,W_j_inf.z)
  names(result)=c("Time","Delta","Covari","Resid","W_j_inf.z")
  
  return(result)
}
#W_z_linkftn()
W_z_linkftn(beta_hat_aft,T_aft,D_aft,Z_aft,"1i")

W_t.z=function(b,Time,Delta,Covari,weight,test){
  if(test=="omni"){
    return(W_t.z_omni(b,Time,Delta,Covari,weight))
  }
  if(test=="ftnform"){
    return(W_z_ftnform(b,Time,Delta,Covari,weight))
  }
  if(test=="linkftn"){
    return(W_z_linkftn(b,Time,Delta,Covari,weight))
  }
  if(test=="aft"){
    return(print("NOT YET..."))
  }
}
#W_t.z()
W_t.z(beta_hat_aft,T_aft,D_aft,Z_aft,"1i","omni")

#-------------------------------------------------------------
#-------------------------SAMPLE PATH-------------------------
#-------------------------------------------------------------
What_t=function(b,std,Time,Delta,Covari,weight,test,tol){
  #b=beta_hat_cox;std=std_hat_cox;Time=T_cox;Delta=D_cox;Covari=Z_cox;weight=given_weight;test=given_test;tol=given_tol;
  #b=beta_hat_aft;std=std_hat_aft;Time=T_aft;Delta=D_aft;Covari=Z_aft;weight=given_weight;test=given_test;tol=given_tol;
  Covari=as.matrix(Covari,nrow=n)
  # Covari is n by J matrix consited of the covariates
  
  n=length(Time) # the number of individuals
  j=ncol(Covari) # the number of parameters
  
  e_i_beta=log(Time)-Covari*b
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=Covari[order_resid]
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  if (test=="omni"){
    if (weight=="a"){w_i.z=Covari*(Covari<=median(Covari))}
    if (weight=="b"){w_i.z=Covari}
    if (weight=="c"){w_i.z=1*(Covari<=median(Covari))}
    if (weight=="d"){w_i.z=1}
  }
  
  if (test=="ftn.form"){
    Covari_w=Covari[order(Covari)]
    
    w_i.z=list(NA)
    for(j in 1:n){
      w_i.z[[j]]=(Covari_w>=Covari_w[j])*1
    }
    #w_i.z
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
  
  Mhat_i_s_t=mapply("-",N_i_t,lapply(lapply(
    Y_i_t,"*",dAhat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_s_t
  
  dN_i_t=lapply(N_i_t,function(x){diff(c(0,x))})
  #dN_i_t
  
  Q_t=S_0_t/n # Gehan's weight
  #Q_t Q_t=1 #
  
  U_t=Reduce('+',lapply(mapply("*",dN_i_t,lapply(
    lapply(Covari,'-',E_t),"*",Q_t), SIMPLIFY = FALSE),cumsum))
  #U_t
  
  S_w_s_t=Reduce('+',mapply(
    "*", Y_i_t,w_i.z, SIMPLIFY = FALSE))
  #S_w_s_t
  
  E_w_s_t=S_w_s_t/S_0_t
  #E_w_s_t
  
  dMhat_i_s_t=lapply(Mhat_i_s_t,function(x){diff(c(0,x))})
  #dMhat_i_s_t
  
  #-----------------------------------------------------------
  #----------------------kernel Smoothing---------------------
  #-----------------------------------------------------------
  
  #-----------------------------g0----------------------------
  Ghat_0_t=1-exp(-Ahat_0_t)
  #Ghat_0_t
  
  dGhat_0_t=diff(c(0,Ghat_0_t))
  #dGhat_0_t
  
  ghat_0_t=(ksmooth(e_i_beta,dGhat_0_t,"normal",
                    bandwidth = 1.06*sd(dGhat_0_t)*n^(-0.2),x.points=e_i_beta)$y)
  #ghat_0_t
  
  fhat_Y_t=Reduce('+',lapply(w_i.z*Covari,'*',ghat_0_t*Time))/n
  #fhat_Y_t
  
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
  
  fhat_N_t=Reduce('+',lapply(Delta*w_i.z*Covari,'*',fhat_0_t*Time))/n
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
    
    N_i_t_U=list(NA)
    for(j in 1:n){
      N_i_t_U[[j]]=(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]
    }
    #N_i_t_U
    
    Y_i_t_U=list(NA)
    for(j in 1:n){
      Y_i_t_U[[j]]=(e_i_beta_U<=e_i_beta_U[j])*1
    }
    #Y_i_t_U
    
    N_d_t_U=Reduce('+',N_i_t_U)
    #N_d_t_U
    
    S_0_t_U=Reduce('+',Y_i_t_U)
    #S_0_t_U
    
    S_1_t_U=Reduce('+',mapply(
      "*",Y_i_t_U, Covari_U, SIMPLIFY = FALSE))
    #S_1_t_U
    
    E_t_U=S_1_t_U/S_0_t_U
    #E_t
    
    J_t_U=(S_0_t_U>0)*1
    #J_t_U
    
    dN_d_t_U=diff(c(0,N_d_t_U))
    #dN_d_t_U
    
    Ahat_0_t_U=cumsum((J_t_U/S_0_t_U)*dN_d_t_U)
    #Ahat_0_t_U
    
    dAhat_0_t_U=diff(c(0,Ahat_0_t_U))
    #dAhat_0_t _U
    
    Mhat_i_s_t=mapply("-", N_i_t_U,lapply(
      lapply(Y_i_t_U,"*",dAhat_0_t_U),cumsum), SIMPLIFY = FALSE)
    #Mhat_i_s_t _U
    
    dN_i_t_U=lapply(N_i_t_U,function(x){diff(c(0,x))})
    #dN_i_t_U
    
    Q_t_U=S_0_t_U/n # Gehan's weight
    #Q_t_U Q_t_U=1 #
    
    U_t_U=Reduce('+',lapply(mapply("*",dN_i_t_U,lapply(
      lapply(Covari_U,'-',E_t_U),"*",Q_t_U), SIMPLIFY = FALSE),cumsum))
    #U_t_U
    
    #U_t_U_order=U_t_U[order(Time_U)]
    #U_t_U_order
    
    U_inf_U=U_t_U[n]
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  tolerance=tol+1 #initial value
  
  while (tolerance>tol){
    
    G_i=rnorm(n)
    #G_i
    
    U_w_G_t=Reduce('+',lapply(mapply("*",mapply(
      "*",dMhat_i_s_t,lapply(lapply(w_i.z,'-',E_w_s_t),"*",Q_t)
      , SIMPLIFY = FALSE),G_i,SIMPLIFY = FALSE),cumsum))
    #U_w_G_t
    
    U_G_t=Reduce('+',lapply(mapply("*",mapply(
      "*",dMhat_i_s_t,lapply(lapply(Covari,'-',E_t),"*",Q_t)
      , SIMPLIFY = FALSE),G_i,SIMPLIFY = FALSE),cumsum))
    #U_G_t
    
    #U_G_t_order=U_G_t[order(Time)]
    #U_G_t_order
    
    U_G_inf=U_G_t[n]
    #U_G_t.inf
    
    beta_hat_s_list=optimize(function(beta){abs(U_beta(beta_U=beta)-U_G_inf)},
                             c(b-5*std,b+5*std),
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
  
  S_1_t_s=Reduce('+',mapply("*", Y_i_t_s, Covari_s,
                                   SIMPLIFY = FALSE))
  #S_1_t_s
  
  E_t_s=S_1_t_s/S_0_t_s
  #E_t_s
  
  J_t_s=(S_0_t_s>0)*1
  #J_t_s
  
  dN_d_t_s=diff(c(0,N_d_t_s))
  #dN_d_t_s
  
  Ahat_0_t_s=cumsum((J_t_s/S_0_t_s)*dN_d_t_s)
  #Ahat_0_t_s
  
  dAhat_0_t_s=diff(c(0,Ahat_0_t_s))
  #dAhat_0_t_s
  
  order_Time=order(Time)
  #order_Time
  
  order_Time_s=order(Time_s)
  #order_Time_s
  
  #order_Covari=order(Covari)
  
  F.T.=(1/sqrt(n))*U_w_G_t
  S.T.=sqrt(n)*(fhat_N_t+cumsum(fhat_Y_t*dAhat_0_t))*(b-beta_hat_s)
  T.T.=(1/sqrt(n))*cumsum(S_w_s_t*diff(c(0,Ahat_0_t-Ahat_0_t_s)))
  
  What_t=F.T.-S.T.-T.T.
  #What_t
  
  #return(G_i)
  #return(beta_hat_s_list)
  #return(AA)
  #return(BB)
  #return(CC)
  #return(What_t)
  #return(list(What_t=What_t,beta_hat_s=beta_hat_s))
  
  if (test=="omni"){
    return(What_t[order_Time])
    }
  if (test=="ftn.form"){return(What_t)}
  
}
#What_t()

sample_path=function(path,b,std,Time,Delta,Covari,weight,test,tol){
  #path=path;b=beta_hat_cox;std=std_hat_cox;Time=T_cox;Delta=D_cox;Covari=Z_cox;weight=given_weight;test=given_test;tol=given_tol;
  #path=path;b=beta_hat_aft;std=std_hat_aft;Time=T_aft;Delta=D_aft;Covari=Z_aft;weight=given_weight;test=given_test;tol=given_tol;
  
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
#system.time(sample_path(path,beta_hat_aft,std_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,given_tol))
#--------------------------CENSORING--------------------------
result_aft=sample_path(path,beta_hat_aft,std_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,given_tol)
result_cox=sample_path(path,beta_hat_cox,std_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,given_tol)

result_aft$p_value
result_aft$std_p_value

result_cox$p_value
result_cox$std_p_value

# PLOT : W_aft vs What_aft
Figure1_W_aft=plotting(result_aft,0,50);Figure1_W_aft

# PLOT : std.W_aft vs std.What_aft
Figure1_std.W_aft=plotting(result_aft,1,50);Figure1_std.W_aft

# PLOT : W_cox vs What_cox
Figure1_W_cox=plotting(result_cox,0,50);Figure1_W_cox

# PLOT : W_cox vs What_cox
Figure1_std.W_cox=plotting(result_cox,1,50);Figure1_std.W_cox



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






bb=seq(-400,length.out = 801)/100
plot(dgev(bb, loc=0, scale=1, shape=0, log = FALSE),type = "l",col="red")
par(new=TRUE)
plot(dgev(bb, loc=0, scale=1, shape=0.5, log = FALSE),type = "l",col="blue")
par(new=TRUE)
plot(dgev(bb, loc=0, scale=1, shape=-0.5, log = FALSE),type = "l",col="green")

#------------------------DATA GENERATION----------------------
n=200
beta_0=1
Z=rnorm(n,1,1)

#---------------------WEIBULL DISTRIBUTION--------------------
# censoring ratio : about 20%
mu.T_wb=0;sigma.T_wb=0.5;k.T_wb=1;
mu.C_wb=2;sigma.C_wb=0.5;k.C_wb=1; 

# data generation from weibull distribution i.e k=1
log.T_wb=beta_0%*%Z+rgev(n,loc=mu.T_wb,scale=sigma.T_wb,shape=k.T_wb) # true failure time
log.C_wb=beta_0%*%Z+rgev(n,loc=mu.C_wb,scale=sigma.C_wb,shape=k.C_wb) # censoring time
log.X_wb=log.C_wb*(log.T_wb>log.C_wb)+log.T_wb*(log.T_wb<=log.C_wb) # observed failure time
D_wb=0*(log.T_wb>log.C_wb)+1*(log.T_wb<=log.C_wb) # delta 0:censored & 1:observed

# ordering
D_wb=D_wb[order(log.X_wb)]
Z_wb=Z[order(log.X_wb)]
log.X_wb=log.X_wb[order(log.X_wb)]

length(which(D_wb==0))

#------------WEIBULL DISTRIBUTION(FUNCTIONAL FORM)------------
Z_f=matrix(c(Z,Z^2), nrow = n, ncol = 2)
gamma_0=0.1
beta=c(beta_0,gamma_0)

# censoring ratio : about 20%
mu.T_f=0;sigma.T_f=0.5;k.T_f=1;
mu.C_f=2;sigma.C_f=0.5;k.C_f=1; 

# data generation from weibull distribution i.e k=1
log.T_wb_f=Z_f%*%beta+rgev(n,mu.T_f,sigma.T_f,k.T_f) # true failure time
log.C_wb_f=Z_f%*%beta+rgev(n,mu.C_f,sigma.C_f,k.C_f) # censoring time
log.X_wb_f=log.C_wb_f*(log.T_wb_f>log.C_wb_f)+log.T_wb_f*(log.T_wb_f<=log.C_wb_f) # observed failure time
D_wb_f=0*(log.T_wb_f>log.C_wb_f)+1*(log.T_wb_f<=log.C_wb_f) # delta 0:censored & 1:observed

# ordering
D_wb_f=D_wb[order(log.X_wb_f)]
Z_wb_f=Z[order(log.X_wb_f)]
log.X_wb_f=log.X_wb[order(log.X_wb_f)]

length(which(D_wb_f==0))

#------------GNERALIZED EXTREME VALUE DISTRIBUTION------------
# censoring ratio : about 20%
mu.T_gev=0;sigma.T_gev=0.5;k.T_gev=0.1;
mu.C_gev=2;sigma.C_gev=0.5;k.C_gev=0.1; 

# data generation from weibull distribution i.e k=1
log.T_gev=beta_0%*%Z+rgev(n,mu.T_gev,sigma.T_gev,k.T_gev) # true failure time
log.C_gev=beta_0%*%Z+rgev(n,mu.C_gev,sigma.C_gev,k.C_gev) # censoring time
log.X_gev=log.C_gev*(log.T_gev>log.C_gev)+log.T_gev*(log.T_gev<=log.C_gev) # observed failure time
D_gev=0*(log.T_gev>log.C_gev)+1*(log.T_gev<=log.C_gev) # delta 0:censored & 1:observed

# ordering
D_gev=D_gev[order(log.X_gev)]
Z_gev=Z[order(log.X_gev)]
log.X_gev=log.X_gev[order(log.X_gev)]

length(which(D_gev==0))
log.X_gev

#-------------------------------------------------------------
#-------------Estimate Beta_hat_wb by using Aftgee------------
#-------------------------------------------------------------
aftsrr_beta_wb=aftsrr(Surv(exp(log.X_wb),D_wb)~Z_wb,method="nonsm")
beta_hat_wb=unlist(summary(aftsrr_beta_wb))$coefficients1;beta_hat_wb
std_hat_wb=unlist(summary(aftsrr_beta_wb))$coefficients2;std_hat_wb

#-------------------------------------------------------------
#------------Estimate Beta_hat_wb_f by using Aftgee-----------
#-------------------------------------------------------------
aftsrr_beta_wb_f=aftsrr(Surv(exp(log.X_wb_f),D_wb_f)~Z_wb_f,method="nonsm")
beta_hat_wb_f=unlist(summary(aftsrr_beta_wb_f))$coefficients1;beta_hat_wb_f
std_hat_wb_f=unlist(summary(aftsrr_beta_wb_f))$coefficients2;std_hat_wb_f

#-------------------------------------------------------------
#------------Estimate Beta_hat_gev by using Aftgee------------
#-------------------------------------------------------------
aftsrr_beta_gev=aftsrr(Surv(exp(log.X_gev),D_gev)~Z_gev,method="nonsm")
beta_hat_gev=unlist(summary(aftsrr_beta_gev))$coefficients1;beta_hat_gev
std_hat_gev=unlist(summary(aftsrr_beta_gev))$coefficients2;std_hat_gev







