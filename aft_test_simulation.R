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
#-------------------------OMNIBUS TEST------------------------
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

#-------------------------------------------------------------
#-----------------------FUNCTIONAL FORM-----------------------
#-------------------------------------------------------------
simulation_fform=function(simulation,n,path,alpha,tol){
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
    T_ln_aft=as.vector(exp(-beta_0*Z)*qlnorm(runif(n),5,1))
    C_ln_aft=as.vector(exp(-beta_0*Z)*qlnorm(runif(n),6.5,1))
    X_ln_aft=C_ln_aft*(T_ln_aft>C_ln_aft)+T_ln_aft*(T_ln_aft<=C_ln_aft)
    D_ln_aft=0*(T_ln_aft>C_ln_aft)+1*(T_ln_aft<=C_ln_aft)
    Z_ln_aft=Z
    
    T_ln_aft_f=as.vector(exp(-beta_0*Z-gamma_0*(Z^2))*qlnorm(runif(n),5,1))
    C_ln_aft_f=as.vector(exp(-beta_0*Z-gamma_0*(Z^2))*qlnorm(runif(n),6.5,1))
    X_ln_aft_f=C_ln_aft_f*(T_ln_aft_f>C_ln_aft_f)+T_ln_aft_f*(T_ln_aft_f<=C_ln_aft_f)
    D_ln_aft_f=0*(T_ln_aft_f>C_ln_aft_f)+1*(T_ln_aft_f<=C_ln_aft_f)
    Z_ln_aft_f=Z
    
    #------------Estimate Beta_hat_ln_aft_f by using Aftgee-----------
    aftsrr_beta_ln_aft=aftsrr(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,method="nonsm")
    beta_hat_ln_aft=-as.vector(aftsrr_beta_ln_aft$beta);beta_hat_ln_aft
    std_hat_ln_aft=diag(aftsrr_beta_ln_aft$covmat$ISMB);std_hat_ln_aft
    
    aftsrr_beta_ln_aft_f=aftsrr(Surv(X_ln_aft_f,D_ln_aft_f)~Z_ln_aft_f,method="nonsm")
    beta_hat_ln_aft_f=-as.vector(aftsrr_beta_ln_aft_f$beta);beta_hat_ln_aft_f
    std_hat_ln_aft_f=diag(aftsrr_beta_ln_aft_f$covmat$ISMB);std_hat_ln_aft_f
    
    # result_ln_aft
    result_ln_aft=sample_path_fform(path,beta_hat_ln_aft,std_hat_ln_aft,
                                      X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
    
    # result_ln_aft_f
    result_ln_aft_f=sample_path_fform(path,beta_hat_ln_aft_f,std_hat_ln_aft_f,
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
    
    p_mean=rbind(c(result_ln_aft$p_value,result_ln_aft$std.p_value),
                 c(result_ln_aft_f$p_value,result_ln_aft_f$std.p_value))
    colnames(p_mean)=c("W","std.W")
    rownames(p_mean)=c("p_ln_aft_mean","p_ln_aft_f_mean")
    #p_mean
    
    p_alpha=(p_mean>=alpha)*1
    colnames(p_alpha)=c("W","std.W")
    rownames(p_alpha)=c("p_ln_aft_alpha","p_ln_aft_f_alpha")
    #p_alpha
    
    p_value=list(p_mean,p_alpha)
    #p_value
    
    #result[[k]]=list(result_ln_aft,result_ln_cox,p_value)
    result[[k]]=list(p_value)
  }
  return(result)
}
#simulation_fform

prob.table_fform=function(simul_result){
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
#prob.table_fform

date()
simulation_result_fform1=simulation_fform(simulation,n,path,alpha,given_tol)
prob.table_fform(simulation_result_fform1)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_fform1_n500p150sim100")
date()
simulation_result_fform2=simulation_fform(simulation,n,path,alpha,given_tol)
prob.table_fform(simulation_result_fform2)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_fform2_n500p150sim100")
date()
simulation_result_fform3=simulation_fform(simulation,n,path,alpha,given_tol)
prob.table_fform(simulation_result_fform3)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_fform3_n500p150sim100")
date()
simulation_result_fform4=simulation_fform(simulation,n,path,alpha,given_tol)
prob.table_fform(simulation_result_fform4)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_fform4_n500p150sim100")
date()
simulation_result_fform5=simulation_fform(simulation,n,path,alpha,given_tol)
prob.table_fform(simulation_result_fform5)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_fform5_n500p150sim100")
date()
simulation_result_fform=c(simulation_result_fform1,simulation_result_fform2,
                          simulation_result_fform3,simulation_result_fform4,
                          simulation_result_fform5)
prob.table_fform(simulation_result_fform)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_fform_n500p150sim500")
date()

#-------------------------------------------------------------
#------------------------LINK FUNCTION------------------------
#-------------------------------------------------------------
simulation_linkf=function(simulation,n,path,alpha,tol){
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
    Z1=matrix(abs(rnorm(n,3,1)),nrow=n)
    Z2=matrix(rnorm(n,1,1),nrow=n)
    
    #-------------------LOG NORMAL DISTRIBUTION-------------------
    T_ln_aft=as.vector(exp(-beta_0*Z1-gamma_0*Z2)*qlnorm(runif(n),5,1))
    C_ln_aft=as.vector(exp(-beta_0*Z1-gamma_0*Z2)*qlnorm(runif(n),6.5,1))
    X_ln_aft=C_ln_aft*(T_ln_aft>C_ln_aft)+T_ln_aft*(T_ln_aft<=C_ln_aft)
    D_ln_aft=0*(T_ln_aft>C_ln_aft)+1*(T_ln_aft<=C_ln_aft)
    Z1_ln_aft=Z1
    Z2_ln_aft=Z2
    Z_ln_aft=cbind(Z1_ln_aft,Z2_ln_aft)
    
    T_ln_aft_l=as.vector(exp(-beta_0*sqrt(Z1)-gamma_0*(Z2^2))*qlnorm(runif(n),5,1))
    C_ln_aft_l=as.vector(exp(-beta_0*sqrt(Z1)-gamma_0*(Z2^2))*qlnorm(runif(n),6.5,1))
    X_ln_aft_l=C_ln_aft_l*(T_ln_aft_l>C_ln_aft_l)+T_ln_aft_l*(T_ln_aft_l<=C_ln_aft_l)
    D_ln_aft_l=0*(T_ln_aft_l>C_ln_aft_l)+1*(T_ln_aft_l<=C_ln_aft_l)
    Z1_ln_aft_l=Z1
    Z2_ln_aft_l=Z2
    Z_ln_aft_l=cbind(Z1_ln_aft_l,Z2_ln_aft_l)
    
    #------------Estimate Beta_hat_ln_aft by using Aftgee-----------
    aftsrr_beta_ln_aft=aftsrr(Surv(X_ln_aft,D_ln_aft)~Z1_ln_aft+Z2_ln_aft,method="nonsm")
    beta_hat_ln_aft=-as.vector(aftsrr_beta_ln_aft$beta);beta_hat_ln_aft
    std_hat_ln_aft=diag(aftsrr_beta_ln_aft$covmat$ISMB);std_hat_ln_aft
    
    aftsrr_beta_ln_aft_l=aftsrr(Surv(X_ln_aft_l,D_ln_aft_l)~Z1_ln_aft_l+Z2_ln_aft_l,method="nonsm")
    beta_hat_ln_aft_l=-as.vector(aftsrr_beta_ln_aft_l$beta);beta_hat_ln_aft_l
    std_hat_ln_aft_l=diag(aftsrr_beta_ln_aft_l$covmat$ISMB);std_hat_ln_aft_l
    
    # result_ln_aft
    result_ln_aft=sample_path_linkf(path,beta_hat_ln_aft,std_hat_ln_aft,
                                    X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
    
    result_ln_aft_l=sample_path_linkf(path,beta_hat_ln_aft_l,std_hat_ln_aft_l,
                                      X_ln_aft_l,D_ln_aft_l,Z_ln_aft_l,given_tol)
    
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
                 c(result_ln_aft_l$p_value,result_ln_aft_l$std.p_value))
    colnames(p_mean)=c("W","std.W")
    rownames(p_mean)=c("p_ln_aft_mean","p_ln_aft_l_mean")
    #p_mean
    
    p_alpha=(p_mean>=alpha)*1
    colnames(p_alpha)=c("W","std.W")
    rownames(p_alpha)=c("p_ln_aft_alpha","p_ln_aft_l_alpha")
    #p_alpha
    
    p_value=list(p_mean,p_alpha)
    #p_value
    
    #result[[k]]=list(result_ln_aft,result_ln_cox,p_value)
    result[[k]]=list(p_value)
  }
  return(result)
}
#simulation_linkf

prob.table_linkf=function(simul_result){
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
#prob.table_linkf

date()
simulation_result_linkf1=simulation_linkf(simulation,n,path,alpha,given_tol)
prob.table_linkf(simulation_result_linkf1)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_linkf1_n500p150sim100")
date()
simulation_result_linkf2=simulation_linkf(simulation,n,path,alpha,given_tol)
prob.table_linkf(simulation_result_linkf2)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_linkf2_n500p150sim100")
date()
simulation_result_linkf3=simulation_linkf(simulation,n,path,alpha,given_tol)
prob.table_linkf(simulation_result_linkf3)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_linkf3_n500p150sim100")
date()
simulation_result_linkf4=simulation_linkf(simulation,n,path,alpha,given_tol)
prob.table_linkf(simulation_result_linkf4)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_linkf4_n500p150sim100")
date()
simulation_result_linkf5=simulation_linkf(simulation,n,path,alpha,given_tol)
prob.table_linkf(simulation_result_linkf5)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_linkf5_n500p150sim100")
date()
simulation_result_linkf=c(simulation_result_linkf1,simulation_result_linkf2,
                          simulation_result_linkf3,simulation_result_linkf4,
                          simulation_result_linkf5)
prob.table_linkf(simulation_result_linkf)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_linkf_n500p150sim500")
date()



