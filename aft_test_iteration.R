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

iteration=200
n=200
path=200
alpha=0.05

given_tol=0.1

#-------------------------------------------------------------
#-------------------------OMNIBUS TEST------------------------
#-------------------------------------------------------------
iteration_function_omni=function(iteration,n,path,alpha,tol){
  #iteration=iteration;n=n;path=path;alpha=alpha;tol=given_tol;
  
  result=list(NA)
  
  for(k in 1:iteration){
    if(k%%1==0) {
      cat("Iteration",k,"\n")
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
    beta_hat_ln_cox=-as.vector(coxsrr_beta_ln_cox$beta);beta_hat_ln_cox
    std_hat_ln_cox=diag(coxsrr_beta_ln_cox$covmat$ISMB);std_hat_ln_cox
    
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
    
    result[[k]]=list(result_ln_aft,result_ln_cox,p_value)
  }
  return(result)
}
#iteration_function_omni

prob.table_omni=function(iter_result){
  iter=length(iter_result)
  
  p_mean_set=list(NA)
  p_alpha_set=list(NA)
  
  for(k in 1:iter){
    p_mean_set[[k]]=iter_result[[k]][[3]][[1]]
    p_alpha_set[[k]]=iter_result[[k]][[3]][[2]]
  }
  
  p_mean=Reduce("+",p_mean_set)/iter
  p_alpha=Reduce("+",p_alpha_set)/iter
  
  return(list(p_mean,p_alpha))
}
#prob.table_omni

date()
iteration_result_omni1=iteration_function_omni(iteration,n,path,alpha,given_tol)
prob.table_omni(iteration_result_omni1)
date()
iteration_result_omni2=iteration_function_omni(iteration,n,path,alpha,given_tol)
prob.table_omni(iteration_result_omni2)
date()
iteration_result_omni3=iteration_function_omni(iteration,n,path,alpha,given_tol)
prob.table_omni(iteration_result_omni3)
date()
iteration_result_omni4=iteration_function_omni(iteration,n,path,alpha,given_tol)
prob.table_omni(iteration_result_omni4)
date()
iteration_result_omni5=iteration_function_omni(iteration,n,path,alpha,given_tol)
prob.table_omni(iteration_result_omni5)
date()
iteration_result_omni=c(iteration_result_omni1,iteration_result_omni2,
                        iteration_result_omni3,iteration_result_omni4,
                        iteration_result_omni5)
prob.table_omni(iteration_result_omni)
date()

#-------------------------------------------------------------
#-----------------------FUNCTIONAL FORM-----------------------
#-------------------------------------------------------------
iteration_function_ftnform=function(iteration,n,path,alpha,tol){
  #iteration=iteration;n=n;path=path;alpha=alpha;tol=given_tol;
  
  result=list(NA)
  
  for(k in 1:iteration){
    if(k%%1==0) {
      cat("Iteration",k,"\n")
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
#iteration_function_ftnform

prob.table_ftnform=function(iter_result){
  iter=length(iter_result)
  
  p_mean_set=list(NA)
  p_alpha_set=list(NA)
  
  for(k in 1:iter){
    p_mean_set[[k]]=iter_result[[k]][[2]][[1]]
    p_alpha_set[[k]]=iter_result[[k]][[2]][[2]]
  }
  
  p_mean=Reduce("+",p_mean_set)/iter
  p_alpha=Reduce("+",p_alpha_set)/iter
  
  return(list(p_mean,p_alpha))
}
#prob.table_ftnform

date()
iteration_result_ftnform1=iteration_function_ftnform(iteration,n,path,alpha,given_tol)
prob.table_ftnform(iteration_result_ftnform1)
date()
iteration_result_ftnform2=iteration_function_ftnform(iteration,n,path,alpha,given_tol)
prob.table_ftnform(iteration_result_ftnform2)
date()
iteration_result_ftnform3=iteration_function_ftnform(iteration,n,path,alpha,given_tol)
prob.table_ftnform(iteration_result_ftnform3)
date()
iteration_result_ftnform4=iteration_function_ftnform(iteration,n,path,alpha,given_tol)
prob.table_ftnform(iteration_result_ftnform4)
date()
iteration_result_ftnform5=iteration_function_ftnform(iteration,n,path,alpha,given_tol)
prob.table_ftnform(iteration_result_ftnform5)
date()
iteration_result_ftnform=c(iteration_result_ftnform1,iteration_result_ftnform2,
                           iteration_result_ftnform3,iteration_result_ftnform4,
                           iteration_result_ftnform5)
prob.table_ftnform(iteration_result_ftnform)
date()

#-------------------------------------------------------------
#------------------------LINK FUNCTION------------------------
#-------------------------------------------------------------
iteration_function_linkftn=function(iteration,n,path,alpha,tol){
  #iteration=iteration;n=n;path=path;alpha=alpha;tol=given_tol;
  
  result=list(NA)
  
  for(k in 1:iteration){
    if(k%%1==0) {
      cat("Iteration",k,"\n")
    }
    # -------------------------------------------------------------
    # ------------------------DATA GENERATE------------------------
    # -------------------------------------------------------------
    # n=200
    beta_0=1
    gamma_0=0.1
    Z1=matrix(rnorm(n,3,1),nrow=n)
    Z2=matrix(rnorm(n,5,2),nrow=n)
    
    #-------------------LOG NORMAL DISTRIBUTION-------------------
    T_ln_aft=as.vector(exp(-beta_0*exp(Z1)-gamma_0*log(Z2))*qlnorm(runif(n),5,1))
    C_ln_aft=as.vector(exp(-beta_0*exp(Z1)-gamma_0*log(Z2))*qlnorm(runif(n),6.5,1))
    X_ln_aft=C_ln_aft*(T_ln_aft>C_ln_aft)+T_ln_aft*(T_ln_aft<=C_ln_aft)
    D_ln_aft=0*(T_ln_aft>C_ln_aft)+1*(T_ln_aft<=C_ln_aft)
    Z1_ln_aft=Z1
    Z2_ln_aft=Z2
    Z_ln_aft=cbind(Z1_ln_aft,Z2_ln_aft)
    
    #------------Estimate Beta_hat_ln_aft by using Aftgee-----------
    aftsrr_beta_ln_aft=aftsrr(Surv(X_ln_aft,D_ln_aft)~Z1_ln_aft+Z2_ln_aft,method="nonsm")
    beta_hat_ln_aft=-as.vector(aftsrr_beta_ln_aft$beta);beta_hat_ln_aft
    std_hat_ln_aft=diag(aftsrr_beta_ln_aft$covmat$ISMB);std_hat_ln_aft
    
    # result_ln_aft
    result_ln_aft=sample_path_linkftn(path,beta_hat_ln_aft,std_hat_ln_aft,
                                      X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
    
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
    
    p_mean=rbind(c(result_ln_aft$p_value,result_ln_aft$std.p_value))
    colnames(p_mean)=c("W","std.W")
    rownames(p_mean)=c("p_ln_aft_mean")
    #p_mean
    
    p_alpha=(p_mean>=alpha)*1
    colnames(p_alpha)=c("W","std.W")
    rownames(p_alpha)=c("p_ln_aft_alpha")
    #p_alpha
    
    p_value=list(p_mean,p_alpha)
    #p_value
    
    result[[k]]=list(result_ln_aft,p_value)
  }
  return(result)
}
#iteration_function_linkftn

prob.table_linkftn=function(iter_result){
  iter=length(iter_result)
  
  p_mean_set=list(NA)
  p_alpha_set=list(NA)
  
  for(k in 1:iter){
    p_mean_set[[k]]=iter_result[[k]][[2]][[1]]
    p_alpha_set[[k]]=iter_result[[k]][[2]][[2]]
  }
  
  p_mean=Reduce("+",p_mean_set)/iter
  p_alpha=Reduce("+",p_alpha_set)/iter
  
  return(list(p_mean,p_alpha))
}
#prob.table_linkftn

date()
iteration_result_linkftn1=iteration_function_linkftn(iteration,n,path,alpha,given_tol)
prob.table_linkftn(iteration_result_linkftn1)
date()
iteration_result_linkftn2=iteration_function_linkftn(iteration,n,path,alpha,given_tol)
prob.table_linkftn(iteration_result_linkftn2)
date()
iteration_result_linkftn3=iteration_function_linkftn(iteration,n,path,alpha,given_tol)
prob.table_linkftn(iteration_result_linkftn3)
date()
iteration_result_linkftn4=iteration_function_linkftn(iteration,n,path,alpha,given_tol)
prob.table_linkftn(iteration_result_linkftn4)
date()
iteration_result_linkftn5=iteration_function_linkftn(iteration,n,path,alpha,given_tol)
prob.table_linkftn(iteration_result_linkftn5)
date()
iteration_result_linkftn=c(iteration_result_linkftn1,iteration_result_linkftn2,
                           iteration_result_linkftn3,iteration_result_linkftn4,
                           iteration_result_linkftn5)
prob.table_linkftn(iteration_result_linkftn)
date()



