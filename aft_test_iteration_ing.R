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

given_weight="c"
#(weight=="a"){w_i=Covari*(Covari<=median(Covari))}
#(weight=="b"){w_i=Covari}  
#(weight=="c"){w_i=1*(Covari<=median(Covari))}
#(weight=="d"){w_i=1}

given_test="omni"
# test="link.ftn"
# test="ftn.form"
# test="....???"

given_tol=0.1
# given_tol

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
    gamma_0=0.1
    Z=matrix(rnorm(n,3,1),nrow=n)
    
    #-------------------LOG NORMAL DISTRIBUTION-------------------
    T_ln_aft=as.vector(exp(-beta_0*Z)*qlnorm(runif(n),5,1))
    C_ln_aft=as.vector(exp(-beta_0*Z)*qlnorm(runif(n),6.5,1))
    X_ln_aft=C_ln_aft*(T_ln_aft>C_ln_aft)+T_ln_aft*(T_ln_aft<=C_ln_aft)
    D_ln_aft=0*(T_ln_aft>C_ln_aft)+1*(T_ln_aft<=C_ln_aft)
    Z_ln_aft=Z
    
    T_ln_cox=as.vector(qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5,1,lower.tail = FALSE))
    C_ln_cox=as.vector(qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5.7,1,lower.tail = FALSE))
    X_ln_cox=C_ln_cox*(T_ln_cox > C_ln_cox)+T_ln_cox*(T_ln_cox<=C_ln_cox)
    D_ln_cox=0*(T_ln_cox > C_ln_cox)+1*(T_ln_cox <= C_ln_cox)
    Z_ln_cox=Z
    
    #------------Estimate Beta_hat_ln_aft_f by using Aftgee-----------
    aftsrr_beta_ln_aft=aftsrr(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,method="nonsm")
    beta_hat_ln_aft=-unlist(summary(aftsrr_beta_ln_aft))$coefficients1;beta_hat_ln_aft
    std_hat_ln_aft=unlist(summary(aftsrr_beta_ln_aft))$coefficients2;std_hat_ln_aft
    
    aftsrr_beta_ln_cox=aftsrr(Surv(X_ln_cox,D_ln_cox)~Z_ln_cox,method="nonsm")
    beta_hat_ln_cox=-unlist(summary(aftsrr_beta_ln_cox))$coefficients1;beta_hat_ln_cox
    std_hat_ln_cox=unlist(summary(aftsrr_beta_ln_cox))$coefficients2;std_hat_ln_cox
    
    
    # result_ln_aft
    result_ln_aft=sample_path_omni(path,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
    
    # result_ln_cox
    result_ln_cox=sample_path_omni(path,beta_hat_ln_cox,std_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
    
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
    
    p_exact=rbind(c(result_ln_aft$p_value,result_ln_aft$std.p_value),
                  c(result_ln_cox$p_value,result_ln_cox$std.p_value))
    colnames(p_exact)=c("W","std.W")
    rownames(p_exact)=c("p_ln_aft_exact","p_ln_cox_exact")
    #p_exact
    
    p_alpha=(p_exact>=alpha)*1
    colnames(p_alpha)=c("W","std.W")
    rownames(p_alpha)=c("p_ln_aft_alpha","p_ln_cox_alpha")
    #p_alpha
    
    p_value=list(p_exact,p_alpha)
    #p_value
    
    result[[k]]=list(result_ln_aft,result_ln_cox,p_value)
  }
  return(result)
}
date()
iteration_result_omni1=iteration_function_omni(iteration,n,path,alpha,given_tol)
date()

prob.table=function(iter_result){
  iter=length(iter_result)
  
  p_exact_set=list(NA)
  p_alpha_set=list(NA)
  
  for(k in 1:iter){
    p_exact_set[[k]]=iter_result[[k]][[3]][[1]]
    p_alpha_set[[k]]=iter_result[[k]][[3]][[2]]
  }
  
  p_exact=Reduce("+",p_exact_set)/iter
  p_alpha=Reduce("+",p_alpha_set)/iter
  
  return(list(p_exact,p_alpha))
}
prob.table(iteration_result_omni1)

#given_to1=100000
#iteration_result1 # iteration 100

#given_to1=0.1
#iteration_result2 # iteration 100
#iteration_result3 # iteration 100 
#iteration_result4 # iteration 400 
#iteration_result5 # iteration 400 

#iteration 1000
#iteration_result=c(iteration_result2,iteration_result3,iteration_result4,iteration_result5)
#prob.table(iteration_result)
#iteration_result[[1]]