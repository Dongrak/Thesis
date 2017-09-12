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

library(ggplot2)
library(survival)
library(aftgee)

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

iteration_function=function(iteration,n,path,alpha,weight,test,tol){
  #iteration=iteration;n=n;path=path;alpha=alpha;weight=given_weight;test=given_test;tol=given_tol;
  
  result=list(NA)
  
  for(k in 1:iteration){
    if(k%%1==0) {
      cat("Iteration",k,"\n")
    }
  # -------------------------------------------------------------
  # ------------------------DATA GENERATE------------------------
  # -------------------------------------------------------------
  n=200
  beta_0=1
  gamma_0=0.1
  Z=matrix(rnorm(n,3,1),nrow=n)
  
  #---------------------WEIBULL DISTRIBUTION--------------------
  # T ~ Generalized Gamma(alpha=1,beta,sigma) i.e. Weibull ~ (beta,sigma)
  # alpha : shape // beta : scale // sigma : rate parameter of generaized gamma
  alpha_wb_T=1;alpha_wb_C=1     # alpha must be one when weibull
  beta_wb_T=5;beta_wb_C=6.5       #
  sigma_wb_T=5;sigma_wb_C=5     #
  
  V_wb_T=log(beta_wb_T)+log(qgamma(runif(n),alpha_wb_T,1))/sigma_wb_T
  V_wb_C=log(beta_wb_C)+log(qgamma(runif(n),alpha_wb_C,1))/sigma_wb_C
  T_wb=exp(-Z%*%beta_0+V_wb_T)
  C_wb=exp(-Z%*%beta_0+V_wb_C)
  X_wb=C_wb*(T_wb>C_wb)+T_wb*(T_wb<=C_wb) # observed failure time
  D_wb=0*(T_wb>C_wb)+1*(T_wb<=C_wb)       # delta 0:censored & 1:observed
  
  D_wb=D_wb[order(X_wb)]
  Z_wb=Z[order(X_wb)]
  X_wb=X_wb[order(X_wb)]
  #length(which(D_wb==0))/n
  
  #---------------GNERALIZED GUMBELL DISTRIBUTION---------------
  # T ~ Generalized Gamma(alpha,beta,sigma) if. V ~ Generalized Gumbell(beta,sigma)
  # alpha : shape // beta : scale // sigma : rate parameter of generaized gamma
  alpha_gg_T=100;alpha_gg_C=100    # alpha must be greater than one whne not weibull
  beta_gg_T=5;beta_gg_C=5.1       # gamma distribution
  sigma_gg_T=5;sigma_gg_C=5        #
  
  V_gg_T=log(beta_gg_T)+log(qgamma(runif(n),alpha_gg_T,1))/sigma_gg_T
  V_gg_C=log(beta_gg_C)+log(qgamma(runif(n),alpha_gg_C,1))/sigma_gg_C
  T_gg=exp(-Z%*%beta_0+V_gg_T)
  C_gg=exp(-Z%*%beta_0+V_gg_C)
  X_gg=C_gg*(T_gg>C_gg)+T_gg*(T_gg<=C_gg) # observed failure time
  D_gg=0*(T_gg>C_gg)+1*(T_gg<=C_gg)       # delta 0:censored & 1:observed
  
  D_gg=D_gg[order(X_gg)]
  Z_gg=Z[order(X_gg)]
  X_gg=X_gg[order(X_gg)]
  #length(which(D_gg==0))/n
  
  #------------WEIBULL DISTRIBUTION(FUNCTIONAL FORM)------------
  Z_f=matrix(c(Z,Z^2), nrow = n, ncol = 2)
  beta_f=c(beta_0,gamma_0)
  
  alpha_wb_f_T=1;alpha_wb_f_C=1     # alpha must be one when weibull
  beta_wb_f_T=10;beta_wb_f_C=12     #
  sigma_wb_f_T=10;sigma_wb_f_C=12   #
  
  V_wb_f_T=log(beta_wb_f_T)+log(qgamma(runif(n),alpha_wb_f_T,1))/sigma_wb_f_T
  V_wb_f_C=log(beta_wb_f_C)+log(qgamma(runif(n),alpha_wb_f_C,1))/sigma_wb_f_C
  T_wb_f=exp(-Z_f%*%beta_f+V_wb_f_T)
  C_wb_f=exp(-Z_f%*%beta_f+V_wb_f_C)
  X_wb_f=C_wb_f*(T_wb_f>C_wb_f)+T_wb_f*(T_wb_f<=C_wb_f) # observed failure time
  D_wb_f=0*(T_wb_f>C_wb_f)+1*(T_wb_f<=C_wb_f)       # delta 0:censored & 1:observed
  
  D_wb_f=D_wb_f[order(X_wb_f)]
  Z_wb_f=Z[order(X_wb_f)]
  X_wb_f=X_wb_f[order(X_wb_f)]
  #length(which(D_wb_f==0))/n
  
  #-------------------------------------------------------------
  #-------------Estimate Beta_hat_wb by using Aftgee------------
  #-------------------------------------------------------------
  aftsrr_beta_wb=aftsrr(Surv(X_wb,D_wb)~Z_wb,method="nonsm")
  beta_hat_wb=-unlist(summary(aftsrr_beta_wb))$coefficients1;beta_hat_wb
  std_hat_wb=unlist(summary(aftsrr_beta_wb))$coefficients2;std_hat_wb
  
  #-------------------------------------------------------------
  #------------Estimate Beta_hat_gev by using Aftgee------------
  #-------------------------------------------------------------
  aftsrr_beta_gg=aftsrr(Surv(X_gg,D_gg)~Z_gg,method="nonsm")
  beta_hat_gg=-unlist(summary(aftsrr_beta_gg))$coefficients1;beta_hat_gg
  std_hat_gg=unlist(summary(aftsrr_beta_gg))$coefficients2;std_hat_gg
  
  # result_wb
  result_wb=sample_path(path,beta_hat_wb,std_hat_wb,X_wb,D_wb,Z_wb,given_weight,given_test,given_tol)
  
  # result_gg
  result_gg=sample_path(path,beta_hat_gg,std_hat_gg,X_gg,D_gg,Z_gg,given_weight,given_test,given_tol)
  
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
  
  p_exact=rbind(c(result_wb$p_value,result_wb$std_p_value),
                c(result_gg$p_value,result_gg$std_p_value))
  colnames(p_exact)=c("W","std.W")
  rownames(p_exact)=c("p_wb_exact","p_gg_exact")
  #p_exact
  
  P_alpha=(p_exact>=alpha)*1
  colnames(P_alpha)=c("W","std.W")
  rownames(P_alpha)=c("p_wb_alpha","p_gg_alpha")
  #P_alpha
  
  p_value=list(p_exact,P_alpha)
  #p_value
  
  result[[k]]=list(result_wb,result_gg,p_value)
  }
  return(result)
}
date()
iteration_result1=iteration_function(iteration,n,path,alpha,given_weight,given_test,given_tol)
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
prob.table(iteration_result1)

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