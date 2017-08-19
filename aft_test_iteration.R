iteration=1
n=200
path=200
alpha=0.95

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
  
  result=list(NA)
  
  for(k in 1:iteration){
    if(k%%1==0) {
      cat("Iteration",k,"\n")
    }
  # -------------------------------------------------------------
  # ------------------------DATA GENERATE------------------------
  # -------------------------------------------------------------
  id=c(1:n) # identification
  beta_0=1 # beta_0
  Z=rnorm(n,3,1) # covariate
  
  # -----------------------AFT DATA GENERATE---------------------
  T_s_aft=exp(-beta_0*Z)*qlnorm(runif(n),5,1) # lognormal baseline hazard
  C_aft=exp(-beta_0*Z)*qlnorm(runif(n),6.5,1) # censoring time for aft model
  T_aft=C_aft*(T_s_aft>C_aft)+T_s_aft*(T_s_aft<=C_aft) # observed time for aft model
  D_aft=0*(T_s_aft>C_aft)+1*(T_s_aft<=C_aft) # delta
  D_aft=D_aft[order(T_aft)]
  Z_aft=Z[order(T_aft)]
  T_aft=T_aft[order(T_aft)]
  
  # -----------------------COX DATA GENERATE---------------------
  T_s_cox=qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5,1,lower.tail = FALSE) # lognormal baseline hazard
  C_cox=qlnorm((1-runif(n))^(1/exp(beta_0*Z)),5.7,1,lower.tail = FALSE) # censoring time for cox model
  T_cox=C_cox*(T_s_cox > C_cox)+T_s_cox*(T_s_cox<=C_cox) # observed time for cox model
  D_cox=0*(T_s_cox > C_cox)+1*(T_s_cox <= C_cox) # delta
  D_cox=D_cox[order(T_cox)]
  Z_cox=Z[order(T_cox)]
  T_cox=T_cox[order(T_cox)]
  
  # -------------------------------------------------------------
  # ------------Estimate Beta_hat_aft by using Aftgee------------
  # -------------------------------------------------------------
  aftsrr_beta_aft=aftsrr(Surv(T_aft,D_aft)~Z_aft,method="nonsm")
  beta_hat_aft=-unlist(summary(aftsrr_beta_aft))$coefficients1;beta_hat_aft
  std_hat_aft=unlist(summary(aftsrr_beta_aft))$coefficients2;std_hat_aft
  
  # -------------------------------------------------------------
  # ------------Estimate Beta_hat_cox by using Aftgee------------
  # -------------------------------------------------------------
  aftsrr_beta_cox=aftsrr(Surv(T_cox,D_cox)~Z_cox,method="nonsm")
  beta_hat_cox=-unlist(summary(aftsrr_beta_cox))$coefficients1;beta_hat_cox
  std_hat_cox=unlist(summary(aftsrr_beta_cox))$coefficients2;std_hat_cox
  
  # result_aft
  result_aft=sample_path(path,beta_hat_aft,std_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,given_tol)
  
  # result_cox
  result_cox=sample_path(path,beta_hat_cox,std_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,given_tol)
  
  p_exact=rbind(c(result_aft$p_value,result_aft$std.p_value),
                c(result_cox$p_value,result_cox$std.p_value))
  colnames(p_exact)=c("W","std.W")
  rownames(p_exact)=c("p_aft_exact","p_cox_exact")
  #p_exact
  
  P_alpha=(p_exact>=alpha)*1
  colnames(P_alpha)=c("W","std.W")
  rownames(P_alpha)=c("p_aft_alpha","p_cox_alpha")
  #P_alpha
  
  p_value=list(p_exact,P_alpha)
  #p_value
  
  result[[k]]=list(result_aft,result_cox,p_value)
  }
  return(result)
}
iteration_result3=iteration_function(iteration,n,path,alpha,given_weight,given_test,given_tol)

#iteration_result # iteration 150
#iteration_result1 # iteration 200

iteration_result_2=c(iteration_result,iteration_result1)

length(iteration_result_2)

prob.table=function(iter_result){
  iter=length(iter_result)
  sum.prob.list=list(0)
  for(k in 1:iter){
    sum.prob.list=mapply("+",sum.prob.list,iter_result[[k]][[3]])
  }
  mean.prob.list=sum.prob.list/iter
  return(mean.prob.list)
}
prob.table(iteration_result_2)
