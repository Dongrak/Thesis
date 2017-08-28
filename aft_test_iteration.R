iteration=1
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
  # 0.45이면 accetp <= W가 55번이나 튀어나간거다!
  # p_alpha는 acceptance rate을 구하는 것이다! 
  #-----------------------------------------------------------
  
  p_exact=rbind(c(result_aft$p_value,result_aft$std_p_value),
                c(result_cox$p_value,result_cox$std_p_value))
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

#iteration_result2 # iteration 150
#iteration_result1 # iteration 200

#iteration_result3=c(iteration_result1,iteration_result2)

#prob.table(iteration_result3)
