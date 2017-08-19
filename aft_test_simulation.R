n_vector=c(100,200,500,1000)

#n_vector=c(10,20,50,100)
#sim_vector=c(10,15,20,25)

for(k in 1:length(n_vector)){
  cat("n_vector : ",k,"\n")
}

table_function=function(n_vector){
  result_aft=list()
  result_cox=list()
  
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
  #------------Estimate Beta_hat_aft by using Aftgee------------
  #-------------------------------------------------------------
  aftsrr_beta_aft=aftsrr(Surv(T_aft,D_aft)~Z_aft,method="nonsm")
  beta_hat_aft=-unlist(summary(aftsrr_beta_aft))$coefficients1;beta_hat_aft
  std_hat_aft=unlist(summary(aftsrr_beta_aft))$coefficients2;std_hat_aft
  
  #-------------------------------------------------------------
  #------------Estimate Beta_hat_cox by using Aftgee------------
  #-------------------------------------------------------------
  aftsrr_beta_cox=aftsrr(Surv(T_cox,D_cox)~Z_cox,method="nonsm")
  beta_hat_cox=-unlist(summary(aftsrr_beta_cox))$coefficients1;beta_hat_cox
  std_hat_cox=unlist(summary(aftsrr_beta_cox))$coefficients2;std_hat_cox
  
  #-------------------------------------------------------------
  #-------------------------------------------------------------
  #-------------------------------------------------------------
  
  # result_aft
  result_aft=smaple_path(path,beta_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,given_tol)
  
  # result_cox
  result_cox=smaple_path(path,beta_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,given_tol)
  
  result_p=rbind(c(result_aft$p_value,result_aft$std.p_value),
                 c(result_cox$p_value,result_cox$std.p_value))
  
  
}


#n_vector=c(100,200,500,1000,2000)
#sim_vector=c(200,500,1000,2000)

table_fuction=function(result){
hh=matrix(nrow=length(n_vector),ncol=length(sim_vector))
rownames(hh)=c("n=100","n=200","n=500","n=1000")
colnames(hh)=c("sim=100","sim=150","sim=200","sim=250")
for(k in 1:length(n_vector)){
  for(j in 1:length(sim_vector)){
    n=n_vector[k]
    sim=sim_vector[j]
    
    aa=result[[k]]
    i=4+j
    bb=aa[[i]]
    cc=sprintf("%1.3f",bb[3,1])
    dd=sprintf("%1.3f",bb[3,2])
    ee=paste(cc,dd)
    hh[k,j]=ee
    
  }
}
return(hh)
}

table_aft=table_fuction(result_aft);table_aft
table_cox=table_fuction(result_cox);table_cox

#result_aft[[5]][2]
#result_cox[[5]][2]

