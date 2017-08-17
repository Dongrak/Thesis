alpha=0.9
beta=0.55

iteration=1

result=list()
result_aft=list()
result_cox=list()
result_p=list()

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
  
  # dataset_What()
  dataset_What_aft=simulation_What(sim,beta_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,given_tol)
  dataset_What_cox=simulation_What(sim,beta_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,given_tol)
  
  # dataset_W()
  dataset_W_aft=simulation_W(beta_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,dataset_What_aft)
  dataset_W_cox=simulation_W(beta_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,dataset_What_cox)
  
  # dataset_What50_aft=dataset_What_aft[1:(n*50),]
  # dataset_What50_cox=dataset_What_cox[1:(n*50),]
  
  # PLOT : W_aft vs What_aft
  # Figure1_W_aft=
  # ggplot()+
  # geom_line(data=dataset_What_aft,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
  # geom_line(data=dataset_W_aft,aes(x=t_i,y=W),colour="tomato")
  # Figure1_W_aft
  
  # PLOT : W_cox vs What_cox
  # Figure1_W_cox=
  # ggplot()+
  # geom_line(data=dataset_What_cox,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
  # geom_line(data=dataset_W_cox,aes(x=t_i,y=W),colour="tomato")
  # Figure1_W_cox
  
  # PLOT : std.W_aft vs std.What_aft
  # Figure1_std.W_aft=
  # ggplot()+
  # geom_line(data=dataset_What_aft,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
  # geom_line(data=dataset_W_aft,aes(x=t_i,y=std.W),colour="tomato")
  # Figure1_std.W_aft
  
  # PLOT : W_cox vs What_cox
  # Figure1_std.W_cox=
  # ggplot()+
  # geom_line(data=dataset_What_cox,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
  # geom_line(data=dataset_W_cox,aes(x=t_i,y=std.W),colour="tomato")
  # Figure1_std.W_cox
  
  kol_typ_test_aft=kolmogorov(dataset_W_aft,dataset_What_aft)
  kol_typ_test_cox=kolmogorov(dataset_W_cox,dataset_What_cox)
  
  p_aft=kol_typ_test_aft[3,]
  p_cox=kol_typ_test_cox[3,]
  
  P_result_exact=rbind(p_aft,p_cox)
  
  P_result_alpha_beta=rbind((p_aft>=alpha)*1,(p_cox>beta)*1)
  rownames(P_result_alpha_beta)=c("p_aft_alpha","p_cox_beta")
  
  P_result_binary=rbind((p_aft==1)*1,(p_cox>0)*1)
  rownames(P_result_binary)=c("p_aft_1","p_cox_0")
  
  result_aft[[k]]=list(#Figure1_W_aft.aft,Figure1_std.W_aft.aft,
    dataset_W_aft,dataset_What_aft,kol_typ_test_aft)
  
  result_cox[[k]]=list(#Figure1_W_cox.cox,Figure1_std.W_cox.cox,
    dataset_W_cox,dataset_What_cox,kol_typ_test_cox)
  
  result_p[[k]]=list(P_result_exact,P_result_alpha_beta,P_result_binary)
  
  result[[k]]=list(result_aft[[k]],result_cox[[k]],result_p[[k]])
}

result_n_exact=0
result_n_alpha_beta=0
result_n_binary=0

for(k in 1:iteration){
  result_n_exact=result_n_exact+result_p[[k]][[1]]
  result_n_alpha_beta=result_n_alpha_beta+result_p[[k]][[2]]
  result_n_binary=result_n_binary+result_p[[k]][[3]]
}

P_result_exact=result_n_exact/iteration;P_result_exact
P_result_alpha_beta=result_n_alpha_beta/iteration;P_result_alpha_beta
P_result_binary=result_n_binary/iteration;P_result_binary

for(k in 1:iteration){
  print(result[[k]][[3]][[1]])
  #print(result[[k]][[3]][[2]])
  #print(result[[k]][[3]][[3]])
}
