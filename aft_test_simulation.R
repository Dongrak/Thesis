n_vector=c(100,200,500,1000)
sim_vector=c(100,150,200,250)

#n_vector=c(10,20,50,100)
#sim_vector=c(10,15,20,25)

result_aft=list()
result_cox=list()

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
  
  # dataset_What()
  dataset_What_aft=simulation_What(sim,beta_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,given_tol)
  dataset_What_cox=simulation_What(sim,beta_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,given_tol)
  
  # dataset_W()
  dataset_W_aft=simulation_W(beta_hat_aft,T_aft,D_aft,Z_aft,given_weight,given_test,dataset_What_aft)
  dataset_W_cox=simulation_W(beta_hat_cox,T_cox,D_cox,Z_cox,given_weight,given_test,dataset_What_cox)
  
  dataset_What50_aft=dataset_What_aft[1:(n*50),]
  dataset_What50_cox=dataset_What_cox[1:(n*50),]
  
  # PLOT : W_aft vs What_aft
  Figure1_W_aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_aft
  
  # PLOT : W_cox vs What_cox
  Figure1_W_cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=W),colour="tomato")
  #Figure1_W_cox
  
  # PLOT : std.W_aft vs std.What_aft
  Figure1_std.W_aft=
    ggplot()+
    geom_line(data=dataset_What50_aft,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_aft,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_aft
  
  # PLOT : W_cox vs What_cox
  Figure1_std.W_cox=
    ggplot()+
    geom_line(data=dataset_What50_cox,aes(x=t_i,y=std.What,group=group),colour="grey",alpha=0.5)+
    geom_line(data=dataset_W_cox,aes(x=t_i,y=std.W),colour="tomato")
  #Figure1_std.W_cox

  dataset_What1_aft=dataset_What_aft[1:(n*sim_vector[1]),]
  dataset_What1_cox=dataset_What_cox[1:(n*sim_vector[1]),]
  
  dataset_What2_aft=dataset_What_aft[1:(n*sim_vector[2]),]
  dataset_What2_cox=dataset_What_cox[1:(n*sim_vector[2]),]
  
  dataset_What3_aft=dataset_What_aft[1:(n*sim_vector[3]),]
  dataset_What3_cox=dataset_What_cox[1:(n*sim_vector[3]),]
  
  dataset_What4_aft=dataset_What_aft[1:(n*sim_vector[4]),]
  dataset_What4_cox=dataset_What_cox[1:(n*sim_vector[4]),]
  
  result_aft[[k]]=list(Figure1_W_aft,Figure1_std.W_aft,
                           dataset_W_aft,dataset_What_aft,
                           kolmogorov(dataset_W_aft,dataset_What1_aft),
                           kolmogorov(dataset_W_aft,dataset_What2_aft),
                           kolmogorov(dataset_W_aft,dataset_What3_aft),
                           kolmogorov(dataset_W_aft,dataset_What4_aft))
  
  result_cox[[k]]=list(Figure1_W_cox,Figure1_std.W_cox,
                           dataset_W_cox,dataset_What_cox,
                           kolmogorov(dataset_W_cox,dataset_What1_cox),
                           kolmogorov(dataset_W_cox,dataset_What2_cox),
                           kolmogorov(dataset_W_cox,dataset_What3_cox),
                           kolmogorov(dataset_W_cox,dataset_What4_cox))
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

