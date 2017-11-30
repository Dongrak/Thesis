#-------------------------------------------------------------
#---------------------------SETTING---------------------------
#-------------------------------------------------------------
#rm(list=ls());gc();

memory.limit(16*2^20)

options(max.print=999999)
options(error=NULL)

# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("survival")
# install.packages("aftgee")
# install.packages("doParallel")
# install.packages("Rcpp")
# install.packages("RcppArmadillo")
# install.packages("ENmisc")
# install.packages("plotly")

library(ggplot2)
library(gridExtra)
library(survival)
library(aftgee)
library(doParallel)
# library(Rcpp)
# library(RcppArmadillo)
# library(ENmisc)
# library(plotly)

simulation=100
n=250
path=200
alpha=0.05

beta_0=1

gamma_0=0.1
# gamma_0=0.2
# gamma_0=0.3
# gamma_0=0.4
# gamma_0=0.5

given_tol=1

#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------
afttest_link=function(path,b,std,Time,Delta,Covari,tol){
  # path=200;b=beta_hat_ln_aft;std=std_hat_ln_aft;Time=X_ln_aft;Delta=D_ln_aft;Covari=Z_ln_aft;tol=given_tol;
  # path=200;b=beta_hat_ln_cox;std=std_hat_ln_cox;Time=X_ln_cox;Delta=D_ln_cox;Covari=Z_ln_cox;tol=given_tol;
  # path=200;b=c(1.3,1.1);Covari=c(Z_ln_aft,Z_ln_aft^2-4*Z_ln_aft);
  
  n=length(Time) # the number of subjects
  p=length(b) # the number of parameters
  
  Covari=matrix(Covari,nrow=n)
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  order_Covari=apply(Covari,2,function(x){order(x)}) 
  pi_i_z=sapply(1:n,function(j){apply(sapply(1:p,function(i){c(rep(0,(which(
    order_Covari[,i]==j)-1)),rep(1,(n+1-which(order_Covari[,i]==j))))}),1,prod)},simplify=F)
  
  N_i_t=sapply(1:n,function(j){(e_i_beta>=e_i_beta[j])*Delta[j]},simplify=F)
  #N_i_t
  
  Y_i_t=sapply(1:n,function(j){(e_i_beta<=e_i_beta[j])*1},simplify=F)
  #Y_i_t
  
  N_d_t=Reduce('+',N_i_t)
  #N_d_t
  
  S_0_t=Reduce('+',Y_i_t)
  #S_0_t
  
  S_1_t=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,
                          as.list(data.frame(t(Covari))),SIMPLIFY=FALSE))
  #S_1_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t 
  
  Mhat_i_t=mapply("-",N_i_t,lapply(lapply(
    Y_i_t,'*',dLambdahat_0_t),cumsum),SIMPLIFY=FALSE)
  #Mhat_i_t
  
  Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
  #Mhat_i_inf
  
  obs_path=Reduce('+',mapply('*',pi_i_z,Mhat_i_inf,SIMPLIFY=FALSE))/sqrt(n)
  #obs_path
  
  S_pi_t.z=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,pi_i_z,SIMPLIFY=FALSE))
  #S_pi_t.z
  
  dMhat_i_t=lapply(Mhat_i_t,function(x){diff(c(0,x))})
  #dMhat_i_t
  
  #-----------------------------------------------------------
  #----------------------kernel Smoothing---------------------
  #-----------------------------------------------------------
  
  #-----------------------------g0----------------------------
  Ghat_0_t=1-exp(-Lambdahat_0_t)
  #Ghat_0_t
  
  dGhat_0_t=diff(c(0,Ghat_0_t))
  #dGhat_0_t
  
  ghat_0_t=(ksmooth(e_i_beta,dGhat_0_t,"normal",
                    bandwidth = 1.06*sd(dGhat_0_t)*n^(-0.2),x.points=e_i_beta)$y)
  #ghat_0_t
  
  ghat_t.z=sapply(1:p,function(j){Reduce('+',lapply(mapply('*',pi_i_z,Covari[,j],
                                                           SIMPLIFY=FALSE),function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n},simplify=F)
  #ghat_t.z
  
  #-----------------------------f0----------------------------
  Fhat_0_e=1-cumprod(1-(Delta/S_0_t))
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
  
  fhat_t.z=sapply(1:p,function(j){Reduce('+',lapply(mapply('*',pi_i_z,Delta*Covari[,j],
                                                           SIMPLIFY=FALSE),function(x,y){t(x%*%t(y))},ghat_0_t*Time))/n},simplify=F)
  #fhat_t.z
  
  fhat_inf.z=lapply(fhat_t.z,function(x){x[n,]})
  
  #-----------------------------------------------------------
  #--------Find Beta_hat_star by using optimize function------
  #-----------------------------------------------------------
  U_beta=function(beta_U){
    #beta_U=b;
    
    Time_U=Time;Delta_U=Delta;Covari_U=Covari;
    
    e_i_beta_U=as.vector(log(Time_U)+Covari_U%*%beta_U)
    
    order_resid_U=order(e_i_beta_U)
    
    Time_U=Time_U[order_resid_U]
    Covari_U=matrix(Covari_U[order_resid_U,],nrow=n)
    Delta_U=Delta_U[order_resid_U]
    e_i_beta_U=e_i_beta_U[order_resid_U]
    
    N_i_t_U=sapply(1:n,function(j){(e_i_beta_U>=e_i_beta_U[j])*Delta_U[j]},simplify=F)
    #N_i_t_U
    
    Y_i_t_U=sapply(1:n,function(j){(e_i_beta_U<=e_i_beta_U[j])*1},simplify=F)
    #Y_i_t_U
    
    dN_i_t_U=lapply(N_i_t_U,function(x){diff(c(0,x))})
    #dN_i_t_U
    
    S_0_t_U=Reduce('+',Y_i_t_U)
    #S_0_t_U
    
    S_1_t_U=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t_U,as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))
    #S_1_t_U
    
    U_inf_U=apply(S_0_t_U*Reduce('+',mapply(function(x,y){x%*%t(y)},dN_i_t_U,
                                            as.list(data.frame(t(Covari_U))),SIMPLIFY=FALSE))-S_1_t_U*Reduce('+',dN_i_t_U),2,sum)/n
    #U_inf_U
    
    return(U_inf_U)
  }
  #U_beta()
  
  #------------------------SAMPLE PATH-----------------------
  
  app_path=list(NA)
  
  co=detectCores(logical=FALSE)-1 # number of core if logical is False else it means thread
  registerDoParallel(co)
  cl=makeCluster(co)
  app_path=foreach(k=1:path,.inorder=FALSE) %dopar% {
    
    # path_check=ceiling(path/2)
    # for (k inC 1:path){
    # if(k%%path_check==0) {
    #   cat("Sample Path",k,"\n")
    # }
    
    # tolerance=tol+1 #initial value
    
    # while (tolerance>tol){
    
    phi_i=rnorm(n)
    #phi_i
    
    U_pi_phi_inf.z=apply(S_0_t*Reduce('+',mapply('*',mapply(function(x,y){x%*%t(y)},dMhat_i_t,
                                                            pi_i_z,SIMPLIFY=FALSE),phi_i,SIMPLIFY=FALSE))-S_pi_t.z*Reduce('+',mapply('*',
                                                                                                                                     dMhat_i_t,phi_i,SIMPLIFY=FALSE)),2,sum)/n
    #U_pi_phi_inf.z
    
    U_phi_inf=apply(S_0_t*Reduce('+',mapply('*',mapply(function(x,y){x%*%t(y)},dMhat_i_t,
                                                       as.list(data.frame(t(Covari))),SIMPLIFY=FALSE),phi_i,SIMPLIFY=FALSE))-S_1_t*
                      Reduce('+',mapply('*',dMhat_i_t,phi_i,SIMPLIFY=FALSE)),2,sum)/n
    #U_phi_inf
    
    if(p==1){
      beta_hat_s_list=optimize(function(BETA){sum((U_beta(BETA)-U_phi_inf)^2)},
                               c(b-2*std,b+2*std),tol = 1e-16)
      #beta_hat_s_list
      
      beta_hat_s=beta_hat_s_list$minimum
      #beta_hat_s
      
      tolerance=beta_hat_s_list$objective
      #tolerance
    }
    if(p>1){
      beta_hat_s_list=optim(b,function(BETA){sum((U_beta(BETA)-U_phi_inf)^2)})
      #beta_hat_s_list
      
      beta_hat_s=beta_hat_s_list$par
      #beta_hat_s
      
      tolerance=beta_hat_s_list$value
      #tolerance
    }
    # }
    
    e_i_beta_s=as.vector(log(Time)+Covari%*%beta_hat_s)
    
    order_resid_s=order(e_i_beta_s)
    
    Delta_s=Delta[order_resid_s]
    e_i_beta_s=e_i_beta_s[order_resid_s]
    
    N_i_t_s=sapply(1:n,function(j){(e_i_beta_s>=e_i_beta_s[j])*Delta_s[j]},simplify=F)
    #N_i_t_s
    
    Y_i_t_s=sapply(1:n,function(j){(e_i_beta_s<=e_i_beta_s[j])*1},simplify=F)
    #Y_i_t_s
    
    N_d_t_s=Reduce('+',N_i_t_s)
    
    S_0_t_s=Reduce('+',Y_i_t_s)
    #S_0_t_s
    
    J_t_s=(S_0_t_s>0)*1
    #J_t_s
    
    dN_d_t_s=diff(c(0,N_d_t_s))
    #dN_d_t_s
    
    Lambdahat_0_t_s=cumsum((J_t_s/S_0_t_s)*dN_d_t_s)
    #Lambdahat_0_t_s
    
    F.T.=U_pi_phi_inf.z/sqrt(n)
    S.T.=sqrt(n)*Reduce('+',mapply('*',mapply('+',fhat_inf.z,mapply(function(x){apply(x,2,sum)},
                                                                    lapply(ghat_t.z,'*',dLambdahat_0_t),SIMPLIFY=FALSE),SIMPLIFY=FALSE),
                                   (b-beta_hat_s),SIMPLIFY=FALSE))
    T.T.=apply((S_pi_t.z*diff(c(0,Lambdahat_0_t-Lambdahat_0_t_s))),2,sum)/sqrt(n)
    
    app_path=F.T.-S.T.-T.T.
    # app_path[[k]]=F.T.-S.T.-T.T.
    #app_path
  }
  stopCluster(cl)
  closeAllConnections()
  
  std.boot=apply(mapply(function(x){as.vector(x)},app_path),1,sd)
  # std.boot
  
  app_std.path=lapply(app_path,function(x){x/std.boot})
  # app_std.path
  
  obs_std.path=obs_path/std.boot
  # obs_std.path
  
  #-----------------------MAXIMUM VALUE-----------------------
  max_app_path=unlist(lapply(app_path,function(x){max(abs(x))}))
  # max_app_path
  
  max_obs_path=max(abs(obs_path))
  # max_obs_path
  
  max_app_std.path=unlist(lapply(app_std.path,function(x){max(abs(x))}))
  # max_app_std.path
  
  max_obs_std.path=max(abs(obs_std.path))
  # max_obs_std.path
  
  #--------------------------P VALUE--------------------------
  p_value=length(which((max_app_path>max_obs_path)*1==1))/path
  # p_value
  
  std.p_value=length(which((max_app_std.path>max_obs_std.path)*1==1))/path
  # std.p_value
  
  # result=list(Time,Delta,Covari,e_i_beta,std.boot,
  #             app_path,app_std.path,
  #             obs_path,obs_std.path,
  #             p_value,std.p_value)
  # 
  # names(result)=c("Time","Delta","Covari","Resid","std.boot",
  #                 "app_path","app_std.path",
  #                 "obs_path","obs_std.path",
  #                 "p_value","std.p_value")
  
  result=list(p_value,std.p_value);names(result)=c("p_value","std.p_value");
  # result
  
  rm(list=(ls()[ls()!="result"]));gc();
  return(result)
}
#afttest_link()

#-------------------------------------------------------------
#------------------------LINK FUNCTION------------------------
#-------------------------------------------------------------
simulation_link=function(simulation,n,path,alpha,tol){
  #simulation=simulation;n=n;path=path;alpha=alpha;tol=given_tol;
  
  result=list(NA)
  
  # co=detectCores(logical=FALSE)-1 # number of core if logical is False else it means thread
  # registerDoParallel(co)
  # cl=makeCluster(co) 
  # result=foreach(k=1:simulation,.packages=c('aftgee','survival','doParallel'),.export='afttest_link',.inorder=FALSE) %dopar% {
  for (k in 1:simulation) {
    if(k%%1==0) {
      cat("simulation",k,"\n")
    }
  
    # -------------------------------------------------------------
    # ------------------------DATA GENERATE------------------------
    # -------------------------------------------------------------
    Z=matrix(runif(n,0,5),nrow=n)

    #-------------------LOG NORMAL DISTRIBUTION-------------------
    T_ln_aft=as.vector(exp(-gamma_0*Z)*qlnorm(runif(n),5,1))
    C_ln_aft=as.vector(exp(-gamma_0*Z)*qlnorm(runif(n),6.5,1))
    X_ln_aft=C_ln_aft*(T_ln_aft>C_ln_aft)+T_ln_aft*(T_ln_aft<=C_ln_aft)
    D_ln_aft=0*(T_ln_aft>C_ln_aft)+1*(T_ln_aft<=C_ln_aft)
    Z_ln_aft=Z

    T_ln_aft_l=as.vector((1+gamma_0*Z)*qlnorm(runif(n),5,1))
    C_ln_aft_l=as.vector((1+gamma_0*Z)*qlnorm(runif(n),6.5,1))
    X_ln_aft_l=C_ln_aft_l*(T_ln_aft_l>C_ln_aft_l)+T_ln_aft_l*(T_ln_aft_l<=C_ln_aft_l)
    D_ln_aft_l=0*(T_ln_aft_l>C_ln_aft_l)+1*(T_ln_aft_l<=C_ln_aft_l)
    Z_ln_aft_l=Z
    
    #------------Estimate Beta_hat_ln_aft by using Aftgee-----------
    aftsrr_beta_ln_aft=aftsrr(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,method="nonsm")
    beta_hat_ln_aft=-as.vector(aftsrr_beta_ln_aft$beta);beta_hat_ln_aft
    std_hat_ln_aft=diag(aftsrr_beta_ln_aft$covmat$ISMB);std_hat_ln_aft
    
    aftsrr_beta_ln_aft_l=aftsrr(Surv(X_ln_aft_l,D_ln_aft_l)~Z_ln_aft_l,method="nonsm")
    beta_hat_ln_aft_l=-as.vector(aftsrr_beta_ln_aft_l$beta);beta_hat_ln_aft_l
    std_hat_ln_aft_l=diag(aftsrr_beta_ln_aft_l$covmat$ISMB);std_hat_ln_aft_l
    
    # result_ln_aft
    result_ln_aft=afttest_link(path,beta_hat_ln_aft,std_hat_ln_aft,
                                    X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
    
    result_ln_aft_l=afttest_link(path,beta_hat_ln_aft_l,std_hat_ln_aft_l,
                                      X_ln_aft_l,D_ln_aft_l,Z_ln_aft_l,given_tol)
    
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
    
    # result[[k]]=list(result_ln_aft,result_ln_aft_l,p_value)
    result[[k]]=list(p_value)
    # result=list(p_value)
  }
        # stopClustCer(cl)
  
  rm(list=(ls()[ls()!="result"]));gc();
  return(result)
}
#simulation_link

prob.table_link=function(simul_result){
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
#prob.table_link


date()
simulation_result_link1=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link1)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link2=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link2)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",2*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link3=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link3)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",3*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link4=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link4)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",4*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link5=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link5)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",5*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link6=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link6)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",6*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link7=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link7)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",7*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link8=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link8)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",8*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link9=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link9)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",9*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link0=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link0)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",10*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()
simulation_result_link=c(simulation_result_link1,simulation_result_link2,
                         simulation_result_link3,simulation_result_link4,
                         simulation_result_link5,simulation_result_link6,
                         simulation_result_link7,simulation_result_link8,
                         simulation_result_link9,simulation_result_link0)
prob.table_link(simulation_result_link)
save.image(paste0(getwd(),"/result_link_1","_n",n,"p",path,"s",10*simulation,"g",gsub(".","",gamma_0,fixed = TRUE)))
date()



