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
# install.packages("Rcpp")
# install.packages("RcppArmadillo")
# install.packages("ENmisc")
# install.packages("plotly")

library(ggplot2)
library(gridExtra)
library(survival)
library(aftgee)
library(Rcpp)
library(RcppArmadillo)
# library(ENmisc)
# library(plotly)

simulation=50
n=250
path=200
alpha=0.05

given_tol=0.1

#-------------------------------------------------------------
#-----------------------TEST STATISTICS-----------------------
#-------------------------------------------------------------

#-------------------------------------------------------------
#------------------------LINK FUNCTION------------------------
#-------------------------------------------------------------
simulation_link=function(simulation,n,path,alpha,tol){
  #simulation=simulation;n=n;path=path;alpha=alpha;tol=given_tol;
  
  result=list(NA)
  
  for(k in 1:simulation){
    if(k%%1==0) {
      cat("simulation",k,"\n")
    }
    # -------------------------------------------------------------
    # ------------------------DATA GENERATE------------------------
    # -------------------------------------------------------------
    # n=500
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
    
    #result[[k]]=list(result_ln_aft,result_ln_cox,p_value)
    result[[k]]=list(p_value)
  }
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
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_link_n250p150sim500")
date()
simulation_result_link2=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link2)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_link_n250p150sim100")
date()
simulation_result_link3=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link3)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_link_n250p150sim100")
date()
simulation_result_link4=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link4)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_link_n250p150sim100")
date()
simulation_result_link5=simulation_link(simulation,n,path,alpha,given_tol)
prob.table_link(simulation_result_link5)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_link_n250p150sim100")
date()
simulation_result_link=c(simulation_result_link1,simulation_result_link2,
                          simulation_result_link3,simulation_result_link4,
                          simulation_result_link5)
prob.table_link(simulation_result_link)
save.image("C:\\Users\\WOOJUNG\\Desktop\\simulation_result\\simulation_result_link_n250p150sim500")
date()


