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

#-------------------------------------------------------------
#---------------------------EXAMPLE---------------------------
#-------------------------------------------------------------
given_tol=100

k=length(which(is.na(pbc$trt)==0))
ID=pbc$id[1:k]
D_pbc=pbc$status[1:k] # paper
D_pbc[which(D_pbc==1)]=0
D_pbc[which(D_pbc==2)]=1
X_pbc=pbc$time[1:k]

Z1_pbc=pbc$bili[1:k]
lZ1_pbc=log(Z1_pbc)
Z2_pbc=pbc$protime[1:k]
lZ2_pbc=log(Z2_pbc)
Z3_pbc=pbc$albumin[1:k]
lZ3_pbc=log(Z3_pbc)
Z4_pbc=pbc$age[1:k]
lZ4_pbc=log(Z4_pbc)
Z5_pbc=pbc$edema[1:k]

###############################################################################
Z_pbc=c(lZ1_pbc,lZ2_pbc,lZ3_pbc,Z4_pbc,Z5_pbc)

###
# aftsrr_beta_pbc=aftsrr(Surv(X_pbc,D_pbc)~lZ1_pbc+lZ2_pbc+lZ3_pbc+Z4_pbc+Z5_pbc,method="nonsm")
# beta_hat_pbc=-as.vector(aftsrr_beta_pbc$beta);beta_hat_pbc
# se_hat_pbc=diag(aftsrr_beta_pbc$covmat$ISMB);se_hat_pbc

###
# result_pbc_omni=afttest_omni(1000,beta_hat_pbc,se_hat_pbc,
#                              X_pbc,D_pbc,Z_pbc,given_tol)
# result_pbc_omni$p_value
# result_pbc_omni$std.p_value
# plotting_omni(result_pbc_omni,1,50,"unstd")
# plotting_omni(result_pbc_omni,1,50,"std")

###
# result_pbc_form=afttest_form(1000,beta_hat_pbc,se_hat_pbc,
#                              X_pbc,D_pbc,Z_pbc,given_tol,form=4)
result_pbc_form$p_value
result_pbc_form$Beta
result_pbc_form$SE
plotting_form(result_pbc_form,50)

cairo_ps("afttest_form_real_f.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
plotting_form(result_pbc_form,50)
dev.off()

###
# result_pbc_link=afttest_link(1000,beta_hat_pbc,se_hat_pbc,
#                              X_pbc,D_pbc,Z_pbc,given_tol)
# result_pbc_link$p_value
# plotting_link(result_pbc_link,1,50)

###############################################################################
# Z_pbc_log=c(lZ1_pbc,lZ2_pbc,lZ3_pbc,lZ4_pbc,Z5_pbc)

# aftsrr_beta_pbc_log=aftsrr(Surv(X_pbc,D_pbc)~lZ1_pbc+lZ2_pbc+lZ3_pbc+lZ4_pbc+Z5_pbc,method="nonsm")
# beta_hat_pbc_log=-as.vector(aftsrr_beta_pbc_log$beta);beta_hat_pbc_log
# se_hat_pbc_log=diag(aftsrr_beta_pbc_log$covmat$ISMB);se_hat_pbc_log

###
# result_pbc_omni_log=afttest_omni(1000,beta_hat_pbc_log,se_hat_pbc_log,
#                                  X_pbc,D_pbc,Z_pbc_log,given_tol)
# result_pbc_omni_log$p_value
# result_pbc_omni_log$std.p_value
# plotting_omni(result_pbc_omni_log,1,50,"unstd")
# plotting_omni(result_pbc_omni_log,1,50,"std")

 ###
# result_pbc_form_log=afttest_form(1000,beta_hat_pbc_log,se_hat_pbc_log,
#                                  X_pbc,D_pbc,Z_pbc_log,given_tol,form=4)
result_pbc_form_log$p_value
result_pbc_form_log$Beta
result_pbc_form_log$SE
plotting_form(result_pbc_form_log,50)

cairo_ps("afttest_form_real_f_log.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
plotting_form(result_pbc_form_log,50)
dev.off()

###
# result_pbc_link_log=afttest_link(1000,beta_hat_pbc_log,se_hat_pbc_log,
#                                  X_pbc,D_pbc,Z_pbc_log,given_tol)
# result_pbc_link_log$p_value
# plotting_link(result_pbc_link_log,1,50)
