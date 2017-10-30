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

#-------------------------------------------------------------
#------------------------WEIGHT&TOLERANCE---------------------
#-------------------------------------------------------------
path=200

given_tol=0.1

#------------------------DATA GENERATION----------------------
n=200
beta_0=1
gamma_0=0.1
Z=matrix(rnorm(n,3,1),nrow=n)

#---------------------WEIBULL DISTRIBUTION--------------------
# T ~ Generalized Gamma(alpha=1,beta,sigma) i.e. Weibull ~ (beta,sigma)
# alpha : shape // beta : scale // sigma : rate parameter of generaized gamma
alpha_wb_T=1;alpha_wb_C=1     # alpha must be one when weibull
beta_wb_T=100;beta_wb_C=101       #
sigma_wb_T=100;sigma_wb_C=100     #

V_wb_T=log(beta_wb_T)+log(qgamma(runif(n),alpha_wb_T,1))/sigma_wb_T
V_wb_C=log(beta_wb_C)+log(qgamma(runif(n),alpha_wb_C,1))/sigma_wb_C
T_wb=exp(-Z%*%beta_0+V_wb_T)
C_wb=exp(-Z%*%beta_0+V_wb_C)
X_wb=C_wb*(T_wb>C_wb)+T_wb*(T_wb<=C_wb) # observed failure time
D_wb=0*(T_wb>C_wb)+1*(T_wb<=C_wb)
Z_wb=Z # delta 0:censored & 1:observed

#---------------GNERALIZED GUMBELL DISTRIBUTION---------------
# T ~ Generalized Gamma(alpha,beta,sigma) if. V ~ Generalized Gumbell(beta,sigma)
# alpha : shape // beta : scale // sigma : rate parameter of generaized gamma
alpha_gg_T=100;alpha_gg_C=110    # alpha must be greater than one whne not weibull
beta_gg_T=0.1;beta_gg_C=0.1      # gamma distribution
sigma_gg_T=1;sigma_gg_C=1        #

V_gg_T=log(beta_gg_T)+log(qgamma(runif(n),alpha_gg_T,1))/sigma_gg_T
V_gg_C=log(beta_gg_C)+log(qgamma(runif(n),alpha_gg_C,1))/sigma_gg_C
T_gg=exp(-Z%*%beta_0+V_gg_T)
C_gg=exp(-Z%*%beta_0+V_gg_C)
X_gg=C_gg*(T_gg>C_gg)+T_gg*(T_gg<=C_gg) # observed failure time
D_gg=0*(T_gg>C_gg)+1*(T_gg<=C_gg)       # delta 0:censored & 1:observed
Z_gg=Z 

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
Z_wb_f=Z_f

#-------------------------------------------------------------
#-------------Estimate Beta_hat_wb by using Aftgee------------
#-------------------------------------------------------------
aftsrr_beta_wb=aftsrr(Surv(X_wb,D_wb)~Z_wb,method="nonsm")
beta_hat_wb=-as.vector(aftsrr_beta_wb$beta);beta_hat_wb
std_hat_wb=diag(aftsrr_beta_wb$covmat$ISMB);std_hat_wb

#-------------------------------------------------------------
#------------Estimate Beta_hat_gev by using Aftgee------------
#-------------------------------------------------------------
aftsrr_beta_gg=aftsrr(Surv(X_gg,D_gg)~Z_gg,method="nonsm")
beta_hat_gg=-as.vector(aftsrr_beta_gg$beta);beta_hat_gg
std_hat_gg=diag(aftsrr_beta_gg$covmat$ISMB);std_hat_gg

#-------------------------------------------------------------
#------------Estimate Beta_hat_wb_f by using Aftgee-----------
#-------------------------------------------------------------
aftsrr_beta_wb_f=aftsrr(Surv(X_wb_f,D_wb_f)~Z_wb_f,method="nonsm")
beta_hat_wb_f=-as.vector(aftsrr_beta_wb_f$beta);beta_hat_wb_f
std_hat_wb_f=diag(aftsrr_beta_wb_f$covmat$ISMB);std_hat_wb_f

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

T_ln_aft_f=as.vector(exp(-beta_0*Z-gamma_0*(Z^2))*qlnorm(runif(n),5,1))
C_ln_aft_f=as.vector(exp(-beta_0*Z-gamma_0*(Z^2))*qlnorm(runif(n),6.5,1))
X_ln_aft_f=C_ln_aft_f*(T_ln_aft_f>C_ln_aft_f)+T_ln_aft_f*(T_ln_aft_f<=C_ln_aft_f)
D_ln_aft_f=0*(T_ln_aft_f>C_ln_aft_f)+1*(T_ln_aft_f<=C_ln_aft_f)
Z_ln_aft_f=Z

#------------Estimate Beta_hat_wb_f by using Aftgee-----------
aftsrr_beta_ln_aft=aftsrr(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,method="nonsm")
beta_hat_ln_aft=-as.vector(aftsrr_beta_ln_aft$beta);beta_hat_ln_aft
std_hat_ln_aft=diag(aftsrr_beta_ln_aft$covmat$ISMB);std_hat_ln_aft

aftsrr_beta_ln_cox=aftsrr(Surv(X_ln_cox,D_ln_cox)~Z_ln_cox,method="nonsm")
beta_hat_ln_cox=-as.vector(aftsrr_beta_ln_cox$beta);beta_hat_ln_cox
std_hat_ln_cox=diag(aftsrr_beta_ln_cox$covmat$ISMB);std_hat_ln_cox

aftsrr_beta_ln_aft_f=aftsrr(Surv(X_ln_aft_f,D_ln_aft_f)~Z_ln_aft_f,method="nonsm")
beta_hat_ln_aft_f=-as.vector(aftsrr_beta_ln_aft_f$beta);beta_hat_ln_aft_f
std_hat_ln_aft_f=diag(aftsrr_beta_ln_aft_f$covmat$ISMB);std_hat_ln_aft_f