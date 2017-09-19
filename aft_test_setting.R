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
#install.packages("ENmisc")

library(ggplot2)
library(survival)
library(aftgee)
library(ENmisc)

#-------------------------------------------------------------
#------------------------WEIGHT&TOLERANCE---------------------
#-------------------------------------------------------------
path=200

given_tol=1

given_weight="c"
#(weight=="a"){w_i=Covari*(Covari<=median(Covari))}
#(weight=="b"){w_i=Covari}  
#(weight=="c"){w_i=1*(Covari<=median(Covari))}
#(weight=="d"){w_i=1}

given_test="omni"
# given_test=="link.ftn"
# given_test=="ftn.form"
# given_test=="....???"

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
D_wb=0*(T_wb>C_wb)+1*(T_wb<=C_wb)       # delta 0:censored & 1:observed

D_wb=D_wb[order(X_wb)]
Z_wb=Z[order(X_wb)]
X_wb=X_wb[order(X_wb)]
X_wb
length(which(D_wb==0))/n

#---------------GNERALIZED GUMBELL DISTRIBUTION---------------
# T ~ Generalized Gamma(alpha,beta,sigma) if. V ~ Generalized Gumbell(beta,sigma)
# alpha : shape // beta : scale // sigma : rate parameter of generaized gamma
alpha_gg_T=100;alpha_gg_C=110    # alpha must be greater than one whne not weibull
beta_gg_T=0.1;beta_gg_C=0.1       # gamma distribution
sigma_gg_T=1;sigma_gg_C=1       #

V_gg_T=log(beta_gg_T)+log(qgamma(runif(n),alpha_gg_T,1))/sigma_gg_T
V_gg_C=log(beta_gg_C)+log(qgamma(runif(n),alpha_gg_C,1))/sigma_gg_C
T_gg=exp(-Z%*%beta_0+V_gg_T)
C_gg=exp(-Z%*%beta_0+V_gg_C)
X_gg=C_gg*(T_gg>C_gg)+T_gg*(T_gg<=C_gg) # observed failure time
D_gg=0*(T_gg>C_gg)+1*(T_gg<=C_gg)       # delta 0:censored & 1:observed

D_gg=D_gg[order(X_gg)]
Z_gg=Z[order(X_gg)]
X_gg=X_gg[order(X_gg)]
X_gg
length(which(D_gg==0))/n

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

#-------------------------------------------------------------
#------------Estimate Beta_hat_wb_f by using Aftgee-----------
#-------------------------------------------------------------
aftsrr_beta_wb_f=aftsrr(Surv(X_wb_f,D_wb_f)~Z_wb_f,method="nonsm")
beta_hat_wb_f=-unlist(summary(aftsrr_beta_wb_f))$coefficients1;beta_hat_wb_f
std_hat_wb_f=unlist(summary(aftsrr_beta_wb_f))$coefficients2;std_hat_wb_f

