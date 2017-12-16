path_vector=c(200,500,1000,1500,2000)#,2500,3000)

date()
a200=afttest_omni(path_vector[1],beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
a200$p_value
a200$std.p_value
save.image(paste0(getwd(),"/a200"))
# plotting_omni(a200,"rank",30)
# plot_ly(z=a200$app_path[[1]],type="surface")
# plot_ly(z=a200$obs_std.path,type="surface")

date()
b500=afttest_omni(path_vector[2],beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
b500$p_value
b500$std.p_value
save.image(paste0(getwd(),"/b500"))
# plotting_omni(b500,"rank",50)

date()
c1000=afttest_omni(path_vector[3],beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
c1000$p_value
c1000$std.p_value
save.image(paste0(getwd(),"/c1000"))
# plotting_omni(c1000,"rank",50)

date()
d1500=afttest_omni(path_vector[4],beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
d1500$p_value
d1500$std.p_value
save.image(paste0(getwd(),"/d1500"))
# plotting_omni(d1500,"rank",50)

date()
e2000=afttest_omni(path_vector[5],beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
e2000$p_value
e2000$std.p_value
save.image(paste0(getwd(),"/e2000"))
# plotting_omni(e2000,"rank",50)

# date()
# f2500=afttest_omni(path_vector[6],beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
# f2500$p_value
# f2500$std.p_value
# save.image(paste0(getwd(),"/1201data"))
# # plotting_omni(f2500,"rank",50)
# 
# date()
# g3000=afttest_omni(path_vector[7],beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
# g3000$p_value
# g3000$std.p_value
# # plotting_omni(g3000,"rank",50)

date()
p_vector=c(a200$p_value,b500$p_value,c1000$p_value,d1500$p_value,
           e2000$p_value)#,f2500$p_value,g3000$p_value)
std.p_vector=c(a200$std.p_value,b500$std.p_value,c1000$std.p_value,d1500$std.p_value,
               e2000$std.p_value)#,f2500$std.p_value,g3000$std.p_value)
var_list=list(a200$std.boot,b500$std.boot,c1000$std.boot,d1500$std.boot,
              e2000$std.boot)#,f2500$std.boot,g3000$std.boot)
save.image(paste0(getwd(),"/1201data"))

############################################################################
############################################################################
n=250
beta_0=1
Z=matrix(rexp(n,1/10),nrow=n)
Z_log=log(Z)

T_log=as.vector(exp(-beta_0*Z_log+rnorm(n,5,1)))
C_log=as.vector(exp(-beta_0*Z_log+rnorm(n,6.5,1)))
X_log=C_log*(T_log>C_log)+T_log*(T_log<=C_log)
D_log=0*(T_log>C_log)+1*(T_log<=C_log)

############################################################################
aftsrr_beta=aftsrr(Surv(X_log,D_log)~Z,method="nonsm")
beta_hat=-as.vector(aftsrr_beta$beta);beta_hat
std_hat=diag(aftsrr_beta$covmat$ISMB);std_hat

result=afttest_form(1000,beta_hat,std_hat,X_log,D_log,Z,0.1)
result$p_value
plotting_form(result,1,50)

cairo_ps("afttest_form.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
plotting_form(result,1,100)
dev.off()

############################################################################
aftsrr_beta_log=aftsrr(Surv(X_log,D_log)~Z_log,method="nonsm")
beta_hat_log=-as.vector(aftsrr_beta_log$beta);beta_hat_log
std_hat_log=diag(aftsrr_beta_log$covmat$ISMB);std_hat_log

result_log=afttest_form(1000,beta_hat_log,std_hat_log,X_log,D_log,Z_log,0.1)
result_log$p_value
plotting_form(result_log,1,50)

cairo_ps("afttest_form_log.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
plotting_form(result_log,1,100)
dev.off()

##################################################################################
# n=250
# beta_0=1
# gamma_0=0.1
# 
# Z1=matrix(rexp(n,1/10),nrow=n)
# Z2=matrix(runif(n,0,10),nrow=n)
# 
# T_l=as.vector((1+beta_0*Z1+gamma_0*Z2+rnorm(n,5,1)))
# C_l=as.vector((1+beta_0*Z1+gamma_0*Z2+rnorm(n,6.5,1)))
# X_l=C_l*(T_l>C_l)+T_l*(T_l<=C_l)
# D_l=0*(T_l>C_l)+1*(T_l<=C_l)
# Z1_l=Z1
# Z2_l=Z2
# Z_l=cbind(Z1_l,Z2_l)
# 
# aftsrr_beta_l=aftsrr(Surv(X_l,D_l)~Z_l,method="nonsm")
# beta_hat_l=-as.vector(aftsrr_beta_l$beta);beta_hat_l
# std_hat_l=diag(aftsrr_beta_l$covmat$ISMB);std_hat_l
# 
# result_l=afttest_link(200,beta_hat_l,std_hat_l,X_l,D_l,Z_l,0.1)
# result_l$p_value
# plotting_link(result_l,1,100)
# 
# cairo_ps("afttest_link_l.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# plotting_link(result_l,1,100)
# dev.off()
# 
# aftsrr_beta_lz=aftsrr(Surv(X_l,D_l)~Z_l,method="nonsm")
# beta_hat_lz=-as.vector(aftsrr_beta_l$beta);beta_hat_lz
# std_hat_lz=diag(aftsrr_beta_l$covmat$ISMB);std_hat_lz
# 
# result_lz=afttest_link(200,beta_hat_lz,std_hat_lz,X_l,D_l,Z_l,0.1)
# result_lz$p_value
# plotting_link(result_lz,1,100)
# 
# length(result_lz$std.boot)
# 
# cairo_ps("afttest_link_lz.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# plotting_link(result_lz,1,100)
# dev.off()
##################################################################################
# n=250
# beta_0=1
# gamma_0=0.1
# 
# Z1=matrix(rexp(n,1/10),nrow=n)
# Z2=matrix(runif(n,0,10),nrow=n)
# 
# T_aft=as.vector(exp(-beta_0*Z1-gamma_0*Z2+rnorm(n,5,1)))
# C_aft=as.vector(exp(-beta_0*Z1-gamma_0*Z2+rnorm(n,6.5,1)))
# X_aft=C_aft*(T_aft>C_aft)+T_aft*(T_aft<=C_aft)
# D_aft=0*(T_aft>C_aft)+1*(T_aft<=C_aft)
# Z1_aft=Z1
# Z2_aft=Z2
# Z_aft=cbind(Z1_aft,Z2_aft)
# 
# aftsrr_beta_aft=aftsrr(Surv(X_aft,D_aft)~Z_aft,method="nonsm")
# beta_hat_aft=-as.vector(aftsrr_beta_aft$beta);beta_hat_aft
# std_hat_aft=diag(aftsrr_beta_aft$covmat$ISMB);std_hat_aft
# 
# result_aft_z=afttest_link(200,beta_hat_aft,std_hat_aft,X_aft,D_aft,Z_aft,0.1)
# result_aft_z$p_value
# plotting_link(result_aft_z,1,50)
# 
# cairo_ps("afttest_link_aft_z.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# plotting_link(result_aft_z,1,100)
# dev.off()
# 
# result_aft=afttest_link(200,beta_hat_aft,std_hat_aft,X_aft,D_aft,Z_aft,0.1)
# result_aft$p_value
# plotting_link(result_aft,1,50)
# 
# cairo_ps("afttest_link_aft.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# plotting_link(result_aft,1,100)
# dev.off()


# ############################################################################
# path=200;b=beta_hat_f;std=std_hat_f;Time=X_f;Delta=D_f;Covari=logZ;form=1
# 
# n=length(Time) # the number of subjects
# p=length(b) # the number of parameters
# 
# Covari=matrix(Covari,nrow=n)
# 
# e_i_beta=as.vector(log(Time)+Covari%*%b)
# 
# order_resid=order(e_i_beta)
# 
# Time=Time[order_resid]
# Covari=matrix(Covari[order_resid,],nrow=n)
# Delta=Delta[order_resid]
# e_i_beta=e_i_beta[order_resid] 
# 
# e_i_beta_form=as.vector(Covari[,form]*b[form])
# order_e_i_beta_form=order(e_i_beta_form)
# 
# # weight function
# order_Covari_form=order(Covari[,form])
# pi_i_z=sapply(1:n,function(j){c(rep(0,(which(order_Covari_form==j)-1)),
#                                 rep(1,(n+1-which(order_Covari_form==j))))},simplify=F)
# 
# N_i_t=sapply(1:n,function(j){(e_i_beta>=e_i_beta[j])*Delta[j]},simplify=F)
# #N_i_t
# 
# Y_i_t=sapply(1:n,function(j){(e_i_beta<=e_i_beta[j])*1},simplify=F)
# #Y_i_t
# 
# N_d_t=Reduce('+',N_i_t)
# #N_d_t
# 
# S_0_t=Reduce('+',Y_i_t)
# #S_0_t
# 
# S_1_t=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,
#                         as.list(data.frame(t(Covari))),SIMPLIFY=FALSE))
# #S_1_t
# 
# J_t=(S_0_t>0)*1
# #J_t
# 
# dN_d_t=diff(c(0,N_d_t))
# #dN_d_t
# 
# Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
# #Lambdahat_0_t
# 
# dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
# #dLambdahat_0_t 
# 
# Mhat_i_t=mapply("-",N_i_t,lapply(lapply(
#   Y_i_t,'*',dLambdahat_0_t),cumsum),SIMPLIFY=FALSE)
# #Mhat_i_t
# 
# Mhat_i_inf=unlist(lapply(Mhat_i_t,function(x){x[n]}))
# #Mhat_i_inf
# 
# obs_path=Reduce('+',mapply('*',pi_i_z,Mhat_i_inf,SIMPLIFY=FALSE))/sqrt(n)
# #obs_path
# 
# plot(obs_path,type="s")
# 
# obs_path_logZ=obs_path[order(logZ)]
# logZ_Z=Z[order(logZ)]
# 
# plot(logZ_Z,obs_path_logZ,type="s")
# 
# # obs_path_Z=obs_path[order(Z)]
# # Z_Z=Z[order(Z)]
# # 
# # plot(Z_Z,obs_path_Z,type="s")


##################################################################################
n=250
beta_0=1
Z=matrix(rexp(n,1/10),nrow=n)
Z_log=log(Z)
given_tol=1

T_log=as.vector(exp(-beta_0*Z_log+rnorm(n,5,1)))
C_log=as.vector(exp(-beta_0*Z_log+rnorm(n,6.5,1)))
X_log=C_log*(T_log>C_log)+T_log*(T_log<=C_log)
D_log=0*(T_log>C_log)+1*(T_log<=C_log)

########################------------------------------------------
# aftsrr_beta_omni=aftsrr(Surv(X_log,D_log)~Z_log,method="nonsm")
# beta_hat_omni=-as.vector(aftsrr_beta_omni$beta);beta_hat_omni
# std_hat_omni=diag(aftsrr_beta_omni$covmat$ISMB);std_hat_omni
# 
# result_omni=afttest_omni(1000,beta_hat_omni,std_hat_omni,X_log,D_log,Z_log,given_tol)
# result_omni$p_value
# result_omni$std.p_value
# # cairo_ps("afttest_omni_std.eps",onefile=F,
# #          height=4,width=8, fallback_resolution = 600)
# # plotting_omni(result_omni,1,100,"std")
# # dev.off()

########################------------------------------------------
result_omni_no=afttest_omni(1000,beta_hat_omni,std_hat_omni,X_log,D_log,Z_log,1000000000)
result_omni_no$p_value
result_omni_no$std.p_value

cairo_ps("afttest_omni_no_unstd.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
plotting_omni(result_omni_no,1,100,"unstd")
dev.off()

########################------------------------------------------
# aftsrr_beta_omni_false=aftsrr(Surv(X_log,D_log)~Z,method="nonsm")
# beta_hat_omni_false=-as.vector(aftsrr_beta_omni_false$beta);beta_hat_omni_false
# std_hat_omni_false=diag(aftsrr_beta_omni_false$covmat$ISMB);std_hat_omni_false
# 
# date()
# result_omni_false=afttest_omni(1000,beta_hat_omni_false,std_hat_omni_false,X_log,D_log,Z,given_tol)
# date()
# result_omni_false$p_value
# result_omni_false$std.p_value
# 
# # cairo_ps("afttest_omni_false_unstd.eps",onefile=F,
# #          height=4,width=8, fallback_resolution = 600)
# # plotting_omni(result_omni_false,1,100,"unstd")
# # dev.off()


