order_Covari_form=order(Covari[,form])
pi_i_z=sapply(1:n,function(j){c(rep(0,(which(order_Covari_form==j)-1)),
                                rep(1,(n+1-which(order_Covari_form==j))))},simplify=F)

order_Covari=apply(Covari,2,function(x){order(x)}) 

k=167
length(which(pi_i_z1[[k]]==pi_i_z[[k]]))


pi_i_z=as.list(data.frame(t(matrix(unlist(pi_i_z),nrow=n))))

aa=matrix(c(1,3,9,7,5,2,8,6),nrow=4);aa
which(aa==max(aa),arr.ind=TRUE)
aa[1:2,]

aa=matrix(c(2:5,3,1,4,5),nrow=4)
bb=as.list(data.frame(t(aa)))
cc=matrix(c(4:7,-1:2),nrow=4)
dd=lapply(bb,function(x){t(x-t(cc))})
lapply(matrix(c(2:5,3,1,4,5),nrow=4),function(x){x})
length(lapply(Covari,function(x){t(x-t(E_t))}))
lapply(matrix(c(2:5,3:6),nrow=4),function(x){x-matrix(c(4:7,-1:2),nrow=4)})


pi_ij_z=list(NA)
pi_i_z=list(NA)
for(j in 1:p){
  Covari_order_j=as.vector(Covari[,j])[order(Covari[,j])]
  for(i in 1:n){
    pi_i_z[[i]]=((Covari_order_j>=Covari_order_j[i]))[order_resid]
  }
  pi_ij_z[[j]]=pi_i_z
}

aa=matrix(c((1:10),c(rep(1,4),rep(0,6))),10,2)


mapply(function(x,y){x%*%t(y)},list(c(0:3),c(4:7)),list(c(1:4),c(-1:2)),SIMPLIFY=FALSE)



aa=list(matrix(0:3,nrow=2),matrix(((0:3)*2),nrow=2));aa

apply(matrix(0:5,nrow=3),2,cumsum)
rowSums()

aa=list(matrix(0:3,nrow=2),matrix(((0:3)*2),nrow=2));aa
cc=list(matrix(((-1:2)*2),nrow=2),matrix((-1:2),nrow=2));cc
dd=list(matrix(((-2:1)),nrow=2),matrix((3:6)/2,nrow=2));cc
ee=list(aa,cc,dd);ee

Reduce("pmax",lapply(ee[[3]],function(x){abs(x)}))


Reduce("pmax",mapply(function(x){abs(x[[1]])},ee,SIMPLIFY=FALSE))

Reduce("pmax",mapply(function(x){abs(x[[1]])},ee,SIMPLIFY=FALSE))
lapply(ee[[1]],function(x){x/bb[[1]]})
matrix(apply(mapply(function(x){as.vector(x)},ee[[1]]),1,sd),nrow=2)
is.matrix(dd[[2]])

lapply(ee,function(x,y){x[[2]]/y[[2]]},aa)
mapply(function(x){mean(x[[1]])},dd,SIMPLIFY=FALSE)
mapply(function(x){mapply("pmax",x,SIMPLIFY=FALSE)},dd,SIMPLIFY=FALSE)
Reduce("pmax",dd[[1]])
lapply(dataset_What,function(x,y){x[[j]]/y[[j]]},std.boot)
lapply(dd,function(x){pmax(x[[j]])})

##############################################################
# omnibus test
# maxmax_What=c(NA)
# for(i in 1:30){
#   aa=What_omni(beta_hat_ln_cox,std_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,10000)
#   maxmax_What[i]=max(abs(aa))
#   plot(aa[,100],type="s",col="grey",ylim=c(-1,1));par(new=TRUE)
# }
# bb=W_omni(beta_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox)$obs_stat_omni
# maxmax_W=max(abs(bb))
# plot(bb[,100],type="s",col="red",ylim=c(-1,1))
# abline(h=0)
# length(which((maxmax_What>=maxmax_W)==1))/30
##############################################################
# functional form
# maxmax_What=c(NA)
# for(i in 1:30){
#   aa=What_form(beta_hat_ln_cox,std_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,10000,form=1)
#   maxmax_What[i]=max(abs(aa))
#   plot(aa,type="s",col="grey",ylim=c(-5,5));par(new=TRUE)
# }
# bb=W_form(beta_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,form=1)$obs_stat_form
# maxmax_W=max(abs(bb))
# plot(bb,type="s",col="red",ylim=c(-5,5))
# abline(h=0)
# length(which((maxmax_What>=maxmax_W)==1))/30
##############################################################
# link function
# maxmax_What=c(NA)
# for(i in 1:30){
#   aa=What_link(beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,10000)
#   maxmax_What[i]=max(abs(aa))
#   plot(aa,type="s",col="grey",ylim=c(-5,5));par(new=TRUE)
# }
# bb=W_link(beta_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft)$obs_stat_link
# maxmax_W=max(abs(bb))
# plot(bb,type="s",col="red",ylim=c(-5,5))
# abline(h=0)
# length(which((maxmax_What>=maxmax_W)==1))/30
##############################################################

##############################################################
# omnibus test
#####aft
# path1=30
# result_omni_aft=sample_path_omni(path1,beta_hat_ln_aft,std_hat_ln_aft,
# X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
# for(i in 1:path1){
# plot(result_omni_aft$dataset_What[[i]][,(n/2)],ylim=c(-3,3),type="s",
# col="grey");par(new=TRUE)
# }
# plot(result_omni_aft$dataset_W[,(n/2)],ylim=c(-3,3),type="s",col="red")
# for(i in 1:path1){
# plot(result_omni_aft$dataset_std.What[[i]][,(n/2)],ylim=c(-3,3),type="s",
# col="grey");par(new=TRUE)
# }
# plot(result_omni_aft$dataset_std.W[,(n/2)],ylim=c(-3,3),type="s",col="red")
# plotting_omni(result_omni_aft,"real",30)
# result_omni_aft$p_value
# result_omni_aft$std.p_value

#####cox
# path1=30
# result_omni_cox=sample_path_omni(path1,beta_hat_ln_cox,std_hat_ln_cox,
# X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
# for(i in 1:path1){
# plot(result_omni_aft$dataset_What[[i]][,(n/2)],ylim=c(-3,3),type="s",
# col="grey");par(new=TRUE)
# }
# plot(result_omni_cox$dataset_W[,(n/2)],ylim=c(-3,3),type="s",col="red")
# for(i in 1:path1){
# plot(result_omni_cox$dataset_std.What[[i]][,(n/2)],ylim=c(-3,3),type="s",
# col="grey");par(new=TRUE)
# }
# plot(result_omni_cox$dataset_std.W[,(n/2)],ylim=c(-3,3),type="s",col="red")
# plotting_omni(result_omni_cox,"real",30)
# result_omni_cox$p_value
# result_omni_cox$std.p_value

##############################################################
# functional form
#####aft
# path1=30
# result_form_aft=sample_path_form(path1,beta_hat_ln_aft,std_hat_ln_aft,
# X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
# for(i in 1:path1){
#   plot(result_form_aft$dataset_std.What[[i]][,1],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_form_aft$dataset_std.W[,1],ylim=c(-3,3),type="s",col="red")
# for(i in 1:path1){
#   plot(result_form_aft$dataset_std.What[[i]][,1],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_form_aft$dataset_std.W[,1],ylim=c(-3,3),type="s",col="red")
# plotting_form(result_form_aft,"real",30)
# result_form_aft$p_value
# result_form_aft$std.p_value

#####cox
# path1=30
# result_form_cox=sample_path_form(path1,beta_hat_ln_cox,std_hat_ln_cox,
# X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
# for(i in 1:path1){
#   plot(result_form_cox$dataset_What[[i]][,1],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_form_cox$dataset_W[,1],ylim=c(-3,3),type="s",col="red")
# for(i in 1:path1){
#   plot(result_form_cox$dataset_std.What[[i]][,1],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_form_cox$dataset_std.W[,1],ylim=c(-3,3),type="s",col="red")
# plotting_form(result_form_cox,"real",30)
# result_form_cox$p_value
# result_form_cox$std.p_value

##############################################################
# link function
#####aft
# path1=30
# result_link_aft=sample_path_link(path1,beta_hat_ln_aft,std_hat_ln_aft,
# X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
# for(i in 1:path1){
#   plot(result_link_aft$dataset_What[[i]],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_link_aft$dataset_W,ylim=c(-3,3),type="s",col="red")
# for(i in 1:path1){
#   plot(result_link_aft$dataset_std.What[[i]],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_link_aft$dataset_std.W,ylim=c(-3,3),type="s",col="red")
# plotting_form(result_form_aft,"real",30)
# result_form_aft$p_value
# result_form_aft$std.p_value

#####cox
# path1=30
# result_link_cox=sample_path_link(path1,beta_hat_ln_cox,std_hat_ln_cox,
# X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
# for(i in 1:path1){
#   plot(result_link_cox$dataset_What[[i]],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_link_cox$dataset_W,ylim=c(-3,3),type="s",col="red")
# for(i in 1:path1){
#   plot(result_link_cox$dataset_std.What[[i]],ylim=c(-3,3),type="s",
#   col="grey");par(new=TRUE)
# }
# plot(result_link_cox$dataset_std.W,ylim=c(-3,3),type="s",col="red")
# plotting_link(result_link_cox,"real",30)
# result_link_cox$p_value
# result_link_cox$std.p_value

##############################################################
##############################################################
##############################################################
############################미완성############################
##############################################################
##############################################################
##############################################################
W_aft=function(b,Time,Delta,Covari){
  #b=beta_hat_gg;Time=X_gg;Delta=D_gg;Covari=Z_gg
  #b=beta_hat_wb;Time=X_wb;Delta=D_wb;Covari=Z_wb
  
  e_i_beta=as.vector(log(Time)+Covari%*%b)
  
  order_resid=order(e_i_beta)
  
  Time=Time[order_resid]
  Covari=matrix(Covari[order_resid,],nrow=n)
  Delta=Delta[order_resid]
  e_i_beta=e_i_beta[order_resid]
  
  # weight function
  pi_ij_z=list(NA)
  for(j in 1:p){
    pi_ij_z[[j]]=Covari[,j]
  }
  
  N_i_t=list(NA)
  for(i in 1:n){
    N_i_t[[i]]=(e_i_beta>=e_i_beta[i])*Delta[i]
  }
  #N_i_t
  
  Y_i_t=list(NA)
  for(i in 1:n){
    Y_i_t[[i]]=(e_i_beta<=e_i_beta[i])*1
  }
  #Y_i_t
  
  N_d_t=Reduce('+',N_i_t)
  #N_d_t
  
  S_0_t=Reduce('+',Y_i_t)
  #S_0_t
  
  S_1_t=Reduce('+',mapply(function(x,y){x%*%t(y)},Y_i_t,as.list(data.frame(t(Covari))),SIMPLIFY=FALSE))
  #S_1_t
  
  E_t=S_1_t/S_0_t
  #E_t
  
  J_t=(S_0_t>0)*1
  #J_t
  
  dN_d_t=diff(c(0,N_d_t))
  #dN_d_t
  
  Lambdahat_0_t=cumsum((J_t/S_0_t)*dN_d_t)
  #Lambdahat_0_t
  
  dLambdahat_0_t=diff(c(0,Lambdahat_0_t))
  #dLambdahat_0_t
  
  Mhat_i_t=mapply("-", N_i_t,lapply(lapply(
    Y_i_t,"*",dLambdahat_0_t),cumsum), SIMPLIFY = FALSE)
  #Mhat_i_t
  
  pi_ij_z.Mhat_i_t=list(NA)
  for(j in 1:p){
    pi_ij_z.Mhat_i_t[[j]]=lapply(Mhat_i_t,function(x){pi_ij_z[[j]]%*%t(x)})
  }
  #pi_ij_z.Mhat_i_t
  # z by t matrix
  
  obs_stat=list(NA)
  for(j in 1:p){
    obs_stat[[j]]=Reduce('+',pi_ij_z.Mhat_i_t[[j]])/sqrt(n)
  }
  #obs_stat
  
  result=list(Time,Delta,Covari,e_i_beta,obs_stat)
  names(result)=c("Time","Delta","Covari","Resid","obs_stat")
  
  return(result)
}
#W_aft()

What_aft=function(b,std,Time,Delta,Covari,tol){}
#What_aft()

sample_path_aft=function(path,b,std,Time,Delta,Covari,tol){}
# result_aft=sample_path_aft(1,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,0.1)

plotting_aft=function(result,path){}

plotting_std.aft=function(result,path){}

# eps파일로 저장하는 코드
cairo_ps("afttest_link_std.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
afttestplot(afttest_link,standardization="standardized",path=path1)
dev.off()


dataset_cox=data.frame(X_ln_cox,D_ln_cox,Z_ln_cox)

path1=50
tol1=0.1

afttest_omni=afttest(Surv(X_ln_cox,D_ln_cox)~Z_ln_cox,dataset_cox,"omni",path1,tol1)
afttestplot(afttest_omni,standardization="unstandardized",path=path1)
afttestplot(afttest_omni,standardization="standardized",path=path1)

afttest_form=afttest(Surv(X_ln_cox,D_ln_cox)~Z_ln_cox,dataset_cox,"form",path1,tol1)
afttestplot(afttest_form,standardization="unstandardized",path=path1)
afttestplot(afttest_form,standardization="standardized",path=path1)

afttest_link=afttest(Surv(X_ln_cox,D_ln_cox)~Z_ln_cox,dataset_cox,"link",path1,tol1)
afttestplot(afttest_link,standardization="unstandardized",path=path1)
afttestplot(afttest_link,standardization="standardized",path=path1)





