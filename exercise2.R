aa=afttest_omni(250,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
aa$p_value
aa$std.p_value
# plotting_omni(aa,"rank",50)

aaa=afttest_omni(500,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
aaa$p_value
aaa$std.p_value
# plotting_omni(aaa,"rank",50)

aaaa=afttest_omni(1000,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
aaaa$p_value
aaaa$std.p_value
# plotting_omni(aaaa,"rank",50)

aaaaa=afttest_omni(1500,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
aaaaa$p_value
aaaaa$std.p_value
# plotting_omni(aaaaa,"rank",50)

bb=afttest_omni(250,beta_hat_ln_cox,std_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
bb$p_value
bb$std.p_value
# plotting_omni(bb,"rank",50)

bbb=afttest_omni(500,beta_hat_ln_cox,std_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
bbb$p_value
bbb$std.p_value
# plotting_omni(bbb,"rank",50)

bbbb=afttest_omni(1000,beta_hat_ln_cox,std_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
bbbb$p_value
bbbb$std.p_value
# plotting_omni(bbbb,"rank",50)

bbbbb=afttest_omni(1500,beta_hat_ln_cox,std_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
bbbbb$p_value
bbbbb$std.p_value
# plotting_omni(bbbbb,"rank",50)
# plot_ly(z=bbbbb$app_path[[1]],type="surface")
# plot_ly(z=bbbbb$obs_std.path,type="surface")

cc=afttest_form(250,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
cc$p_value
cc$std.p_value
# plotting_form(cc,"rank",50)

ccc=afttest_form(500,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
ccc$p_value
ccc$std.p_value
# plotting_form(ccc,"rank",50)

cccc=afttest_form(1000,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
cccc$p_value
cccc$std.p_value
# plotting_form(cccc,"rank",50)

ccccc=afttest_form(1500,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
ccccc$p_value
ccccc$std.p_value
# plotting_form(ccccc,"rank",50)


dd=afttest_form(250,beta_hat_ln_aft_f,std_hat_ln_aft_f,X_ln_aft_f,D_ln_aft_f,Z_ln_aft_f,given_tol)
dd$p_value
dd$std.p_value
# plotting_form(dd,"rank",50)

ee=afttest_link(250,beta_hat_ln_aft,std_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
ee$p_value
ee$std.p_value
# plotting_link(ee,"rank",50)

ff=afttest_link(250,beta_hat_ln_aft_f,std_hat_ln_aft_f,X_ln_aft_f,D_ln_aft_f,Z_ln_aft_f,given_tol)
ff$p_value
ff$std.p_value
# plotting_link(ff,"rank",50)

# dataset_aft=data.frame(X_ln_aft,D_ln_aft,Z_ln_aft)
# 
# path1=10
# tol1=0.1
# 
# afttest_omni=afttest(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,dataset_aft,"omni",path1,tol1)
# afttestplot(afttest_omni,path=path1)
# afttestplot(afttest_omni,standardization="standardized",path=path1)
# 
# afttest_fform=afttest(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,dataset_aft,"fform",path1,tol1)
# afttestplot(afttest_fform,standardization="unstandardized",path=path1)
# afttestplot(afttest_fform,standardization="standardized",path=path1)
# 
# afttest_linkf=afttest(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,dataset_aft,"linkf",path1,tol1)
# afttestplot(afttest_linkf,standardization="unstandardized",path=path1)
# afttestplot(afttest_linkf,standardization="standardized",path=path1)
# 
# cairo_ps("afttest_omni_aft.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# afttestplot(afttest_omni,standardization="unstandardized",path=path1)
# dev.off()
# 
# cairo_ps("afttest_omni_std_aft.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# afttestplot(afttest_omni,standardization="standardized",path=path1)
# dev.off()
# 
# cairo_ps("afttest_fform_cox.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# afttestplot(afttest_fform,standardization="unstandardized",path=path1)
# dev.off()
# 
# cairo_ps("afttest_fform_std_cox.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# afttestplot(afttest_fform,standardization="standardized",path=path1)
# dev.off()
# 
# cairo_ps("afttest_linkf_cox.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# afttestplot(afttest_linkf,standardization="unstandardized",path=path1)
# dev.off()
# 
# cairo_ps("afttest_linkf_std_cox.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# afttestplot(afttest_linkf,standardization="standardized",path=path1)
# dev.off()

