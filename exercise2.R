dataset_aft=data.frame(X_ln_aft,D_ln_aft,Z_ln_aft)

path1=10
tol1=0.1

afttest_omni=afttest(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,dataset_aft,"omni",path1,tol1)
afttestplot(afttest_omni,path=path1)
afttestplot(afttest_omni,standardization="standardized",path=path1)

afttest_fform=afttest(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,dataset_aft,"fform",path1,tol1)
afttestplot(afttest_fform,standardization="unstandardized",path=path1)
afttestplot(afttest_fform,standardization="standardized",path=path1)

afttest_linkf=afttest(Surv(X_ln_aft,D_ln_aft)~Z_ln_aft,dataset_aft,"linkf",path1,tol1)
afttestplot(afttest_linkf,standardization="unstandardized",path=path1)
afttestplot(afttest_linkf,standardization="standardized",path=path1)

cairo_ps("afttest_omni_aft.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
afttestplot(afttest_omni,standardization="unstandardized",path=path1)
dev.off()

cairo_ps("afttest_omni_std_aft.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
afttestplot(afttest_omni,standardization="standardized",path=path1)
dev.off()

cairo_ps("afttest_fform_cox.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
afttestplot(afttest_fform,standardization="unstandardized",path=path1)
dev.off()

cairo_ps("afttest_fform_std_cox.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
afttestplot(afttest_fform,standardization="standardized",path=path1)
dev.off()

cairo_ps("afttest_linkf_cox.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
afttestplot(afttest_linkf,standardization="unstandardized",path=path1)
dev.off()

cairo_ps("afttest_linkf_std_cox.eps",onefile=F,
         height=4,width=8, fallback_resolution = 600)
afttestplot(afttest_linkf,standardization="standardized",path=path1)
dev.off()

