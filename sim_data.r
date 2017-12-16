#-------------------------------------------------------------
#---------------------------OMNIBUS---------------------------
#-------------------------------------------------------------

#-----------------------------AFT-----------------------------
# aa=afttest_omni(200,beta_hat_ln_aft,se_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
aa$p_value
aa$std.p_value
plotting_omni(aa,50,"unstd")
plotting_omni(aa,50,"std")

# cairo_ps("afttest_omni_aft_unstd.eps",onefile=F,
#          height=5,width=8, fallback_resolution = 600)
# plotting_omni(aa,50,"unstd")
# dev.off()
# 
# cairo_ps("afttest_omni_aft_std.eps",onefile=F,
#          height=5,width=8, fallback_resolution = 600)
# plotting_omni(aa,50,"std")
# dev.off()

aaa=afttest_omni(200,beta_hat_ln_aft_f,se_hat_ln_aft_f,X_ln_aft_f,D_ln_aft_f,Z_ln_aft_f,given_tol)
aaa$p_value
aaa$std.p_value
plotting_omni(aaa,50,"unstd")
plotting_omni(aaa,50,"std")

# cairo_ps("afttest_omni_aft_unstd_f.eps",onefile=F,
#          height=5,width=8, fallback_resolution = 600)
# plotting_omni(aaa,50,"unstd")
# dev.off()
# 
# cairo_ps("afttest_omni_aft_std_f.eps",onefile=F,
#          height=5,width=8, fallback_resolution = 600)
# plotting_omni(aaa,50,"std")
# dev.off()

#-----------------------------COX-----------------------------
# bb=afttest_omni(200,beta_hat_ln_cox,se_hat_ln_cox,X_ln_cox,D_ln_cox,Z_ln_cox,given_tol)
bb$p_value
bb$std.p_value
plotting_omni(bb,50,"unstd")
plotting_omni(bb,50,"std")

# cairo_ps("afttest_omni_cox_unstd.eps",onefile=F,
#          height=5,width=8, fallback_resolution = 600)
# plotting_omni(bb,50,"unstd")
# dev.off()
# 
# cairo_ps("afttest_omni_cox_std.eps",onefile=F,
#          height=5,width=8, fallback_resolution = 600)
# plotting_omni(bb,50,"std")
# dev.off()

#-------------------------------------------------------------
#---------------------------FTNFORM---------------------------
#-------------------------------------------------------------

#------------------------------1------------------------------
# cc=afttest_form(200,beta_hat_ln_aft,se_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
cc$p_value
plotting_form(cc,50)

# cairo_ps("afttest_form_aft.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# plotting_form(cc,50)
# dev.off()

#------------------------------2------------------------------
# dd=afttest_form(200,beta_hat_ln_aft_f,se_hat_ln_aft_f,X_ln_aft_f,D_ln_aft_f,Z_ln_aft_f,given_tol)
dd$p_value
plotting_form(dd,50)

# cairo_ps("afttest_form_aft_f.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# plotting_form(dd,50)
# dev.off()

#-------------------------------------------------------------
#---------------------------LINKFTN---------------------------
#-------------------------------------------------------------

#------------------------------1------------------------------
# ee=afttest_link(200,beta_hat_ln_aft,se_hat_ln_aft,X_ln_aft,D_ln_aft,Z_ln_aft,given_tol)
ee$p_value
plotting_link(ee,50)

# cairo_ps("afttest_link_aft.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# plotting_link(ee,50)
# dev.off()

#------------------------------2------------------------------
# ff=afttest_link(200,beta_hat_ln_aft_l,se_hat_ln_aft_l,X_ln_aft_l,D_ln_aft_l,Z_ln_aft_l,given_tol)
ff$p_value
plotting_link(ff,50)

# cairo_ps("afttest_link_aft_l.eps",onefile=F,
#          height=4,width=8, fallback_resolution = 600)
# plotting_link(ff,50)
# dev.off()