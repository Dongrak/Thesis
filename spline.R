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

#########################################################################
t=seq(0,0.5,by=0.1)
y_form_100=c(0.0004,0.002,0.020,0.072,0.230,0.361)
y_form_250=c(0.0028,0.033,0.223,0.552,0.791,0.920)

Figure_form=ggplot()+geom_point(aes(t,y_form_100),colour="red")+geom_point(aes(t,y_form_250),colour="blue")+
  geom_line(aes(x,y),colour="red",data=data.frame(spline(t,y_form_100,n=length(y_form_100)*10,method="hyman")),linetype="dashed")+
  geom_line(aes(x,y),colour="blue",data=data.frame(spline(t,y_form_250,n=length(y_form_100)*10,method="hyman")))+
  ylab("Rejection Ratio")+xlab(bquote(gamma))+
  theme(plot.title=element_text(hjust=0.5))
  
cairo_ps("afttest_form_gamma.eps",onefile=F,
         height=3,width=8, fallback_resolution = 600)
Figure_form
dev.off()

#########################################################################
t=seq(0,0.5,by=0.1)
y_link_100=c(0,0.116,0.141,0.154,0.180,0.206)
y_link_250=c(0,0.746,0.764,0.777,0.794,0.798)

Figure_link=ggplot()+geom_point(aes(t,y_link_100),colour="red")+geom_point(aes(t,y_link_250),colour="blue")+
  geom_line(aes(x,y),colour="red",data=data.frame(spline(t,y_link_100,n=length(y_form_100)*10,method="hyman")),linetype="dashed")+
  geom_line(aes(x,y),colour="blue",data=data.frame(spline(t,y_link_250,n=length(y_form_100)*10,method="hyman")))+
  ylab("Rejection Ratio")+xlab(bquote(gamma))+
  theme(plot.title=element_text(hjust=0.5));Figure_link

cairo_ps("afttest_link_gamma.eps",onefile=F,
         height=3,width=8, fallback_resolution = 600)
Figure_link
dev.off()

#########################################################################
t=seq(0,0.5,by=0.1)
y_omni_omni_100=c(0.0026,0.007,0.034,0.066,0.103,0.118)
y_omni_omni_250=c(0.0072,0.093,0.324,0.538,0.664,0.708)

Figure_omni=ggplot()+geom_point(aes(t,y_omni_omni_100),colour="red")+geom_point(aes(t,y_omni_omni_250),colour="blue")+
  geom_line(aes(x,y),colour="red",data=data.frame(spline(t,y_omni_omni_100,n=length(y_form_100)*10,method="hyman")),linetype="dashed")+
  geom_line(aes(x,y),colour="blue",data=data.frame(spline(t,y_omni_omni_250,n=length(y_form_100)*10,method="hyman")))+
  ylab("Rejection Ratio")+xlab(bquote(gamma))+
  theme(plot.title=element_text(hjust=0.5));Figure_omni

cairo_ps("afttest_omni_gamma.eps",onefile=F,
         height=3,width=8, fallback_resolution = 600)
Figure_omni
dev.off()