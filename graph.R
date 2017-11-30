result=aa

quant=c(1,round(n/8),round(n/4),round(3*n/8),round(n/2),round(5*n/8),round(3*n/4),round(7*n/8),n)

Figure=list(NA)
for(k in 1:9){
  
  Q=quant[k]
  
  dataset_What=data.frame()
  
  for (i in 1:50){
    group=i
    A=result$app_path[[i]][,Q]
    AA=data.frame(group,t_i=xaxix,What=A)
    dataset_What=rbind(dataset_What,AA)
  }
  #dataset_What
  
  dataset_W=data.frame(group,t_i=xaxix,W=result$obs_path[,Q])
  #dataset_W
  
  Figure_W=
    ggplot()+
    geom_step(data=dataset_What,aes(x=t_i,y=What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
    geom_step(data=dataset_W,aes(x=t_i,y=W),colour="tomato",lwd=0.25)+
    theme_minimal()
  #Figure_W
  Figure[[k]]=Figure_W
}

# Figure[[1]]

grid.arrange(Figure[[1]],Figure[[2]],Figure[[3]],Figure[[4]],Figure[[5]],
             Figure[[6]],Figure[[7]],Figure[[8]],Figure[[9]],nrow=3)
