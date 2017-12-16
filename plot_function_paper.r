plotting_omni=function(result,path,std){
  
  result_Time=result$Time
  result_log_Time=log(result_Time)
  n=length(result_Time)
  xaxix=(1:n)
  # if (xaxix=="rank"){xaxix=(1:n)[order(result_Time)]}
  # else if (xaxix=="real"){xaxix=result_Time}
  # else (xaxix=result_Time)
  if (std=="std"){
    Figure=list(NA)
    
    for(k in 1:9){
      
      quant=round(quantile(1:n,c(0,0.1,0.25,0.4,0.5,0.6,0.75,0.9,1)))
      
      Q=quant[k]
      
      dataset_std.What=data.frame()
      
      for (i in 1:path){
        group=i
        A=result$app_std.path[[i]][,Q]
        AA=data.frame(group,e_i=xaxix,std.What=A)
        dataset_std.What=rbind(dataset_std.What,AA)
      }
      #dataset_std.What
      
      dataset_std.W=data.frame(group,e_i=xaxix,std.W=result$obs_std.path[,Q])
      #dataset_std.W
      
      Figure_std.W=
        ggplot()+
        geom_step(data=dataset_std.What,aes(x=e_i,y=std.What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
        geom_step(data=dataset_std.W,aes(x=e_i,y=std.W),colour="black",lwd=1.5)+
        ylab("Test Statistic")+xlab("Time Transformed Residuals")+
        ggtitle(paste("Quantile of z",names(quant)[k]))+theme(plot.title=element_text(hjust=0.5))
      #Figure_W
      Figure[[k]]=Figure_std.W
    }
    
    # Figure[[1]]
    
    lay=rbind(c(1,1,1,1),c(1,1,1,1),c(2,3,4,5),c(6,7,8,9))
    
    return(grid.arrange(Figure[[5]],Figure[[1]],Figure[[2]],Figure[[3]],Figure[[4]],
                        Figure[[6]],Figure[[7]],Figure[[8]],Figure[[9]],layout_matrix=lay))
  }
  else{
    Figure=list(NA)
    
    for(k in 1:9){
      
      quant=round(quantile(1:n,c(0,0.1,0.25,0.4,0.5,0.6,0.75,0.9,1)))
      
      Q=quant[k]
      
      dataset_What=data.frame()
      
      for (i in 1:path){
        group=i
        A=result$app_path[[i]][,Q]
        AA=data.frame(group,e_i=xaxix,What=A)
        dataset_What=rbind(dataset_What,AA)
      }
      #dataset_What
      
      dataset_W=data.frame(group,e_i=xaxix,W=result$obs_path[,Q])
      #dataset_W
      
      Figure_W=
        ggplot()+
        geom_step(data=dataset_What,aes(x=e_i,y=What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
        geom_step(data=dataset_W,aes(x=e_i,y=W),colour="black",lwd=1.5)+
        ylab("Test Statistic")+xlab("Time Transformed Residuals")+
        ggtitle(paste("Quantile of z",names(quant)[k]))+theme(plot.title=element_text(hjust=0.5))
      #Figure_W
      Figure[[k]]=Figure_W
    }
    
    # Figure[[1]]
    
    lay=rbind(c(1,1,1,1),c(1,1,1,1),c(2,3,4,5),c(6,7,8,9))
    
    return(grid.arrange(Figure[[5]],Figure[[1]],Figure[[2]],Figure[[3]],Figure[[4]],
                        Figure[[6]],Figure[[7]],Figure[[8]],Figure[[9]],layout_matrix=lay))
  }
}

plotting_form=function(result,path){
  
  # result_Covari=result$Covari
  n=length(result$Covari)/length(result$Beta)
  xaxix=(1:n)
  # if (xaxix=="rank"){xaxix=(1:n)[order(result_Covari)]}
  # else {xaxix=result_Covari}
  
  dataset_What=data.frame()
  
  for (i in 1:path){
    group=i
    A=result$app_path[[i]]
    AA=data.frame(group,e_i=xaxix,What=A)
    dataset_What=rbind(dataset_What,AA)
  }
  #dataset_What
  
  dataset_W=data.frame(group,e_i=xaxix,W=result$obs_path)
  #dataset_W
  
  Figure1_W=
    ggplot()+
    geom_step(data=dataset_What,aes(x=e_i,y=What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
    geom_step(data=dataset_W,aes(x=e_i,y=W),colour="black",lwd=1.5)+
    ylab("Test Statistic")+xlab("Time Transformed Residuals")+
    theme_minimal()
  #Figure1_W
  
  return(Figure1_W)
  
  # dataset_std.What=data.frame()
  # 
  # for (i in 1:path){
  #   group=i
  #   A=result$app_std.path[[i]]
  #   AA=data.frame(group,e_i=xaxix,std.What=A)
  #   dataset_std.What=rbind(dataset_std.What,AA)
  # }
  # #dataset_std.What
  # 
  # dataset_std.W=data.frame(group,e_i=xaxix,std.W=result$obs_std.path)
  # #dataset_std.W
  # 
  # Figure1_std.W=
  #   ggplot()+
  #   geom_step(data=dataset_std.What,aes(x=e_i,y=std.What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
  #   geom_step(data=dataset_std.W,aes(x=e_i,y=std.W),colour="black",lwd=1.5)+
  #   theme_minimal()
  # #Figure1_std.W
  # 
  # return(grid.arrange(Figure1_W,Figure1_std.W,nrow=2))
}

plotting_link=function(result,path){
  
  # result_Covari=result$Covari
  n=length(result$Covari)/length(result$Beta)
  xaxix=(1:n)
  # if (xaxix=="rank"){xaxix=(1:n)[order(result_Covari)]}
  # else {xaxix=result_Covari}
  
  dataset_What=data.frame()
  
  for (i in 1:path){
    group=i
    A=result$app_path[[i]]
    AA=data.frame(group,e_i=xaxix,What=A)
    dataset_What=rbind(dataset_What,AA)
  }
  #dataset_What
  
  dataset_W=data.frame(group,e_i=xaxix,W=result$obs_path)
  #dataset_W
  
  Figure1_W=
    ggplot()+
    geom_step(data=dataset_What,aes(x=e_i,y=What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
    geom_step(data=dataset_W,aes(x=e_i,y=W),colour="black",lwd=1.5)+
    ylab("Test Statistic")+xlab("Time Transformed Residuals")+
    theme_minimal()
  #Figure1_W
  
  return(Figure1_W)
  
  # dataset_std.What=data.frame()
  # 
  # for (i in 1:path){
  #   group=i
  #   A=result$app_std.path[[i]]
  #   AA=data.frame(group,e_i=xaxix,std.What=A)
  #   dataset_std.What=rbind(dataset_std.What,AA)
  # }
  # #dataset_std.What
  # 
  # dataset_std.W=data.frame(group,e_i=xaxix,std.W=result$obs_std.path)
  # #dataset_std.W
  # 
  # Figure1_std.W=
  #   ggplot()+
  #   geom_step(data=dataset_std.What,aes(x=e_i,y=std.What,group=group),colour=adjustcolor("#737c8c", alpha=0.8),alpha=0.5)+
  #   geom_step(data=dataset_std.W,aes(x=e_i,y=std.W),colour="black",lwd=1.5)+
  #   theme_minimal()
  # #Figure1_std.W
  
  # return(grid.arrange(Figure1_W,Figure1_std.W,nrow=2))
}
