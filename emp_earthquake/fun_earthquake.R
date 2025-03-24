#event-based
s2<-function(v){
  sprintf("%02d",v)
}

myG1<-function(df.path,df.delay,df.prob,del.nam,del.tim,del.typ,del.pro,df.miss,t_c,t_end,br,c_p,c_al,c_bt){ #baserate and cause
  apply(cbind(1:samp),1,function(x){
    
    del.pro[del.typ]=dgamma(del.tim[del.typ],c_al[x],c_bt[x])*c_p[x]
    del.pro[!del.typ]=dgamma(del.tim[!del.typ],1,br[x])
    
    for (k in del.nam){
      df.prob[which(df.delay==k)]=del.pro[k]
    }
    
    for (k in 1:length(t_c)){
      df.miss[which(df.miss[,k]==0),k]= (1-c_p[x])+c_p[x]*(1-pgamma(t_end-t_c[k],c_al[x],c_bt[x]))
    }
    sum(rowProds(df.prob)*rowProds(df.miss))
  })
}

myEventBased<-function(sqc,br,c_p,c_al,c_bt){
  t_c=sqc$time[sqc$obj!="E"] %>% sort() 
  names(t_c)=paste("C",s2(c(1:length(t_c))),sep="")
  
  t_e=sqc$time[sqc$obj=="E"] %>% sort()
  names(t_e)=paste("E",s2(c(1:length(t_e))),sep="")
  t_end=max(sqc$time)
  
  lis.cau=list()
  for (k in 1:length(t_e)){
    cau=t_c[t_c<t_e[k]]
    cauTr=which(cau>t_e[k]-8) #trim
    if (length(cauTr)){
      lis.cau[[k]]=c("B",paste("C",s2(cauTr),sep=""))
    }else{
      lis.cau[[k]]=c("B")
    }
  }
  
  df=expand.grid(lis.cau,stringsAsFactors = F)
  idx=df %>% apply(1, function(x){ y=x[x!="B"]
  if (!length(y) || length(unique(y))==length(y)){T}else{F}})
  df.path=df[idx,]

  df.delay=df.path %>% apply(1, function(x){
    idx_e=which(x=="B")
    idx_c=which(x!="B")
    if (length(idx_e)){
      x[idx_e]=paste("E",s2(c(0,idx_e[-length(idx_e)])),"E",s2(idx_e),sep="")
    }
    if (length(idx_c)){
      x[idx_c]=paste(x[idx_c],"E",s2(idx_c),sep="")
    }
    x
  }) %>% t() %>%as.matrix()
  
  df.prob=matrix(NA,nrow = nrow(df.delay),ncol=ncol(df.delay))
  
  del.nam=unique(as.character(df.delay))
  del.tim=c()
  del.pro=rep(NA,length(del.nam))
  names(del.pro)=del.nam
  del.typ=substr(del.nam,1,1)=="C"
  t_all=c(t_c,t_e,"E00"=0)
  for (k in del.nam){
    t_d=round(t_all[substr(k,4,6)]-t_all[substr(k,1,3)],2)
    del.tim=c(del.tim,as.numeric(t_d))
  }
  
  df.miss=matrix(1,nrow=nrow(df.path),ncol =length(t_c))
  for (k in 1:nrow(df.miss)){
    tmp=setdiff(names(t_c),df.path[k,])
    if (length(tmp)){df.miss[k,as.numeric(substr(tmp,2,3))]=0}
  }
  
  g1=myG1(df.path,df.delay,df.prob,del.nam,del.tim,del.typ,del.pro,df.miss,t_c,t_end,br,c_p,c_al,c_bt)
  
  # ix1=sort(g1, index.return=TRUE,decreasing=TRUE)$ix[1:100]
  # g1_par=list("br"=br,"c_p"=c_p,"c_m"=c_m,"c_v"=c_v)
  list("logli"=log(mean(g1)),"G1"=g1)
}

#rate-based
mydec<-function(t,s,r){
  dgamma(t,shape = s,rate=r)/dgamma((s-1)/r,shape = s,rate=r)
}

myG1_poi<-function(dl,df.cau,br,c_p,c_al,c_bt){ #baserate and cause
  x=rep(dl,samp) %>% matrix(nrow=length(dl),ncol=samp)
  y=apply(cbind(1:samp),1,function(x){c_p[x]*colSums(mydec(df.cau+0.5,c_al[x],c_bt[x]))+br[x]})
  dpois(x,y) %>% apply(.,2,prod)
}


myData<-function(dl,mycut){
  table(cut(dl, breaks = mycut)) %>% as.numeric()
}

myRateBased<-function(sqc,br,c_p,c_al,c_bt){
  t_c=sqc$time[sqc$obj!="E"] %>% sort() 
  names(t_c)=paste("C",s2(c(1:length(t_c))),sep="")
  
  t_e=sqc$time[sqc$obj=="E"] %>% sort()
  names(t_e)=paste("E",s2(c(1:length(t_e))),sep="")
  t_end=ceiling(max(sqc$time))
  mycut=c(0:t_end)+0.00001
  
  dl=t_e %>% myData(.,mycut)
  
  df.cau=matrix(NA,nrow=length(t_c),ncol=length(mycut)-1)
  for (k in 1:length(t_c)){
    df.cau[k,]=mycut[-length(mycut)]-t_c[k]
  }
  
  g1=myG1_poi(dl,df.cau,br,c_p,c_al,c_bt)
  
  list("logli"=log(mean(g1)),"G1"=g1)
}