s2<-function(v){sprintf("%02d",v)}

myStiGen<-function(k_b,k_c,w_c,mu,sigm){
  while(1){
    sti.lis=c()
    p=1
    t_c_all=runif(max(k_c),0,t_end)
    t_be_all=runif(max(k_b),0,t_end)
    for (k_b1 in k_b){
      t_be=t_be_all[1:k_b1]%>% sort()
      for (k_c1 in k_c){
        t_c=t_c_all[1:k_c1] %>% sort()
        names(t_c)=paste("C",s2(c(1:length(t_c))),sep="")
        for (w_c1 in w_c){
          wclis=sample(c(T,F),size=k_c1,prob=c(w_c1,1-w_c1),replace = T)
          for (mu1 in mu){
            for (sigm1 in sigm){
              t_ce=t_c+runif(length(t_c),mu1-sigm1,mu1+sigm1)
              t_e=c(t_ce[wclis],t_be) %>% sort()
              names(t_e)=paste("E",s2(c(1:length(t_e))),sep="")
              
              sti.lis[[p]]=list("k_b"=k_b1,"k_c"=k_c1,"w_c"=w_c1,"mu"=mu1,"sigm"=sigm1,
                                "t_c"=t_c,"t_e"=t_e,"t_ce"=t_ce,"t_be"=t_be)
              p=p+1
            }
          }
        }
      }
    }
    if ( max(unlist(sapply(sti.lis,"[[","t_e")))<t_end){break}
  }
  sti.lis
}

myRun<-function(lis){
  li.f=paste("sim/",list.files("sim"),sep="")
  if (paste(lis[["nam"]],lis[["simCode"]],".Rda",sep="") %in% li.f){
    return()
  }
  sti=lis[["sti"]]
  if ("pri" %in% names(lis)){
    for (k in 1:length(sti)){sti[[k]][["pri"]]=lis[["pri"]]}
  }
  mm=lapply(sti, mySim)
  save(mm,file = paste(lis[["nam"]],lis[["simCode"]],".Rda",sep=""))
}
