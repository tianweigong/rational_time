library(dplyr)

sti=data.frame()
for (k in c("Con1","Con2","Con3","Con4")){
  fl=list.files(path=paste('./',k,sep = ""),pattern = "\\inp$")
  for (j in fl){
    tmp=read.table(paste(k,j,sep = "/"),header = T) %>%
      mutate(st=trial.time.event,
             order=sapply(strsplit(st,","),"[",1),
             obj=sapply(strsplit(st,","),"[",3) %>% 
               factor(levels=c(0,1,2,3),labels = c("E","A","B","C")),
             time=sapply(strsplit(st,","),"[",2),
      )%>%
      select(order,obj,time) %>%
      na.omit()%>%
      mutate(Delay=c("Con1"="short","Con2"="short","Con3"="long","Con4"="long")[k],
             PoIE=c("Con1"="high","Con2"="low","Con3"="high","Con4"="low")[k],
             condition=paste(Delay,PoIE,sep = "_"),
             trial_id= paste(condition,sapply(strsplit(j,"[.]"),"[",1),sep = "_")
             )
    sti=rbind(sti,tmp)
  }
}



save(sti,file = "data_earthquake.Rda")

