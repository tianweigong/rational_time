rgamma_int<-function(n,a,b,lb,up){
  re=rgamma(n,a,b)
  while (1) {
    tmp=rgamma(n,a,b)
    re=c(re,tmp) 
    re=re[re>lb &re<up]
    if (length(re)>=n){break}
  }
  re[1:n]
}

#for the delay model
# igamma<-function(dl,a,b){#interval gamma
#   pgamma(dl+0.5,a,b)-pgamma(dl-0.5,a,b)
# }

#for the rate model
myDecay<-function(s,r,t){
  if (length(t)==0){return(0)}
  dgamma(t,shape = s,rate=r)/dgamma((s-1)/r,shape = s,rate=r)
}

myRate<-function(sqc,par){
  l0=par[1]
  l1=par[2]

  li=rep(1,length(sqc))
  for (k in 1:length(sqc)){
    d=sqc[[k]]
    cauif=0
    if (length(d[["gen"]])){
      for (m in 1:length(d[["gen"]])){
        nam.al=paste("alpha.",names(d[["gen"]][m]),names(d[["rate"]]),sep="")
        nam.bt=paste("beta.",names(d[["gen"]][m]),names(d[["rate"]]),sep="")
        cauif=cauif+as.numeric(myDecay(par[nam.al],par[nam.bt],d[["gen"]])) 
      }
    }
    la=(l0+sum(l1*cauif))
    li[k]=dpois(d[["rate"]],la)
    # if(d[["rate"]]){break}
  }
  prod(li)
}
