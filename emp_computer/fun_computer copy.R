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
myDecay<-function(th,t){
  if (length(t)==0){return(0)}
  s=th[1]
  r=th[2]
  dgamma(t,shape = s,rate=r)/dgamma((s-1)/r,shape = s,rate=r)
}

myRate<-function(sqc,par){
  l0=par[1]
  l1=par[2]
  th1=par[c(3,4)]

  li=rep(NA,length(sqc))
  for (k in 1:length(sqc)){
    d=sqc[[k]]
    la=(l0+sum(l1*myDecay(th1,d[["gen"]])))
    li[k]=dpois(d[["rate"]],la)
  }
  prod(li)
}
