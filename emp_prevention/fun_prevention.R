mydec_exp<-function(th,t){#decay
  if (length(t)==0){return(0)}
  exp(-t*th[1])
}

mydec_gam<-function(th,t){
  if (length(t)==0){return(0)}
  r=th[1]
  s=th[2]
  dgamma(t,shape = s,rate=r)/dgamma((s-1)/r,shape = s,rate=r)
}

mydec_gam_cum<-function(th,t){#concave, gamma cumulation function
  if (length(t)==0){return(0)}
  r=th[1]
  s=th[2]
  1-pgamma(t,shape = s,rate = r)
}


myLi<-function(d,v){
  l0=v[1]
  l1=v[2]
  p1=v[3]
  th1=v[c(4,5)]
  th2=v[c(6,7)]
  
  la=(l0+sum(l1*mydec_gen(th1,d[["gen"]])))*prod(1-p1*mydec_pre(th2,d[["pre"]]))
  
  prod((la^d[["rate"]])*(exp(-la))/factorial(d[["rate"]]))
  # sum(log((la^d[["rate"]])*(exp(-la))/factorial(d[["rate"]])))
}
