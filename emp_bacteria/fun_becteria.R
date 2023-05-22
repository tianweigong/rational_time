#rate model
myDecay<-function(th1,t){
  s=th1[1]
  r=th1[2]
  dgamma(t,shape = s,rate=r)/dgamma((s-1)/r,shape = s,rate=r)
}

myLam_g<-function(d,l0,l1,th1){
  la=l0+l1*myDecay(th1,c(1:length(d)))
  
  la=pmin(n_num-cumsum(c(0,d[1:(length(d)-1)])), la)
  
  cbind(d,la) %>% apply(1, function(x){ dpois(x[1],x[2]) }) %>%prod()
}

myLam_n<-function(d,l0){
  la=l0
  
  la=pmin(n_num-cumsum(c(0,d[1:(length(d)-1)])), la)
  
  cbind(d,la) %>% apply(1, function(x){ dpois(x[1],x[2]) }) %>%prod()
}

myLam_p<-function(d,l0,p1,th1){
  la=l0*(1-p1*myDecay(th1,c(1:length(d))))
  
  la=pmin(n_num-cumsum(c(0,d[1:(length(d)-1)])), la)
  
  cbind(d,la) %>% apply(1, function(x){ dpois(x[1],x[2]) }) %>%prod()
}

myType<-function(v){
  
  l0=v[1]
  l1=v[2]
  p1=v[3]
  th1=v[c(4,5)]
  
  
  li_e_g=myLam_g(eg,l0,l1,th1)
  
  li_e_p=myLam_p(eg,l0,p1,th1)#preventative
  
  li_e_n=myLam_n(eg,l0)#non-causal
  
  li_c=myLam_n(cg,l0)
  
  c(li_e_g*li_c,li_e_p*li_c,li_e_n*li_c)
}