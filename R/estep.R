estep=function(x,theta_t,mu_t,pi_t,t){
  jt.df=NULL
  ldf=NULL
  K=length(pi_t)
  for(k in 1:K){
    ldf=rbind(ldf,ldnb(x,theta_t[k,],mu_t[k,]))
  }
  
  for(k in 1:K){
    jt.df=c(jt.df,pi_t[k]*exp(sum((ldf[k,]-mean(ldf))/t)))
  }
  
  pi=jt.df/sum(jt.df)
  return(pi)
}