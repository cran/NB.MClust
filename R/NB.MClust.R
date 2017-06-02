#' NB.MClust Function
#'
#' This function performs model-based clustering on positive integer or continuous data that follow Generalized Negative Binomial distribution. 
#' 
#' @param Count Data matrix of discrete counts.This function groups rows of the data matrix. 
#' @param  K Number of clusters or components specified. It can be a positive integer or a vector of positive integer.
#' @param  ini.shift.mu Initial value in EM algorithm for the shift between clusters in mean.
#' @param  ini.shift.theta Initial value in EM algorithm for the shift between clusters in dispersion. 
#' @param  tau0 Initial value of anealing rates in EM Algorithm. Default and suggested value is 10.
#' @param  rate Stochastic decreasing speed for anealing rate. Default and suggested value is 0.9
#' @param  bic  Whether Bayesian Information should be computed when K is an integer. BIC is forced to be TRUE when K is a vector. 
#' @param  iteration Maximum number of iterations in EM Algorithm, default at 50.  
#' @keywords NB.MClust
#' @return  \item{parameters}{Estimated parameters}
#' 
#'           \item{$prior}{Prior probability that a sample belongs to each cluster}               
#'           \item{$mu}{Mean of each cluster} 
#'           \item{$theta}{Dispersion of each cluster} 
#'           \item{$posterior}{Posterior probability that a sample belongs to each cluster}
#'          
#'          \item{cluster}{Estimated cluster assignment} 
#'          
#'          \item{BIC}{Value of Bayesian Information} 
#'          
#'          \item{K}{Optional or estimated number of clusters, if input K is a vector}
#'  
#' @import MASS
#'         utils  
#' @examples 
#' # Example:
#' 
#' data("Simulated_Count") # A 50x100 integer data frame.
#' 
#' m1=NB.MClust(Simulated_Count,K=2:5)
#' cluster=m1$cluster #Estimated cluster assignment
#' k_hat=m1$K  #Estimated optimal K
#' 
#' @export 


NB.MClust=function(Count,K,ini.shift.mu=0.01,ini.shift.theta=0.01,tau0=10,rate=0.9,bic=TRUE,iteration=100){
  
  if(any(is.na(Count))) stop('Invalid values in count matrix: NA')
  if(any(Count<0)) stop('Invalid values in count matrix: Negative')
  if(any(is.na(K))) stop('K must be specified')
  if(any(K<0)) stop('K must be positive integer(s)')
  Count=as.matrix(Count)
if (length(K)==1){  
  NBMB(Count,K,ini.shift.mu,ini.shift.theta,tau0,rate,bic,iteration)
}else{
  bic=TRUE
  m=vector(mode='list',length=length(K))
  all.bic=NULL
  for(k in K){
    m[[k]]=NBMB(Count,k,ini.shift.mu,ini.shift.theta,tau0,rate,bic=T,iteration)
    all.bic=c(all.bic,m[[k]]$BIC)
  }
  opt.K=K[all.bic==min(all.bic,na.rm=T)]
  opt.m=m[[opt.K]]
  return(opt.m)
 }
}






NBMB=function(Count,K,ini.shift.mu,ini.shift.theta,tau0,rate,bic,iteration){

if(!requireNamespace('MASS')) {
  install.packages('MASS')  
  if(!requireNamespace('MASS',character.only = TRUE)) stop("Package 'MASS' not found")}
  
ini.prior=1/K
Count<<-Count
N=nrow(Count)
G=ncol(Count)
mle=apply(Count,2,function(x){
  m0=glm.nb(x~1)
  return(c(theta0=m0$theta,mu0=m0$coefficients))
}
)

mu_tmp=NULL
for(k in 1:K){
  mu_tmp=rbind(mu_tmp,(1+ini.shift.mu*(k-1))*exp(mle[2,]))}

  pi_tmp=rep(1/K,K)


theta.clean=ifelse(mle[1,]>1000,mean(mle[1,][mle[1,]<=1000]),mle[1,])

theta_tmp=NULL
for(k in 1:K){
  theta_tmp=rbind(theta_tmp,(1+ini.shift.theta*(k-1))*theta.clean)}


post=apply(Count,1,function(x)estep(x,theta_tmp,mu_tmp,pi_tmp,tau0))
pi_next=rowMeans(post,na.rm = T)
post=apply(post,2,function(x){ifelse(is.na(x),pi_next,x)})

mu_next=(post%*%Count)/(N*pi_next)

itr=0
while(max(abs(mu_tmp-mu_next))>0.1 & mean(apply(post,2,min))>1e-3 ){
  
  mu_tmp=mu_next
  pi_tmp=pi_next
  tau0=rate*tau0
  
  post=apply(Count,1,function(x)estep(x,theta_tmp,mu_tmp,pi_tmp,tau0))
  pi_next=rowMeans(post,na.rm =T)
  post=apply(post,2,function(x){ifelse(is.na(x),1/K,x)})
  
  mu_next=(post%*%Count)/(N*pi_next)
  itr=itr+1
  if(itr>iteration)break
}

mu_tmp=mu_next
pi_tmp=pi_next
tau0=rate*tau0

post=apply(Count,1,function(x)estep(x,theta_tmp,mu_tmp,pi_tmp,tau0))
pi_next=rowMeans(post,na.rm =T)
if(sum(is.na(post))>0){
  naid=which(colSums(is.na(post))>0)
  post[,naid]=pi_next
  
}
cluster=apply(post,2,function(x){(1:K)[max(x)==x]})
BIC=NA
if(bic==T)
{llk=NULL
for(k in 1:K){
  llk=rbind(llk,apply(Count,1,function(x){pi_next[k]*sum(ldnb(x,theta_tmp[k,],mu_tmp[k,]))}))
}
llk=ifelse(is.infinite(llk),-700,llk)
BIC=-2*sum(llk,na.rm=T)/N+K*log(N)/N
}
return(list(parameters=list(prior=pi_next,mu=mu_tmp,theta=theta_tmp,posterior=post),cluster=cluster,BIC=BIC,K=K))

}