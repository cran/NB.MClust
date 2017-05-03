#' dnb, ldnb Functions
#'
#' These functions allow you to compute (log-)density of generalized Negative Binomial distribution. 
#' @name nb.density
#' @param x A positive numeric scalor or vector. Decimals and integers are both allowed.   
#' @param  theta Value of dispersion.
#' @param  mu Value of mean.
#' @return \item{dnb}{Density of generalized Negative Binomial} 
#'        \item{ldnb}{Log-density of generalized Negative Binomial} 
NULL
 
#' @rdname nb.density
#' @examples
#' ldnb(x=10.4,theta=3.2,mu=5)
#' @export 
   ldnb=function(x,theta,mu){lgamma(x+theta)-(lgamma(theta)+lfactorial(x))+x*(log(mu)-log(theta+mu))+theta*(log(theta)-log(theta+mu))} 



#' @rdname nb.density
#' @examples
#' dnb(x=10.4,theta=3.2,mu=5)
#' @export 
   dnb=function(x,theta,mu){exp(ldnb(x,theta,mu))}





