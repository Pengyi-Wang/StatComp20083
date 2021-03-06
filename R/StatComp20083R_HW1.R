#' @title Benchmark1 R and Rcpp functions.
#' @name benchmarks1
#' @description The functions used in homework 1
#' @import bootstrap DAAG boot RANN energy Ball stats4 microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom RANN nn2
#' @importFrom boot boot
#' @importFrom energy eqdist.etest
#' @importFrom Ball bd.test
#' @importFrom stats4 mle
#' @importFrom stats runif rnorm rchisq rt rbeta qchisq qt qnorm qf pnorm pt var cov sd cor quantile lm uniroot t.test rpois
#' @importFrom graphics hist lines abline
#' @useDynLib StatComp20083
#' @examples
#' \dontrun{
#' HW1_1()
#' HW1_2()
#' HW1_4()
#' }
NULL

#' @title The question 1 of homework 1
#' @return nothing, but can get a picture.
#' @examples
#' \dontrun{
#' HW1_1()
#' }
#' @export
HW1_1=function(){
  a=b=2
  set.seed(980904)
  n=1000
  u=runif(n)
  x=b/(1-u)^(1/a)
  hist(x, prob = TRUE, breaks=seq(b,max(x)+1,1), main = expression(f(x)==a*b^a/x^(a+1)))
  y=seq(b,max(x)+1,1)
  lines(y,(a*b^a)/(y^(a+1)))
}

#' @title The question 2 of homework 1
#' @return nothing, but can get a picture.
#' @examples
#' \dontrun{
#' HW1_2()
#' }
#' @export
HW1_2=function(){
  n=10000
  set.seed(98)
  u1=runif(n,min=-1,max=1)
  set.seed(09)
  u2=runif(n,min=-1,max=1)
  set.seed(04)
  u3=runif(n,min=-1,max=1)
  x=1:n
  for(i in 1:n)
  {
    if(abs(u3[i])>=abs(u1[i])&abs(u3[i])>=abs(u2[i])){x[i]=u2[i]} else{x[i]=u3[i]}
  }
  hist(x, prob = TRUE, breaks=seq(-1,1,.01), main = expression(f(x)==3/4*(1-x^2)))
  y=seq(-1,1,.01)
  lines(y,3/4*(1-y^2))
}

#' @title The question 4 of homework 1
#' @return nothing, but can get a picture.
#' @examples
#' \dontrun{
#' HW1_4()
#' }
#' @export
HW1_4=function(){
  r=4
  b=2
  set.seed(980904)
  n=1000
  u=runif(n)
  y=b/(1-u)^(1/r)-b
  hist(y, prob = TRUE, breaks=seq(0,max(y)+1,.1), main = expression(f(y)==r*b^r/(b+y)^(r+1)))
  x=seq(0,max(y)+1,.1)
  lines(x,(r*b^r)/((b+x)^(r+1)))
}