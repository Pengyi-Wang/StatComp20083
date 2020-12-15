#' @title Benchmark2 R and Rcpp functions.
#' @name benchmarks2
#' @description The functions used in homework 2
#' @importFrom Rcpp evalCpp
#' @importFrom stats runif
#' @useDynLib StatComp20083
#' @examples
#' \dontrun{
#' HW2_1()
#' HW2_2()
#' }
NULL

#' @title The question 1 of homework 2
#' @return some numbers of the answer
#' @examples
#' \dontrun{
#' HW2_1()
#' }
#' @export
HW2_1= function() {
  pi=3.141592654
  set.seed(980904)
  n=1e4
  u=runif(n,min=0,max=pi/3)
  answer=mean(pi/3*sin(u))
  answer.sd=sd(pi/3*sin(u))/sqrt(n)
  print(c(answer,1/2,(answer-1/2)/(1/2),round(c(answer-1.96*answer.sd,answer+1.96*answer.sd),6)))
}

#' @title The question 2 of homework 2
#' @return some numbers of the answer
#' @examples
#' \dontrun{
#' HW2_2()
#' }
#' @export
HW2_2= function() {
  n=10000;
  x1=runif(n,min=0,max=1)
  x2=runif(n/2,min=0,max=1)
  T_V=0.9839
  print(c(var(exp(x1)),var(c((exp(x2)+exp(1-x2))/2)),(var(exp(x1))-var(c((exp(x2)+exp(1-x2))/2)))/var(exp(x1)),T_V ) )
}
