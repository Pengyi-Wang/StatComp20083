#' @title Benchmark3 R and Rcpp functions.
#' @name benchmarks3
#' @description The functions used in homework 3
#' @examples
#' \dontrun{
#' HW3_2()
#' HW3_3()
#' HW3_4()
#' }
NULL

#' @title The question 2 of homework 3
#' @return some numbers
#' @examples
#' \dontrun{
#' HW3_2()
#' }
#' @export
HW3_2=function(){
  m=10000
  n=5
  N=50
  theta.hat=matrix(0,N,2)
  theta.hatj=numeric(5)
  g=function(x){(exp(-x)/(1+x^2))*(x>0)*(x<1)}
  for(i in 1:50){
    u=runif(m)    
    x=-log(1-u*(1-exp(-1)))
    fg=g(x)/(exp(-x)/(1-exp(-1)))
    theta.hat[i,1]=mean(fg)
    for(j in 1:n){
      u=runif(m/n,min=(j-1)/5,max=j/5)
      x=-log(exp(-(j-1)/5)-u*(exp(-(j-1)/5)-exp(-j/5)))
      fg=g(x)/(exp(-x)/(exp(-(j-1)/5)-exp(-j/5)))
      theta.hatj[j]=mean(fg)
    }
    theta.hat[i,2]=sum(theta.hatj)
  }
  apply(theta.hat,2,mean)
  apply(theta.hat,2,sd)
}

#' @title The question 3 of homework 3
#' @return some numbers
#' @examples
#' \dontrun{
#' HW3_3()
#' }
#' @export
HW3_3=function(){
  n=1000
  m=1000
  alpha=0.05
  theta1=theta2=numeric(n)
  for(i in 1:n){
    x=exp(rnorm(m))
    theta1[i]=mean(log(x))-qt(1-alpha/2,df=m-1)/sqrt(m)*var(log(x))
    theta2[i]=mean(log(x))+qt(1-alpha/2,df=m-1)/sqrt(m)*var(log(x))  
  }
  print(mean((theta1<0)*(0<theta2)))
}

#' @title The question 4 of homework 3
#' @return some numbers
#' @examples
#' \dontrun{
#' HW3_4()
#' }
#' @export
HW3_4=function(){
  n=10000
  m=20
  alpha=.05
  UCL=theta1=theta2=numeric(n)
  for(i in 1:n){
    x1=rnorm(m,mean=0,sd=2)
    UCL[i]=(m-1)*var(x1)/qchisq(alpha,df=m-1)
    x2=rchisq(m,2)
    theta1[i]=mean(x2)-qchisq(1-alpha/2,df=2*m)/m
    theta2[i]=mean(x2)+qchisq(alpha/2,df=2*m)/m 
  }
  print(c(mean((UCL>4)),mean((theta1<2)*(2<theta2))))
}
