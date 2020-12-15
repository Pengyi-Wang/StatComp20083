#' @title Benchmark4 R and Rcpp functions.
#' @name benchmarks4
#' @description The functions used in homework 4
#' @importFrom Rcpp evalCpp 
#' @importFrom stats rnorm runif qnorm rbeta rt
#' @useDynLib StatComp20083
#' @examples
#' \dontrun{
#' HW4_1_1()
#' HW4_1_2()
#' HW4_2_1()
#' HW4_2_2()
#' HW4_2_3()
#' HW4_3_1()
#' HW4_3_2()
#' }
NULL

#' @title the question 1.1 of homework (Beta)
#' @return a data.frame
#' @examples
#' \dontrun{
#' HW4_1_1()
#' }
#' @export
HW4_1_1=function() {
  alpha=.05
  beta=1:5
  n=c(10,20,30,50,100)
  m=1000
  result=matrix(0,length(n),length(beta))
  f=function(alpha,beta,n,m){
    # n:sample
    cv=qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
    power=numeric(m)
    for(i in 1:m){
      x=rbeta(n,beta,beta)
      p=mean((x-mean(x))^3)/mean((x-mean(x))^2)^1.5
      power[i]=as.integer(abs(p)>=cv)
    }
    return(mean(power))
  }
  for(i in 1:length(beta))
    for(j in 1:length(n)){
      result[j,i]=f(alpha,beta[i],n[j],m)
    }
  df=data.frame(
    sample=n,
    beta1=result[,1],
    beta2=result[,2],
    beta3=result[,3],
    beta4=result[,4],
    beta5=result[,5]
  )
  return(df)
}

#' @title the question 1.2 of homework (t(v))
#' @return a data.frame
#' @examples
#' \dontrun{
#' HW4_1_2()
#' }
#' @export
HW4_1_2=function() {
  alpha=.05
  beta=1:5
  n=c(10,20,30,50,100)
  m=1000
  result=matrix(0,length(n),length(beta))
  f=function(alpha,beta,n,m){
    # n:sample
    cv=qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
    power=numeric(m)
    for(i in 1:m){
      x=rt(n,beta)
      p=mean((x-mean(x))^3)/mean((x-mean(x))^2)^1.5
      power[i]=as.integer(abs(p)>=cv)
    }
    return(mean(power))
  }
  for(i in 1:length(beta))
    for(j in 1:length(n)){
      result[j,i]=f(alpha,beta[i],n[j],m)
    }
  df=data.frame(
    sample=n,
    n1=result[,1],
    n2=result[,2],
    n3=result[,3],
    n4=result[,4],
    n5=result[,5]
  )
  return(df)
}

#' @title the question 2.1 of homework(the Standard normal and the variance is the same)
#' @return a data.frame
#' @examples
#' \dontrun{
#' HW4_2_1()
#' }
#' @export
HW4_2_1=function() {
  alpha=.05
  beta=1:5
  n=c(10,20,30,50,100)
  m=1000
  result=matrix(0,length(n),length(beta))
  f=function(alpha,beta,n,m){
    # n:sample
    cv=qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
    power=numeric(m)
    for(i in 1:m){
      x=rbeta(n,beta,beta)
      p=mean((x-mean(x))^3)/mean((x-mean(x))^2)^1.5
      power[i]=as.integer(abs(p)>=cv)
    }
    return(mean(power))
  }
  for(i in 1:length(beta))
    for(j in 1:length(n)){
      result[j,i]=f(alpha,beta[i],n[j],m)
    }
  df=data.frame(
    sample=n,
    beta1=result[,1],
    beta2=result[,2],
    beta3=result[,3],
    beta4=result[,4],
    beta5=result[,5]
  )
  return(df)
}

#' @title the question 2.2 of homework(the Standard normal and the variance is not the same)
#' @return a data.frame
#' @examples
#' \dontrun{
#' HW4_2_2()
#' }
#' @export
HW4_2_2=function() {
  count5test=function(x,y){
    X=x-mean(x)
    Y=y-mean(y)
    outx=sum(X>max(Y))+sum(X<min(Y))
    outy=sum(Y>max(X))+sum(Y<min(X))
    # return 1 (reject) or 0 (do not reject H0)
    return(as.integer(max(c(outx, outy))>5))
  }
  f=function(mu,singa,n,m){
    mcv=qf(alpha/2,n,n)
    Mcv=qf(1-alpha/2,n,n)
    power5test=numeric(m)
    powerf=numeric(m)
    for(i in 1:m){
      x1=rnorm(n,mu,singa[1])
      x2=rnorm(n,mu,singa[2])
      power5test[i]=count5test(x1,x2)
      powerf[i]=as.integer(var(x2)/var(x1)<=mcv)+as.integer(var(x2)/var(x1)>=Mcv)
    }
    return(c(mean(power5test),mean(powerf)))
  }
  alpha=0.055
  n=c(20,50,100)
  m=1000
  mu=0;singa=c(1,1.5)
  result=matrix(0,2,length(n))
  for(i in 1:length(n)){
    result[,i]=f(mu,singa,n[i],m)
  }
  df=data.frame(
    sample=c("count 5-test","F-test"),
    '20'=result[,1],
    '50'=result[,2],
    '100'=result[,3]
  )
  return(df)
}

#' @title the question 2.3 of homework(T-distribution(n=10) and the variance is the same)
#' @return a data.frame
#' @examples
#' \dontrun{
#' HW4_2_3()
#' }
#' @export
HW4_2_3=function() {
  count5test=function(x,y){
    X=x-mean(x)
    Y=y-mean(y)
    outx=sum(X>max(Y))+sum(X<min(Y))
    outy=sum(Y>max(X))+sum(Y<min(X))
    # return 1 (reject) or 0 (do not reject H0)
    return(as.integer(max(c(outx, outy))>5))
  }
  f=function(tn,n,m){
    mcv=qf(alpha/2,n,n)
    Mcv=qf(1-alpha/2,n,n)
    power5test=numeric(m)
    powerf=numeric(m)
    for(i in 1:m){
      x1=rt(n,tn)
      x2=rt(n,tn)
      power5test[i]=count5test(x1,x2)
      powerf[i]=as.integer(var(x2)/var(x1)<=mcv)+as.integer(var(x2)/var(x1)>=Mcv)
    }
    return(c(mean(power5test),mean(powerf)))
  }
  alpha=0.055
  n=c(20,50,100)
  m=1000
  tn=10
  result=matrix(0,2,length(n))
  for(i in 1:length(n)){
    result[,i]=f(tn,n[i],m)
  }
  df=data.frame(
    sample=c("count 5-test","F-test"),
    '20'=result[,1],
    '50'=result[,2],
    '100'=result[,3]
  )
  return(df)
}

#' @title the question 3.1 of homework(Standard normal)
#' @return a data.frame
#' @examples
#' \dontrun{
#' HW4_3_1()
#' }
#' @export
HW4_3_1=function() {
  f=function(n,m,mcv,Mcv){
    power=numeric(m)
    for(i in 1:m){
      x=rnorm(n)
      y=rnorm(n)
      covxy=cov(x,y)
      #power[i]=as.integer(n*sum((outer(x-mean(x),y-mean(y),"*")/covxy)^3)/n^2/6<mcv)+as.integer(n*sum((outer(x-mean(x),y-mean(y),"*")/covxy)^3)/n^2/6>Mcv)
      power[i]=as.integer(n*sum(((x-mean(x))*(y-mean(y))/covxy)^3)/n^2/6<mcv)+as.integer(n*sum(((x-mean(x))*(y-mean(y))/covxy)^3)/n^2/6>Mcv)
    }
    return(mean(power))
  }
  d=2
  alpha=0.05
  n=c(10,20,30,50,100,500)
  m=1000
  mcv=qchisq(alpha/2,d*(d+1)*(d+2)/6)
  Mcv=qchisq(1-alpha/2,d*(d+1)*(d+2)/6)
  result=numeric(length(n))
  for(i in 1:length(n)){
    result[i]=f(n[i],m,mcv,Mcv)
  }
  df=data.frame(
    sample=n,
    power=result
  )
  return(df)
}

#' @title the question 3.2 of homework(a special distribution)
#' @return a picture
#' @examples
#' \dontrun{
#' HW4_3_2()
#' }
#' @export
HW4_3_2=function() {
  f=function(n,m,mcv,Mcv,e){
    power=numeric(m)
    for(i in 1:m){
      sigmax=sample(c(1,10),replace = TRUE,size = n,prob = c(1-e, e))
      x=rnorm(n,0,sigmax)
      sigmay=sample(c(1,10),replace = TRUE,size = n,prob = c(1-e, e))
      y=rnorm(n,0,sigmay)
      covxy=cov(x,y)
      #power[i]=as.integer(n*sum((outer(x-mean(x),y-mean(y),"*")/covxy)^3)/n^2/6<mcv)+as.integer(n*sum((outer(x-mean(x),y-mean(y),"*")/covxy)^3)/n^2/6>Mcv)
      power[i]=as.integer(n*sum(((x-mean(x))*(y-mean(y))/covxy)^3)/n^2/6<mcv)+as.integer(n*sum(((x-mean(x))*(y-mean(y))/covxy)^3)/n^2/6>Mcv)
    }
    return(mean(power))
  }
  d=2
  alpha=0.05
  n=30
  m=1000
  epsilon=c(seq(0,.15,.01),seq(.15,1,.05))
  mcv=qchisq(alpha/2,d*(d+1)*(d+2)/6)
  Mcv=qchisq(1-alpha/2,d*(d+1)*(d+2)/6)
  result=numeric(length(epsilon))
  for(i in 1:length(epsilon)){
    result[i]=f(n,m,mcv,Mcv,epsilon[i])
  }
  plot(epsilon,result,type="b",xlab=bquote(epsilon),ylim=c(0,1))
  abline(h=.1,lty=3)
  se <- sqrt(result*(1-result)/m)
  lines(epsilon, result+se, lty = 3)
  lines(epsilon, result-se, lty = 3)
}