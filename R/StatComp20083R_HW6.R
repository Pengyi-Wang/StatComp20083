#' @title Benchmark6 R and Rcpp functions.
#' @name benchmarks6
#' @description The functions used in homework 6
#' @import boot RANN energy Ball
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm runif
#' @useDynLib StatComp20083
#' @examples
#' \dontrun{
#' HW6_1()
#' HW6_2_1()
#' HW6_2_2()
#' HW6_2_3()
#' HW6_2_4()
#' HW6_2_5()
#' }
NULL

#' @title the question 1 of homework
#' @return a number about test
#' @examples
#' \dontrun{
#' HW6_1()
#' }
#' @export
HW6_1=function(){
  maxout <- function(x, y) {
    X=x-mean(x)
    Y=y-mean(y)
    outx=sum(X>max(Y))+sum(X<min(Y))
    outy=sum(Y>max(X))+sum(Y<min(X))
    return(max(c(outx, outy)))
  }
  count5test=function(x,y,num){
    X=x-mean(x)
    Y=y-mean(y)
    outx=sum(X>max(Y))+sum(x<min(Y))
    outy=sum(Y>max(X))+sum(Y<min(X))
    # return 1 (reject) or 0 (do not reject H0)
    return(as.integer(max(c(outx, outy))>num))
  }
  n1=30
  n2=50
  m=1000
  set.seed(980904)
  stat <- replicate(m,expr={
    x=rnorm(n1)
    y=rnorm(n2)
    maxout(x, y)
  })
  num=as.integer(quantile(stat, c(.95)))-1
  B=1000
  test=numeric(B)
  for(i in 1:B){
    x=rnorm(n1);y=rnorm(n2)
    test[i]=count5test(x,y,num)
  }
  print(mean(test))
}

#' @title the question 2.1 of homework(a normal distribution with a mean of 0 and a variance of 1 and 1.5)
#' @return some number
#' @examples
#' \dontrun{
#' HW6_2_1()
#' }
#' @export
HW6_2_1=function(){
  Tn=function(z, ix, sizes,k){
    n1=sizes[1];n2=sizes[2];n=n1+n2
    if(is.vector(z))z=data.frame(z,0);
    z=z[ix, ];
    NN=nn2(data=z, k=k+1)
    block1=NN$nn.idx[1:n1,-1]
    block2=NN$nn.idx[(n1+1):n,-1]
    i1=sum(block1<n1+.5);i2=sum(block2>n1+.5)
    return((i1+i2)/(k*n))
  }
  
  n=20
  set.seed(980904)
  x=matrix(rnorm(n*2),nrow=n,ncol=2)
  y=matrix(rnorm(n*2,0,1.5),nrow=n,ncol=2)
  z=rbind(x,y)
  p.value=numeric(3)
  N=c(nrow(x),nrow(y))
  
  boot.obj=boot(data=z,statistic=Tn,R=999,sim="permutation",sizes=N,k=3)
  ts=c(boot.obj$t0,boot.obj$t)
  p.value[1]=mean(ts>=ts[1])
  
  boot.obs=eqdist.etest(z,sizes=N,R=999)
  p.value[2]=boot.obs$p.value
  
  p.value[3]=bd.test(x=x,y=y,R=999)$p.value
  
  print(p.value)
}

#' @title the question 2.2 of homework(a normal distribution with a mean of 0 and 1 and a variance of 1)
#' @return some number
#' @examples
#' \dontrun{
#' HW6_2_2()
#' }
#' @export
HW6_2_2=function(){
  Tn=function(z, ix, sizes,k){
    n1=sizes[1];n2=sizes[2];n=n1+n2
    if(is.vector(z))z=data.frame(z,0);
    z=z[ix, ];
    NN=nn2(data=z, k=k+1)
    block1=NN$nn.idx[1:n1,-1]
    block2=NN$nn.idx[(n1+1):n,-1]
    i1=sum(block1<n1+.5);i2=sum(block2>n1+.5)
    return((i1+i2)/(k*n))
  }
  
  n=20
  set.seed(980904)
  x=matrix(rnorm(n*2),nrow=n,ncol=2)
  y=matrix(rnorm(n*2,1,1),nrow=n,ncol=2)
  z=rbind(x,y)
  p.value=numeric(3)
  N=c(nrow(x),nrow(y))
  
  boot.obj=boot(data=z,statistic=Tn,R=999,sim="permutation",sizes=N,k=3)
  ts=c(boot.obj$t0,boot.obj$t)
  p.value[1]=mean(ts>=ts[1])
  
  boot.obs=eqdist.etest(z,sizes=N,R=999)
  p.value[2]=boot.obs$p.value
  
  p.value[3]=bd.test(x=x,y=y,R=999)$p.value
  
  print(p.value)
}

#' @title the question 2.3 of homework(a t distribution with df=1)
#' @return some number
#' @examples
#' \dontrun{
#' HW6_2_3()
#' }
#' @export
HW6_2_3=function(){
  Tn=function(z, ix, sizes,k){
    n1=sizes[1];n2=sizes[2];n=n1+n2
    if(is.vector(z))z=data.frame(z,0);
    z=z[ix, ];
    NN=nn2(data=z, k=k+1)
    block1=NN$nn.idx[1:n1,-1]
    block2=NN$nn.idx[(n1+1):n,-1]
    i1=sum(block1<n1+.5);i2=sum(block2>n1+.5)
    return((i1+i2)/(k*n))
  }
  
  n=20
  set.seed(980904)
  x=matrix(rt(n*1,1),nrow=n,ncol=1)
  y=matrix(rt(n*1,1),nrow=n,ncol=1)
  z=rbind(x,y)
  p.value=numeric(3)
  N=c(nrow(x),nrow(y))
  
  boot.obj=boot(data=z,statistic=Tn,R=999,sim="permutation",sizes=N,k=3)
  ts=c(boot.obj$t0,boot.obj$t)
  p.value[1]=mean(ts>=ts[1])
  
  boot.obs=eqdist.etest(z,sizes=N,R=999)
  p.value[2]=boot.obs$p.value
  
  p.value[3]=bd.test(x=x,y=y,R=999)$p.value
  
  print(p.value)
}

#' @title the question 2.4 of homework(A particular distribution)
#' @return some number
#' @examples
#' \dontrun{
#' HW6_2_4()
#' }
#' @export
HW6_2_4=function(){
  Tn=function(z, ix, sizes,k){
    n1=sizes[1];n2=sizes[2];n=n1+n2
    if(is.vector(z))z=data.frame(z,0);
    z=z[ix, ];
    NN=nn2(data=z, k=k+1)
    block1=NN$nn.idx[1:n1,-1]
    block2=NN$nn.idx[(n1+1):n,-1]
    i1=sum(block1<n1+.5);i2=sum(block2>n1+.5)
    return((i1+i2)/(k*n))
  }
  
  n=20
  set.seed(980904)
  complex=matrix(sample(c(0,1),size=2*2*n,replace=TRUE),nrow=2,ncol=2*n)
  x=matrix(rnorm(n*2,0,complex[1,]+(1-complex[1,])*100),nrow=n,ncol=2)
  y=matrix(rnorm(n*2,0,complex[2,]+(1-complex[2,])*100),nrow=n,ncol=2)
  z=rbind(x,y)
  p.value=numeric(3)
  N=c(nrow(x),nrow(y))
  
  boot.obj=boot(data=z,statistic=Tn,R=999,sim="permutation",sizes=N,k=3)
  ts=c(boot.obj$t0,boot.obj$t)
  p.value[1]=mean(ts>=ts[1])
  
  boot.obs=eqdist.etest(z,sizes=N,R=999)
  p.value[2]=boot.obs$p.value
  
  p.value[3]=bd.test(x=x,y=y,R=999)$p.value
  
  print(p.value)
}

#' @title the question 2.5 of homework(a normal distribution with a mean of 0 and a variance of 1, and the sample size of x is 20 and the sample size of y is 40)
#' @return some number
#' @examples
#' \dontrun{
#' HW6_2_5()
#' }
#' @export
HW6_2_5=function(){
  Tn=function(z, ix, sizes,k){
    n1=sizes[1];n2=sizes[2];n=n1+n2
    if(is.vector(z))z=data.frame(z,0);
    z=z[ix, ];
    NN=nn2(data=z, k=k+1)
    block1=NN$nn.idx[1:n1,-1]
    block2=NN$nn.idx[(n1+1):n,-1]
    i1=sum(block1<n1+.5);i2=sum(block2>n1+.5)
    return((i1+i2)/(k*n))
  }
  
  n1=20
  n2=40
  set.seed(980904)
  x=matrix(rnorm(n1*2),nrow=n1,ncol=2)
  y=matrix(rnorm(n2*2),nrow=n2,ncol=2)
  z=rbind(x,y)
  p.value=numeric(3)
  N=c(nrow(x),nrow(y))
  
  boot.obj=boot(data=z,statistic=Tn,R=999,sim="permutation",sizes=N,k=3)
  ts=c(boot.obj$t0,boot.obj$t)
  p.value[1]=mean(ts>=ts[1])
  
  boot.obs=eqdist.etest(z,sizes=N,R=999)
  p.value[2]=boot.obs$p.value
  
  p.value[3]=bd.test(x=x,y=y,R=999)$p.value
  
  print(p.value)
}
