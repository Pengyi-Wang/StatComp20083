#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description The functions used in my function (It's a function of my undergraduate paper)
#' @importFrom  stats integrate
#' @useDynLib StatComp20083
#' @examples
#' \dontrun{
#' n=100 
#' m=5    
#' max=1
#' min=0
#' i=1   
#' fangcha=0.01 
#' X=matrix(-1,nrow=n,ncol=m)
#' Y=rep(-1,n);U=Y;V=Y
#' while(Y[n]==-1)
#' {
#'   X[i,]=runif(m,min=min,max=max)
#'   gt=function(x){X[i,1]+X[i,2]*x+X[i,3]*x^2+X[i,4]*x^3+X[i,5]*x^4}
#'   y=as.numeric(integrate(gt,0,1)[1])+rnorm(1,mean=0,sd=fangcha)
#'   u=runif(1,min=-0.5,max=1.15);v=runif(1,min=1.15,max=2.8)
#'   if(u<=y&&y<=v){Y[i]=y;U[i]=u;V[i]=v;i=i+1}
#' }
#' rm(y)
#' Gn=creat(Y,U,V)
#' h1=seq(0.1,2,0.1)
#' h=CV(X,Y,U,V,h1)
#' h1=seq(max(0.01,h-0.1),h+0.1,0.01)
#' h=CV(X,Y,U,V,h1)
#' h1=seq(max(0.001,h-0.01),h+0.01,0.001)
#' h=CV(X,Y,U,V,h1)
#' h1=seq(max(0.0001,h-0.001),h+0.001,0.0001)
#' h=CV(X,Y,U,V,h1)
#' }
NULL

#' @title A kernel function
#' @param x The data values
#' @return a number
#' @examples
#' \dontrun{
#' kernel(1)
#' }
#' @export
Kernel=function(x)
{
  h=1/sqrt(2*3.14159)*exp(-x^2/2)
  return(h)
}

#' @title A measure function of a quartic polynomial
#' @description Converts a functional variable to a scalar
#' @param x1 The coefficient of a quartic polynomial(f(x)=x1[1]+x1[2]x+x1[3]x^2+x1[4]x^3+x1[5]x^4)
#' @param x2 The coefficient of a quartic polynomial(like x1)
#' @return a scalar
#' @examples
#' \dontrun{
#' x1=c(0,1,2,3,4)
#' x2=c(1,2,3,4,5)
#' d(x1,x2)
#' }
#' @export
d=function(x1,x2)
{
  x3=x1-x2
  fff=function(x){(x3[1]+x3[2]*x+x3[3]*x^2+x3[4]*x^3+x3[5]*x^4)^2}  
  return(sqrt(as.numeric(integrate(fff,0,1)[1])))
}

#' @title Deal with the impact of truncation data on Y
#' @param x a matrix(x=(Y,U,V)) Y is respond variable and (U,V) is double truncated about y.
#' @return a vector
#' @examples
#' \dontrun{
#' US=matrix(1,nrow=n,ncol=3)
#' US[,1]=Y;US[,2]=U;US[,3]=V
#' TS=GnS(US)
#' Gn=function(y)
#' {
#'   h=0
#'   for(j in 1:n){if(U[j]<=y&&y<=V[j]){h=h+TS[j]}}
#'   return(h)
#' }
#' }
#' @export
GnS=function(x)    # x=(Y,U,V)
{
  n=nrow(x)
  H=1:n
  h1=1:n
  h2=1:n
  F=1:n
  f=1:n
  I=function(t)
  {
    if(t[2]<=t[1]&t[1]<=t[3]){return(1)}
    else {return(0)}
  }
  for(i in 1:n){h1[i]=1/n}
  for(i in 1:n)
  {
    H[i]=0
    for(j in 1:n){H[i]=H[i]+h1[j]*I(c(x[j,1],x[i,2],x[i,3]))}
  }
  while(abs(max(h1-h2))>=10^(-6))
  {
    h2=h1
    for(j in 1:n)
    {
      f[j]=0
      for(i in 1:n){f[j]=f[j]+1/H[i]}
      f[j]=(1/f[j])*(1/H[j])
    }
    for(i in 1:n)
    {
      F[i]=0
      for(j in 1:n){F[i]=F[i]+f[j]*I(c(x[i,1],x[j,2],x[j,3]))}
    }
    for(j in 1:n)
    {
      h1[j]=0
      for(i in 1:n){h1[j]=h1[j]+1/F[i]}
      h1[j]=(1/h1[j])*(1/F[j])
    }
    for(i in 1:n)
    {
      H[i]=0
      for(j in 1:n){H[i]=H[i]+h1[j]*I(c(x[j,1],x[i,2],x[i,3]))}
    }
  }
  return(f)
}

#' @title Deal with the impact of truncation data on Y
#' @param Y a respond variable.
#' @param U double truncated about y.
#' @param V double truncated about y.
#' @return a function
#' @examples
#' \dontrun{
#' Y=c(1,1,1)
#' U=c(0.1,0.1,0.1)
#' V=c(2,2,2)
#' Gn=creat(Y,U,V,n)
#' }
#' @export
creat=function(Y,U,V){
  n=length(Y)
  US=matrix(1,nrow=n,ncol=3)
  US[,1]=Y;US[,2]=U;US[,3]=V
  TS=GnS(US)
  Gn=function(y)
  {
    h=0
    for(j in 1:n){if(U[j]<=y&&y<=V[j]){h=h+TS[j]}}
    return(h)
  }
  return(Gn)
}


#' @title Select the optimal Bandwidth from a set of data(h1).
#' @description Optimize the Bandwidth for H1 with different accuracy.
#' @param X The coefficient matrix of a quartic polynomial (n rows and 5 columns).
#' @param Y a respond variable.
#' @param U double truncated about y.
#' @param V double truncated about y.
#' @param h1 a set of data.
#' @return a bandwidth from h1.
#' @examples
#' \dontrun{
#' h1=seq(0.1,2,0.1)
#' h=CV(X,Y,U,V,h1)
#' h1=seq(max(0.01,h-0.1),h+0.1,0.01)
#' h=CV(X,Y,U,V,h1)
#' h1=seq(max(0.001,h-0.01),h+0.01,0.001)
#' h=CV(X,Y,U,V,h1)
#' h1=seq(max(0.0001,h-0.001),h+0.001,0.0001)
#' h=CV(X,Y,U,V,h1)
#' }
#' @export
CV=function(X,Y,U,V,h1)
{
  Gn=creat(Y,U,V)
  n=nrow(X)
  mNW=function(x,h,y,x1)
  {
    n1=nrow(x)
    f1=0
    f2=0
    for(i in 1:n1){f1=f1+Kernel(d(x1,x[i,])/h)*y[i]/Gn(y[i])}
    for(i in 1:n1){f2=f2+Kernel(d(x1,x[i,])/h)/Gn(y[i])}
    return(f1/f2)
  }
  CVH=function(x,h,y)
  {
    f3=0
    n2=nrow(x)
    for(i in 1:n2){f3=f3+(y[i]-mNW(x[-i,],h,y[-i],x[i,]))^2}
    return(f3)
  }
  y=rep(1,length(h1))
  for(j in 1:length(h1)){y[j]=CVH(X,h1[j],Y)}	
  place=which(y==min(y))
  h=h1[place]
  return(h)
}

