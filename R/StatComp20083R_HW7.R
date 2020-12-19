#' @title Benchmark7 R and Rcpp functions.
#' @name benchmarks7
#' @description The functions used in homework 7
#' @importFrom graphics par abline lines abline
#' @examples
#' \dontrun{
#' HW7_1()
#' HW7_2()
#' HW7_3_1()
#' HW7_3_2()
#' }
NULL

#' @title the question 1 of homework
#' @examples
#' \dontrun{
#' HW7_1()
#' }
#' @export
HW7_1=function(){
  set.seed(980904)
  
  rw.Metropolis=function(sigma, x0, N) {
    x=numeric(N)
    x[1]=x0
    u=runif(N)
    k=0
    for(i in 2:N){
      y=rnorm(1,x[i-1],sigma)
      if(u[i]<=((1/2*exp(-abs(y)))/(1/2*exp(-abs(x[i-1])))))
        x[i]=y  
      else{
        x[i]=x[i-1]
        k=k+1
      }
    }
    return(list(x=x,k=k))
  }
  
  N=2000
  sigma=c(.05,.5,2,16)
  
  x0=25
  rw1=rw.Metropolis(sigma[1],x0,N)
  rw2=rw.Metropolis(sigma[2],x0,N)
  rw3=rw.Metropolis(sigma[3],x0,N)
  rw4=rw.Metropolis(sigma[4],x0,N)
  
  #number of candidate points rejected
  no.reject=data.frame(sigma=sigma,no.reject=c(rw1$k,rw2$k,rw3$k,rw4$k))
  knitr::kable(no.reject)
  
  par(mfrow=c(2,2))  #display 4 graphs together
  refline=c(log(1/20),-log(1/20))
  rw=cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
  for (j in 1:4) {
    plot(rw[,j], type="l",
         xlab=bquote(sigma == .(round(sigma[j],3))),
         ylab="X", ylim=range(rw[,j]))
    abline(h=refline)
  }
  par(mfrow=c(1,1)) #reset to default
  
  a=c(.05, seq(.1, .9, .1), .95)
  Q=c(log(2*a[1:6]),-log(2*(1-a[7:11])))
  rw=cbind(rw1$x, rw2$x, rw3$x, rw4$x)
  mc=rw[501:N, ]
  Qrw=apply(mc, 2, function(x) quantile(x, a))
  qq=data.frame(round(cbind(Q, Qrw), 3))
  names(qq)=c('True','sigma=0.05','sigma=0.5','sigma=2','sigma=16')
  knitr::kable(qq)
}

#' @title the question 2 of homework
#' @examples
#' \dontrun{
#' HW7_2()
#' }
#' @export
HW7_2=function(){
  Gelman.Rubin=function(psi) {
    # psi[i,j] is the statistic psi(X[i,1:j])
    # for chain in i-th row of X
    psi=as.matrix(psi)
    n=ncol(psi)
    k=nrow(psi)
    
    psi.means=rowMeans(psi)     #row means
    B=n*var(psi.means)        #between variance est.
    psi.w=apply(psi,1,"var")  #within variances
    W=mean(psi.w)               #within est.
    v.hat=W*(n-1)/n + (B/n)     #upper variance est.
    r.hat=v.hat/W             #G-R statistic
    return(r.hat)
  }
  
  rw.Metropolis=function(sigma, x0, N) {
    x=numeric(N)
    x[1]=x0
    u=runif(N)
    k=0
    for(i in 2:N){
      y=rnorm(1,x[i-1],sigma)
      if(u[i]<=((1/2*exp(-abs(y)))/(1/2*exp(-abs(x[i-1])))))
        x[i]=y  
      else{
        x[i]=x[i-1]
        k=k+1
      }
    }
    return(x)
  }
  
  sigma=2     #parameter of proposal distribution
  k= 4          #number of chains to generate
  n=15000      #length of chains
  b=1000       #burn-in length
  
  #choose overdispersed initial values
  x0=c(-10, -5, 5, 10)
  
  #generate the chains
  set.seed(980904)
  X=matrix(0,nrow=k,ncol=n)
  for (i in 1:k)
    X[i, ]=rw.Metropolis(sigma,x0[i],n)
  
  #compute diagnostic statistics
  psi=t(apply(X, 1, cumsum))
  for(i in 1:nrow(psi))
    psi[i,]=psi[i,]/(1:ncol(psi))
  
  for (i in 1:k)
    if(i==1){
      plot((b+1):n,psi[i,(b+1):n],ylim=c(-0.2,0.2),type="l",
           xlab='Index',ylab=bquote(phi))
    }else{
      lines(psi[i,(b+1):n],col=i)
    }
  par(mfrow=c(1,1)) #restore default
  
  #plot the sequence of R-hat statistics
  pen=0
  rhat <- rep(0, n)
  for (j in (b+1):n){
    rhat[j] <- Gelman.Rubin(psi[,1:j])
    if(pen==0&rhat[j]<=1.1) {pen=j;}
  }
  
  plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
  abline(h=1.1, lty=2)
  print(pen)
}

#' @title the question 3.1 of homework
#' @examples
#' \dontrun{
#' HW7_3_1()
#' }
#' @export
HW7_3_1=function(){
  compute=function(n){
    k=n
    f=function(x){
      (1-pt(sqrt(x^2*(k-1)/(k-x^2)),k-1))-(1-pt(sqrt(x^2*k/(k+1-x^2)),k))
    }
    res1=uniroot(f,c(0.1,sqrt(k)-0.1))
    return(res1$root)
  }
  n=c(4:25,100,500,1000)
  root=numeric(length(n))
  for(i in 1:length(n)){
    root[i]=compute(n[i])
  }
  rw=round(cbind(n,-root,numeric(length(root)),root),3)
  qq=data.frame(rw)
  names(qq)=c('k','The first solution','The second solution','The third solution')
  knitr::kable(qq)
}

#' @title the question 3.2 of homework
#' @examples
#' \dontrun{
#' HW7_3_2()
#' }
#' @export
HW7_3_2=function(){
  compute=function(n){
    k=n
    f=function(x){
      (1-pt(sqrt(x^2*(k-1)/(k-x^2)),k-1))-(1-pt(sqrt(x^2*k/(k+1-x^2)),k))
    }
    res1=uniroot(f,c(0.1,sqrt(k)-0.1))
    return(res1$root)
  }
  n=c(4:25,100,500,1000)
  root=numeric(length(n))
  for(i in 1:length(n)){
    root[i]=compute(n[i])
  }
  rw=round(cbind(n,root),3)
  qq=data.frame(rw)
  names(qq)=c('k','solution')
  knitr::kable(qq)
}