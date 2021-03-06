#' @title Benchmark9 R and Rcpp functions.
#' @name benchmarks9
#' @description The functions used in homework 9
#' @useDynLib StatComp20083
#' @examples
#' \dontrun{
#' N=2000
#' sigma=c(.05,.5,2,16)
#' 
#' x0=25
#' rw1=rw.Metropolis(sigma[1],x0,N)
#' rw2=rw.Metropolis(sigma[2],x0,N)
#' rw3=rw.Metropolis(sigma[3],x0,N)
#' rw4=rw.Metropolis(sigma[4],x0,N)
#' 
#' crw1=crw(sigma[1],x0,N)
#' crw2=crw(sigma[2],x0,N)
#' crw3=crw(sigma[3],x0,N)
#' crw4=crw(sigma[4],x0,N)
#' 
#' qqplot(rw1$x,crw1$x)
#' qqplot(rw2$x,crw2$x)
#' qqplot(rw3$x,crw3$x)
#' qqplot(rw4$x,crw4$x)
#' 
#' tm1 <- microbenchmark::microbenchmark(
#' rw3 = rw.Metropolis(2,25,2000),
#' crw3 = crw(2,25,2000)
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' }
NULL

#' @title a random walk Metropolis sampler using R
#' @description a random walk Metropolis sampler using R
#' @param sigma the variances are used for the normal distribution
#' @param x0 Initial value of x
#' @param N the number of samples
#' @return a random sample of size N
#' @examples
#' \dontrun{
#' N=2000
#' sigma=c(.05,.5,2,16)
#' x0=25
#' rw1=rw.Metropolis(sigma[1],x0,N)
#' rw2=rw.Metropolis(sigma[2],x0,N)
#' rw3=rw.Metropolis(sigma[3],x0,N)
#' rw4=rw.Metropolis(sigma[4],x0,N)
#' }
#' @export
rw.Metropolis = function(sigma, x0, N) {
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
