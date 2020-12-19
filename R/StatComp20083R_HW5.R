#' @title Benchmark5 R and Rcpp functions.
#' @name benchmarks5
#' @description The functions used in homework 5
#' @examples
#' \dontrun{
#' HW5_1(law)
#' HW5_2_1()
#' HW5_2_2()
#' HW5_2_3()
#' HW5_2_4()
#' HW5_3(scor)
#' HW5_4(ironslag)
#' }
NULL

#' @title the question 1 of homework
#' @param law a data frame (use bootstrap law)
#' @return two number about bias.jack and se.jack
#' @examples
#' \dontrun{
#' HW5_1(law)
#' }
#' @export
HW5_1=function(law){
  n=length(law$LSAT)
  x=matrix(0,nrow=2,ncol=n)
  x[1,]=law$LSAT;x[2,]=law$GPA;
  theta.hat=cor(law$LSAT,law$GPA)
  theta.jack=numeric(n)
  for(i in 1:n){
    theta.jack[i]=cor(x[1,(1:n)[-i]],x[2,(1:n)[-i]])
  }
  bias.jack=(n-1)*(mean(theta.jack)-theta.hat)
  se.jack=(n-1)/sqrt(n)*sd(theta.jack)
  print(c(bias.jack,se.jack))
}

#' @title the question 2.1 of homework (The standard bootstrap CI based on standard normal)
#' @return some number
#' @examples
#' \dontrun{
#' HW5_2_1()
#' }
#' @export
HW5_2_1=function(){
  alpha=0.05
  x=c(3,5,7,18,43,85,91,98,100,130,230,487)
  B=1e4;set.seed(980904);theta.star=numeric(B)
  theta.hat=mean(x)
  for(b in 1:B){
    x.star=sample(x,replace=TRUE)
    theta.star[b]=mean(x.star)
  }
  se.hat=sd(theta.star)
  print(c(theta.hat-qnorm(1-alpha/2)*se.hat,theta.hat+qnorm(1-alpha/2)*se.hat))
}

#' @title the question 2.2 of homework (The standard bootstrap CI based on basic)
#' @return some number
#' @examples
#' \dontrun{
#' HW5_2_2()
#' }
#' @export
HW5_2_2=function(){
  alpha=0.05
  x=c(3,5,7,18,43,85,91,98,100,130,230,487)
  B=1e4;set.seed(980904);theta.star=numeric(B)
  theta.hat=mean(x)
  for(b in 1:B){
    x.star=sample(x,replace=TRUE)
    theta.star[b]=mean(x.star)
  }
  x.p=as.numeric(quantile(theta.star,probs=c(1-alpha/2,alpha/2)))
  print(c(2*theta.hat-x.p[1],2*theta.hat+x.p[2]))
}

#' @title the question 2.3 of homework (The standard bootstrap CI based on percentile)
#' @return some number
#' @examples
#' \dontrun{
#' HW5_2_3()
#' }
#' @export
HW5_2_3=function(){
  alpha=0.05
  x=c(3,5,7,18,43,85,91,98,100,130,230,487)
  B=1e4;set.seed(980904);theta.star=numeric(B)
  theta.hat=mean(x)
  for(b in 1:B){
    x.star=sample(x,replace=TRUE)
    theta.star[b]=mean(x.star)
  }
  x.p=as.numeric(quantile(theta.star,probs=c(alpha/2,1-alpha/2)))
  print(c(x.p[1],x.p[2]))
}

#' @title the question 2.4 of homework (The standard bootstrap CI based on BCa)
#' @return some number
#' @examples
#' \dontrun{
#' HW5_2_4()
#' }
#' @export
HW5_2_4=function(){
  alpha=0.05
  x=c(3,5,7,18,43,85,91,98,100,130,230,487)
  n=length(x)
  B=1e4
  set.seed(980904)
  theta.star=numeric(B)
  theta.hat=mean(x)
  theta.bar=numeric(n)
  for(b in 1:B){
    x.star=sample(x,replace=TRUE)
    theta.star[b]=mean(x.star)
  }
  for(i in 1:n){
    theta.bar[i]=mean(x[-i])
  }
  z_0=qnorm(1/B*sum(as.integer(theta.star<theta.hat)))
  L=mean(theta.bar)-theta.bar
  alpha.hat=sum(L^3)/(6*sum(L^2)^1.5)
  alpha_1=pnorm(z_0+(z_0+qnorm(alpha/2))/(1-alpha.hat*(z_0+qnorm(alpha/2))))
  alpha_2=pnorm(z_0+(z_0+qnorm(1-alpha/2))/(1-alpha.hat*(z_0+qnorm(1-alpha/2))))
  x.p=as.integer(quantile(theta.star,probs=c(alpha_1,alpha_2)))
  print(c(x.p[1],x.p[2]))
}

#' @title the question 3 of homework
#' @param scor a data frame (use bootstrap scor)
#' @return some number
#' @examples
#' \dontrun{
#' HW5_3(scor)
#' }
#' @export
HW5_3=function(scor){
  score=scor
  n=length(score[,1])
  covmatrix1=covmatrix=matrix(0,nrow=5,ncol=5)
  for(i in 1:5)
    for(j in i:5) covmatrix[i,j]=covmatrix[j,i]=cov(score[,i],score[,j])
  lambdavalues=eigen(covmatrix)$values
  theta.hat=max(lambdavalues)/sum(lambdavalues)
  theta.star=numeric(n)
  for(k in 1:n){
    for(i in 1:5)
      for(j in i:5) covmatrix1[i,j]=covmatrix[j,i]=cov(score[-k,i],score[-k,j])
      theta.star[k]=max(eigen(covmatrix)$values)/sum(eigen(covmatrix)$values)
  }
  bias.jack=(n-1)*(mean(theta.star)-theta.hat)
  se.jack=(n-1)/sqrt(n)*sd(theta.star)
  print(c(bias.jack,se.jack))
}

#' @title the question 4 of homework
#' @param ironslag a data frame (use DAAG ironslag)
#' @return some number
#' @examples
#' \dontrun{
#' HW5_4()
#' }
#' @export
HW5_4=function(ironslag){
  magnetic=ironslag$magnetic
  chemical=ironslag$chemical
  n=length(magnetic)
  e1=e2=e3=e4=numeric(n*(n-1)/2)
  k=1
  for(i in 1:n)
    for(j in (i+1):n){
      if(i==n) break
      y=magnetic[c(-i,-j)]
      x=chemical[c(-i,-j)]
      
      L1=lm(y ~ x)
      yhat1=L1$coef[1]+L1$coef[2]*chemical[c(i,j)]
      e1[k]=sum((magnetic[c(i,j)]-yhat1)^2)
      
      L2=lm(y ~ x+I(x^2))
      yhat2=L2$coef[1]+L2$coef[2]*chemical[c(i,j)]+L2$coef[3]*chemical[c(i,j)]^2
      e2[k]=sum((magnetic[c(i,j)]-yhat2)^2)
      
      L3=lm(log(y) ~ x)
      yhat3=L3$coef[1]+L3$coef[2]*chemical[c(i,j)]
      e3[k]=sum((magnetic[c(i,j)]-exp(yhat3))^2)
      
      L4=lm(log(y) ~ log(x))
      yhat4=L4$coef[1]+L4$coef[2]*log(chemical[c(i,j)])
      e4[k]=sum((magnetic[c(i,j)]-exp(yhat4))^2)
      k=k+1
    }
  print(c(mean(e1),mean(e2),mean(e3),mean(e4)))
}
