## -----------------------------------------------------------------------------
a=b=2
set.seed(980904)
n=1000
u=runif(n)
x=b/(1-u)^(1/a)
hist(x, prob = TRUE, breaks=seq(b,max(x)+1,1), main = expression(f(x)==a*b^a/x^(a+1)))
y=seq(b,max(x)+1,1)
lines(y,(a*b^a)/(y^(a+1)))

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
r=4
b=2
set.seed(980904)
n=1000
u=runif(n)
y=b/(1-u)^(1/r)-b
hist(y, prob = TRUE, breaks=seq(0,max(y)+1,.1), main = expression(f(y)==r*b^r/(b+y)^(r+1)))
x=seq(0,max(y)+1,.1)
lines(x,(r*b^r)/((b+x)^(r+1)))

## -----------------------------------------------------------------------------
pi=3.141592654
set.seed(980904)
n=1e4
u=runif(n,min=0,max=pi/3)
answer=mean(pi/3*sin(u))
answer.sd=sd(pi/3*sin(u))/sqrt(n)
print(c(answer,1/2,(answer-1/2)/(1/2),round(c(answer-1.96*answer.sd,answer+1.96*answer.sd),6)))

## -----------------------------------------------------------------------------
n=10000;
x1=runif(n,min=0,max=1)
x2=runif(n/2,min=0,max=1)
T_V=0.9839
print(c(var(exp(x1)),var(c((exp(x2)+exp(1-x2))/2)),(var(exp(x1))-var(c((exp(x2)+exp(1-x2))/2)))/var(exp(x1)),T_V ) )

## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## ----echo=FALSE---------------------------------------------------------------
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
knitr::kable(df)

## ----echo=FALSE---------------------------------------------------------------
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
knitr::kable(df)

## ----echo=FALSE---------------------------------------------------------------
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
    x1=rnorm(n,mu,singa)
    x2=rnorm(n,mu,singa)
    power5test[i]=count5test(x1,x2)
    powerf[i]=as.integer(var(x2)/var(x1)<=mcv)+as.integer(var(x2)/var(x1)>=Mcv)
  }
  return(c(mean(power5test),mean(powerf)))
}
alpha=0.055
n=c(20,50,100)
m=1000
mu=0;singa=1
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
knitr::kable(df)

## ----echo=FALSE---------------------------------------------------------------
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
knitr::kable(df)

## ----echo=FALSE---------------------------------------------------------------
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
knitr::kable(df)

## ----echo=FALSE---------------------------------------------------------------
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
knitr::kable(df)

## ----echo=FALSE---------------------------------------------------------------
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

## -----------------------------------------------------------------------------
x=matrix(0,2,2)
x[1,1]=6510;x[1,2]=3490;
x[2,1]=6760;x[2,2]=3240;
X=sum(x)*(x[1,1]*x[2,2]-x[1,2]*x[2,1])^2/(sum(x[1,])*sum(x[2,])*sum(x[,1])*sum(x[,2]))
print(c(X,qchisq(0.95,1),as.integer(X>qchisq(0.95,1))))

## ----echo=TRUE----------------------------------------------------------------
library(bootstrap)
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

## ----echo=FALSE---------------------------------------------------------------
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


## ----echo=FALSE---------------------------------------------------------------
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


## ----echo=FALSE---------------------------------------------------------------
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


## ----echo=FALSE---------------------------------------------------------------
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


## ----echo=TRUE----------------------------------------------------------------
library(bootstrap)
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

## -----------------------------------------------------------------------------
library(DAAG);
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




## ----echo=TRUE----------------------------------------------------------------
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

## ----echo=TRUE----------------------------------------------------------------
library(boot)
library(RANN)
library(energy)
library(Ball)
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

## ----echo=TRUE----------------------------------------------------------------
library(boot)
library(RANN)
library(energy)
library(Ball)
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

## ----echo=TRUE----------------------------------------------------------------
library(boot)
library(RANN)
library(energy)
library(Ball)
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

## ----echo=TRUE----------------------------------------------------------------
library(boot)
library(RANN)
library(energy)
library(Ball)
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

## ----echo=TRUE----------------------------------------------------------------
library(boot)
library(RANN)
library(energy)
library(Ball)
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

## ----echo=TRUE----------------------------------------------------------------
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


refline=c(log(1/20),-log(1/20))
rw=cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
            xlab=bquote(sigma == .(round(sigma[j],3))),
            ylab="X", ylim=range(rw[,j]))
        abline(h=refline)
    }

a=c(.05, seq(.1, .9, .1), .95)
Q=c(log(2*a[1:6]),-log(2*(1-a[7:11])))
rw=cbind(rw1$x, rw2$x, rw3$x, rw4$x)
mc=rw[501:N, ]
Qrw=apply(mc, 2, function(x) quantile(x, a))
qq=data.frame(round(cbind(Q, Qrw), 3))
names(qq)=c('True','sigma=0.05','sigma=0.5','sigma=2','sigma=16')
knitr::kable(qq)

## ----echo=TRUE----------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
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

## ----echo=TRUE,warning=FALSE--------------------------------------------------
library(StatComp20083)
HW8_1()


## ----echo=TRUE----------------------------------------------------------------
library(DAAG)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1/disp),
  mpg ~ disp+wt,
  mpg ~ I(1/disp)+wt
)
mpg=mtcars$mpg
disp=mtcars$disp
wt=mtcars$wt
lanswer=foranswer=c(NA)
for(i in 1:length(formulas)){
  foranswer=c(foranswer,as.numeric(lm(as.character(formulas[i]))$coefficients),NA)
}
print(foranswer)

## -----------------------------------------------------------------------------
f=function(x){
  return(as.numeric(lm(x)$coefficients))
}
print( lapply( formulas,f ) )

## -----------------------------------------------------------------------------
set.seed(980904)
trials <- replicate(
  100,
  t.test(rpois(10,10), rpois(7,10)),
  simplify = FALSE
)
sapply(trials,function(x) x$p.value)

## -----------------------------------------------------------------------------
f=function(x) x$p.value
sapply(trials,f)

## -----------------------------------------------------------------------------
datalist=list(mtcars,faithful)
lapply(datalist,function(x) vapply(x,mean,numeric(1)))

## -----------------------------------------------------------------------------
mylapply=function(X,FUN,FUN.VALUE,simplify=FALSE){
  out=Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify==TRUE) return(simplify2array(out))
  unlist(out)
}
mylapply(datalist, mean, numeric(1))

## ----eval=FALSE---------------------------------------------------------------
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  #define PI 3.141592654
#  // [[Rcpp::export]]
#  List crw(double sigma, double x0, int N){
#    std::vector<double> x;
#    std::vector<int> knum;
#    double u[N],y;
#    int k=0;
#    x.push_back(x0);
#    int x1=x0;
#  
#    double U,V;
#    int phase=0;
#    double z;
#  
#    for(int i=0;i<N;i++){u[i]=rand()/(RAND_MAX+1.0);}
#    for(int i=1;i<N;i++){
#      y=rnorm(1,x1,sigma)[0];
#      double u_compare=exp(-abs(y)+abs(x1));
#      if(u[i]<=u_compare )
#        {x.push_back(y);x1=y;}
#      else{
#        x.push_back(x1);
#        k=k+1;
#      }
#    }
#    knum.push_back(k);
#    return List::create(
#      _["x"] = x,
#      _["k"] = knum
#    );
#  }

## -----------------------------------------------------------------------------
library(Rcpp)
set.seed(980904)
dir_cpp="../vignettes/"
sourceCpp(paste0(dir_cpp,"crw.cpp"))

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

crw1=crw(sigma[1],x0,N)
crw2=crw(sigma[2],x0,N)
crw3=crw(sigma[3],x0,N)
crw4=crw(sigma[4],x0,N)

#number of candidate points rejected
no.reject=data.frame(sigma=sigma,no.reject=c(crw1$k,crw2$k,crw3$k,crw4$k))
knitr::kable(no.reject)

refline=c(log(1/20),-log(1/20))
rw=cbind(crw1$x, crw2$x, crw3$x,  crw4$x)
for (j in 1:4) {
    plot(rw[,j], type="l",
        xlab=bquote(sigma == .(round(sigma[j],3))),
        ylab="X", ylim=range(rw[,j]))
    abline(h=refline)
}

a=c(.05, seq(.1, .9, .1), .95)
Q=c(log(2*a[1:6]),-log(2*(1-a[7:11])))
rw=cbind(crw1$x, crw2$x, crw3$x, crw4$x)
mc=rw[501:N, ]
Qrw=apply(mc, 2, function(x) quantile(x, a))
qq=data.frame(round(cbind(Q, Qrw), 3))
names(qq)=c('True','sigma=0.05','sigma=0.5','sigma=2','sigma=16')
knitr::kable(qq)

## -----------------------------------------------------------------------------
qqplot(rw1$x,crw1$x)
qqplot(rw2$x,crw2$x)
qqplot(rw3$x,crw3$x)
qqplot(rw4$x,crw4$x)

## -----------------------------------------------------------------------------
library(microbenchmark)
N=2000
x0=25
ts <- microbenchmark(rw3=rw.Metropolis(2,x0,N),crw3=crw(2,x0,N))
summary(ts)[,c(1,3,5,6)]

