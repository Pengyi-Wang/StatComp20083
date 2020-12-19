#' @title Benchmark8 R and Rcpp functions.
#' @name benchmarks8
#' @description The functions used in homework 7
#' @examples
#' \dontrun{
#' HW8_1()
#' HW8_2_1(mtcars)
#' HW8_2_2(mtcars)
#' HW8_3_1()
#' HW8_3_2()
#' HW8_4_1(mtcars, faithful)
#' HW8_4_2(mtcars, faithful)
#' }
NULL

#' @title the question 1 of homework
#' @examples
#' \dontrun{
#' HW8_1()
#' }
#' @export
HW8_1=function(){
  n=c(444,132,361,63)
  p0=0.1;q0=0.1;
  p1=1;q1=1;
  p2=q2=c(0.1)
  while( (abs(p0-p1)+abs(q0-q1))>0.0000001 ){
    mlogL=function(p=0.1,q=0.1){
      return( -(2*n[3]*log(1-p-q)+n[1]*(log(p*(1-p-q))+p0/(2-p0-2*q0)*log(p/(1-p-q)))+n[2]*(log(q*(1-p-q))+q0/(2-q0-2*p0)*log(q/(1-p-q)))+n[4]*log(p*q) ))
    }
    
    fit=mle(mlogL)
    q1=q0;p1=p0;
    p0=as.numeric(fit@coef)[1]
    q0=as.numeric(fit@coef)[2]
    p2=c(p2,p0)
    q2=c(q2,q0)
  }
  print(c(p0,q0))
  plot(1:length(p2),log(p2))
  plot(1:length(q2),log(q2))
}

#' @title the question 2.1 of homework(for())
#' @param mtcars a data frame (use datasets mtcars)
#' @examples
#' \dontrun{
#' HW8_2_1()
#' }
#' @export
HW8_2_1=function(mtcars){
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
}

#' @title the question 2.2 of homework(lapply())
#' @param mtcars a data frame (use datasets mtcars)
#' @examples
#' \dontrun{
#' HW8_2_2()
#' }
#' @export
HW8_2_2=function(mtcars){
  formulas <- list(
    mpg ~ disp,
    mpg ~ I(1/disp),
    mpg ~ disp+wt,
    mpg ~ I(1/disp)+wt
  )
  mpg=mtcars$mpg
  disp=mtcars$disp
  wt=mtcars$wt
  f=function(x){
    return(as.numeric(lm(x)$coefficients))
  }
  print( lapply( formulas,f ) )
}

#' @title the question 3.1 of homework
#' @examples
#' \dontrun{
#' HW8_3_1()
#' }
#' @export
HW8_3_1=function(){
  set.seed(980904)
  trials <- replicate(
    100,
    t.test(rpois(10,10), rpois(7,10)),
    simplify = FALSE
  )
  sapply(trials,function(x) x$p.value)
}

#' @title the question 3.2 of homework
#' @examples
#' \dontrun{
#' HW8_3_2()
#' }
#' @export
HW8_3_2=function(){
  trials <- replicate(
    100,
    t.test(rpois(10,10), rpois(7,10)),
    simplify = FALSE
  )
  f=function(x) x$p.value
  sapply(trials,f)
}

#' @title the question 4.1 of homework
#' @param mtcars a data frame (use datasets mtcars)
#' @param faithful a data frame (use datasets faithful)
#' @examples
#' \dontrun{
#' HW8_4_1()
#' }
#' @export
HW8_4_1=function(mtcars, faithful){
  datalist=list(mtcars, faithful)
  lapply(datalist, function(x) vapply(x, mean, numeric(1)))
}

#' @title the question 4.2 of homework
#' @param mtcars a data frame (use datasets mtcars)
#' @param faithful a data frame (use datasets faithful)
#' @examples
#' \dontrun{
#' HW8_4_2()
#' }
#' @export
HW8_4_2=function(mtcars, faithful){
  datalist=list(mtcars, faithful)
  mylapply=function(X, FUN, FUN.VALUE, simplify = FALSE){
    out=Map(function(x) vapply(x, FUN, FUN.VALUE), X)
    if(simplify == TRUE) return(simplify2array(out))
    unlist(out)
  }
  mylapply(datalist, mean, numeric(1))
}
