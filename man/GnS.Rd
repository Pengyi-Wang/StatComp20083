% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp20083R.R
\name{GnS}
\alias{GnS}
\title{Deal with the impact of truncation data on Y}
\usage{
GnS(x)
}
\arguments{
\item{x}{a matrix(x=(Y,U,V)) Y is respond variable and (U,V) is double truncated about y.}
}
\value{
a vector
}
\description{
Deal with the impact of truncation data on Y
}
\examples{
\dontrun{
US=matrix(1,nrow=n,ncol=3)
US[,1]=Y;US[,2]=U;US[,3]=V
TS=GnS(US)
Gn=function(y)
{
  h=0
  for(j in 1:n){if(U[j]<=y&&y<=V[j]){h=h+TS[j]}}
  return(h)
}
}
}
