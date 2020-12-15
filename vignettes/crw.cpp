#include <Rcpp.h>
using namespace Rcpp;
#define PI 3.141592654
// [[Rcpp::export]]
List crw(double sigma, double x0, int N){
  std::vector<double> x;
  std::vector<int> knum;
  double u[N],y;
  int k=0;
  x.push_back(x0);
  int x1=x0;
  for(int i=0;i<N;i++){u[i]=rand()/(RAND_MAX+1.0);}
  for(int i=1;i<N;i++){
    y=rnorm(1,x1,sigma)[0];
    double u_compare=exp(-abs(y)+abs(x1));
    if(u[i]<=u_compare )
      {x.push_back(y);x1=y;}
    else{
      x.push_back(x1);
      k=k+1;
    }
  }
  knum.push_back(k);
  return List::create(
    _["x"] = x,
    _["k"] = knum
  );
}