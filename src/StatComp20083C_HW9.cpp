#include <Rcpp.h>
using namespace Rcpp;

//' @title a random walk Metropolis sampler using R
//' @description a random walk Metropolis sampler using R
//' @param sigma the variances are used for the normal distribution
//' @param x0 Initial value of x
//' @param N the number of samples
//' @return a random sample of size
//' @examples
//' \dontrun{
//' N=2000
//' sigma=c(.05,.5,2,16)
//' x0=25
//' crw1=crw(sigma[1],x0,N)
//' crw2=crw(sigma[2],x0,N)
//' crw3=crw(sigma[3],x0,N)
//' crw4=crw(sigma[4],x0,N)
//' }
//' @export 
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
