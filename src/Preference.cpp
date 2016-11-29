#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

//Mean Square Error
double double GaussianPreference(double delta, Eigen::VectorXd parms){
  double sigma = parms(0);
  double res = 1-std::exp(-1.0*std::pow(delta,2.0));
  return(res);
}
