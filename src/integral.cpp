// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;


class GaussianPrefKernel: public Func {
private:
  // double y;
  Eigen::VectorXd vec;
  double band;
  double sigma;
  bool plus;
public:
  GaussianPrefKernel(Eigen::VectorXd vec_, double band_, double sigma_, bool plus_) : vec(vec_), band(band_), sigma(sigma_), plus(plus_) {}
  
  double operator()(const double& x) const
  {
    //Calculate kernel
    const double pinum = 3.14159265359;
    const double pi2 = 1.0/std::sqrt(2*pinum);

    // Kernel Function
    double K = (-(vec.array()-x).square().array()/(2*band*band)).exp().sum();
    K = K * pi2;
    K = (1/(vec.size()*band)) * K;
    
    //  Preference function
    double res = 0;
    for(int j = 0; j < vec.size(); j++){
      for(int i = 0; i < vec.size(); i++){
       if(plus){
         if(vec(j) <= vec(i)){
           double qq = 1 - std::exp(-(std::pow((vec(i)-vec(j)), 2))/(2*sigma));
           res = res + (qq * K);
         }
       }
       else{
         if(vec(j) >= vec(i)){
           double qq = 1 - std::exp(-(std::pow((vec(i)-vec(j)), 2))/(2*sigma));
           res = res + (qq * K);
         }
       }
      }
    }
    return(res);
  }
};


// [[Rcpp::export]]
Rcpp::List integrate_test_RMC()
{
  const double lower = 4.3, upper = 6.7;
  Eigen::Vector3d vec(5.2, 4.3, 6.7);
  double band = 0.5*0.5;
  double sigma = 0.7*0.7;
  bool plus = true;
  
  GaussianPrefKernel f(vec, band, sigma, plus);
  
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}