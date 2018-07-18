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
Rcpp::List integrate_GaussianPref()
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



class UsualPrefKernel: public Func {
private:
  // double y;
  Eigen::VectorXd vec;
  double band;
  bool plus;
public:
  UsualPrefKernel(Eigen::VectorXd vec_, double band_, bool plus_) : vec(vec_), band(band_), plus(plus_) {}

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
            double qq = 1;
            res = res + (qq * K);
          }
        } else{
          if(vec(j) >= vec(i)){
            double qq = 1;
            res = res + (qq * K);
        }
      }
    }
  }
    return(res);
  }

};


// [[Rcpp::export]]
Rcpp::List integrate_UsualPref()
{
  const double lower = 4.3, upper = 6.7;
  Eigen::Vector3d vec(5.2, 4.3, 6.7);
  double band = 0.5*0.5;
  bool plus = true;

  UsualPrefKernel f(vec, band, plus);

  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}


class UShapePrefKernel: public Func {
private:
  // double y;
  Eigen::VectorXd vec;
  double band;
  bool plus;
  double q;
public:
  UShapePrefKernel(Eigen::VectorXd vec_, double band_, bool plus_, double q__) : vec(vec_), band(band_), plus(plus_), q(q__) {}

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
        double deltaji = vec(j) - vec(i);
        if(plus){
          if(deltaji >= q){
            double qq = 1;
            res = res + (qq * K);
          }
        } else{
          if(deltaji <= q){
            double qq = 1;
            res = res + (qq * K);
        }
      }
    }
    }
    return(res);
  }
};


// [[Rcpp::export]]
Rcpp::List integrate_UShapePref()
{
  const double lower = 4.3, upper = 6.7;
  Eigen::Vector3d vec(5.2, 4.3, 6.7);
  double band = 0.5*0.5;
  bool plus = true;
  double q = 0.4;

  UShapePrefKernel f(vec, band, plus, q);

  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}


class LevelPrefKernel: public Func {
private:
  // double y;
  Eigen::VectorXd vec;
  double band;
  bool plus;
  double q;
  double p;
public:
  LevelPrefKernel(Eigen::VectorXd vec_, double band_, bool plus_, double q__, double p__) : vec(vec_), band(band_), plus(plus_), q(q__), p(p__) {}

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
        double deltaji = vec(j) - vec(i);
        if(plus){
          if(deltaji >= p){
            double qq = 1;
            res = res + (qq * K);
          }
          else if(deltaji >= q){
            double qq = 0.5;
            res = res + (qq * K);
          }
        } else{
          if(deltaji <= p){
            double qq = 1;
            res = res + (qq * K);
          }
          else if(deltaji <= q){
            double qq = 0.5;
            res = res + (qq * K);
          }
        }
      }
    }
    return(res);
  }
};


// [[Rcpp::export]]
Rcpp::List integrate_LevelPref()
{
  const double lower = 4.3, upper = 6.7;
  Eigen::Vector3d vec(5.2, 4.3, 6.7);
  double band = 0.5*0.5;
  bool plus = true;
  double q = 0.4;
  double p = 2;

  LevelPrefKernel f(vec, band, plus, q, p);

  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}




class VShapePrefKernel: public Func {
private:
  // double y;
  Eigen::VectorXd vec;
  double band;
  bool plus;
  double p;
public:
  VShapePrefKernel(Eigen::VectorXd vec_, double band_, bool plus_, double p__) : vec(vec_), band(band_), plus(plus_), p(p__) {}

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
    for(int j = 0.0; j < vec.size(); j++){
      for(int i = 0.0; i < vec.size(); i++){
        if(plus){
          double deltaji = vec(j) - vec(i);
          if(deltaji <= 0.0){
            double qq = 0.0;
            res = res + (qq * K);
          } else if(deltaji <= p){
            double qq = deltaji/p;
            res = res + (qq * K);
          } else{
            double qq = 1.0;
            res = res + (qq * K);
          }
        } else{
          double deltaji = vec(j) - vec(i);
          if(deltaji >= 0.0){
            double qq = 0.0;
            res = res + (qq * K);
          } else if(deltaji >= p){
            double qq = deltaji/p;
            res = res + (qq * K);
          } else{
            double qq = 1.0;
            res = res + (qq * K);
          }
        }
      }
    }
    return(res);
  }
};


// [[Rcpp::export]]
Rcpp::List integrate_VShapePref()
{
  const double lower = 4.3, upper = 6.7;
  Eigen::Vector3d vec(5.2, 4.3, 6.7);
  double band = 0.5*0.5;
  bool plus = true;
  double p = 2;

  VShapePrefKernel f(vec, band, plus, p);

  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}

