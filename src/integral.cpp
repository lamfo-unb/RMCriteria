// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;


class GaussianPrefKernel: public Func {
private:
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

/////////////////////////////////////////////////////
////////////////////////////////////////////////////

class UsualPrefKernel: public Func {
private:
  Eigen::VectorXd vec;
  double band;
  bool plus;
  double alt;
public:
  UsualPrefKernel(Eigen::VectorXd vec_, double band_, bool plus_, double alt_) : vec(vec_), band(band_), plus(plus_), alt(alt_) {}

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
      double delta = vec(j) - alt;
        if(plus){
          if(delta >= 0){
            double qq = 1;
            res = res + (qq * K);
          } else {
            double qq = 0;
            res = res + (qq * K);
          }
        } else{
          if(delta < 0){
            double qq = 1;
            res = res + (qq * K);
          } else {
            double qq = 0;
            res = res + (qq * K);
          }
    }
  }
    return(res);
  }
};

/////////////////////////////////////////////////////
////////////////////////////////////////////////////

class UShapePrefKernel: public Func {
private:
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

/////////////////////////////////////////////////////
////////////////////////////////////////////////////

class LevelPrefKernel: public Func {
private:
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

/////////////////////////////////////////////////////
////////////////////////////////////////////////////

class VShapePrefKernel: public Func {
private:
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

/////////////////////////////////////////////////////
////////////////////////////////////////////////////

class VShapeIndPrefKernel: public Func {
private:
  Eigen::VectorXd vec;
  double band;
  bool plus;
  double q;
  double p;

public:
  VShapeIndPrefKernel(Eigen::VectorXd vec_, double band_, bool plus_, double q__, double p__) : vec(vec_), band(band_), plus(plus_), q(q__), p(p__) {}

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
        double deltaji = vec(j) - vec(i);
        if(plus){
          if(deltaji <= q){
            double qq = 0.0;
            res = res + (qq * K);
          } else if(deltaji <= p){
            double qq = (deltaji-q)/(p-q);
            res = res + (qq * K);
          } else{
            double qq = 1.0;
            res = res + (qq * K);
          }
        } else{
          if(deltaji >= q){
            double qq = 0.0;
            res = res + (qq * K);
          } else if(deltaji >= p){
            double qq = (deltaji-q)/(p-q);
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

/////////////////////////////////////////////////////
////////////////////////////////////////////////////

// [[Rcpp::export]]
Eigen::VectorXd integrate_KernelPromethee(Eigen::MatrixXd dados, int prefFunction, Eigen::VectorXd weights, Eigen::VectorXd parms, double band, bool normalize)
{
  int colunas = dados.cols();
  Eigen::VectorXd sol(dados.rows());
  Eigen::VectorXd phi(dados.rows());

  int c = 0;
  for(c = 0; c < colunas; c++){

    double lower = dados.colwise().minCoeff()(c);
    double upper = dados.colwise().maxCoeff()(c);
    //    std::cout << "lower " << c << " " << lower << std::endl;
    //    std::cout << "upper " << c << " " << upper << std::endl;
    Eigen::VectorXd vec = dados.col(c);
    //    std::cout << "vec " << vec << std::endl;
    //    std::cout << "dados.rows " << dados.rows() << std::endl;
    bool plus = true;
    //    std::cout << "dados " << dados.col(c) << std::endl;

    for(int i = 0; i < dados.size(); i++){
      double err_est;
      int err_code;
      if (prefFunction == 0)
      {
//        VectorXd vec_, double band_, double sigma_, bool plus_
        GaussianPrefKernel f(vec, band, parms(1), plus);
        const double res = integrate(f, lower, upper, err_est, err_code);
        sol(i) = sol(i) + res;
        std::cout << "res " << res << std::endl;
      }
      else if (prefFunction == 1)
      {
        UsualPrefKernel f(vec, band, plus, vec(i));
        const double res = integrate(f, lower, upper, err_est, err_code);
        sol(i) = sol(i) + res;
        std::cout << "res " << res << std::endl;
      }
      else if (prefFunction == 2)
      {
        UShapePrefKernel f(vec, band, plus, vec(i));
        const double res = integrate(f, lower, upper, err_est, err_code);
        sol(i) = sol(i) + res;
        std::cout << "res " << res << std::endl;
      }

      else if (prefFunction == 3)
      {
        VShapePrefKernel f(vec, band, plus, vec(i));
        const double res = integrate(f, lower, upper, err_est, err_code);
        sol(i) = sol(i) + res;
        std::cout << "res " << res << std::endl;
      }
      // else if (prefFunction == 4)
      // {
      //   LevelPrefKernel f(vec, band, plus, vec(i));
      //   const double res = integrate(f, lower, upper, err_est, err_code);
      //   sol(i) = sol(i) + res;
      //   std::cout << "res " << res << std::endl;
      // }
      // else if (prefFunction == 5)
      // {
      //   VShapeIndPrefKernel f(vec, band, plus, vec(i));
      //   const double res = integrate(f, lower, upper, err_est, err_code);
      //   sol(i) = sol(i) + res;
      //   std::cout << "res " << res << std::endl;
      // }
    }
  }

  phi = phi.array() + sol.array() * weights(c);
  phi = phi.array()/weights.sum();
  double min = phi.minCoeff();
  double max = phi.maxCoeff();
  if(normalize == true) phi = (phi.array() - min)/(max - min);

  return(phi);
}
