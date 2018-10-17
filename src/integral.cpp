#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#include <math.h>
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

    //Get the number of rows
    int rows = vec.size();

    //Initialize the matDelta
    Eigen::MatrixXd matDelta = Eigen::MatrixXd::Zero(rows, rows);


    //Fill the matDelta
    for(int i = 0; i < rows; i++){
      for(int j = i + 1; j < rows; j++){
        double delta = vec(i) - vec(j);
        matDelta(i, j) = delta;
        matDelta(j, i) = (-1.0)*delta;
      }
    }


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

class UsualPrefKernelPlus: public Func {
private:
  Eigen::VectorXd vec;
  double band;
  bool plus;
  int alt;
public:
  UsualPrefKernelPlus(Eigen::VectorXd vec_, double band_, bool plus_, int alt_) : vec(vec_), band(band_), plus(plus_), alt(alt_) {}

  double operator()(const double& x) const
  {
    // Band para cada critério, phi menos, testes, pref functions e colunas

    //Get the number of rows
    int rows = vec.size();


    //Calculate kernel
    const double pinum = 3.14159265359;
    const double pi2 = 1.0/std::sqrt(2*pinum);
    // Kernel Function
    double K = 0;
    double res = 0;
    // Melhoria: tirar o cálculo do K de dentro da função
    for(int i = 0; i < rows; i++){
      K = (1/(vec.size()*band)) * pi2 * exp(-(1/2)*pow((((vec(alt)-vec(i)) - (vec(alt)-x))/band), 2));

      // Preference function using delta of alternative with x
      double delta = vec(alt) - x;
      if(plus){
        if(delta <= 0){
          double qq = 0;
          res = res + K * qq;
        } else {
          double qq = 1;
          res = res + K * qq;
        }
      } else{
        if(delta >= 0){
          double qq = 0;
          res = res + K * qq;
        } else {
          double qq = 1;
          res = res + K * qq;
        }
      }
    }


    //Create the Flow vector
    // Eigen::VectorXd phiPlus  = Eigen::VectorXd::Zero(rows);
    // Eigen::VectorXd phiMinus = Eigen::VectorXd::Zero(rows);
    //
    // for(int row = 0; row < rows; row++){
    //   phiPlus(row) = matDelta.row(row).sum();
    //   phiMinus(row) = matDelta.col(row).sum();
    // }
//    std::cout << "K = " << K << std::endl;
//    std::cout << "vec(alt) " << alt << " = " << vec(alt) << std::endl;
    return(res);
  }
};




// [[Rcpp::export]]
double integrate_KernelPrometheePlus(Eigen::VectorXd dados, int prefFunction, Eigen::VectorXd weights, Eigen::VectorXd parms, double band, bool normalize, int alt)
{
//  int colunas = dados.cols();
//  Eigen::VectorXd sol(dados.rows());
//  double sol = 0;
  Eigen::VectorXd phi(dados.rows());

//  int c = 0;
//  for(c = 0; c < colunas; c++){

//    double lower = dados.colwise().minCoeff()(c);
//    double upper = dados.colwise().maxCoeff()(c);
    double lower = dados.minCoeff();
    double upper = dados.maxCoeff();
    double res = 0;
//    Eigen::VectorXd vec = dados.col(c);
    Eigen::VectorXd vec = dados;
     bool plus = true;

     // Variables to be used by integral function
      double err_est;
      int err_code;
      if (prefFunction == 0){
        GaussianPrefKernel f(vec, band, parms(1), plus);
        res = integrate(f, lower, upper, err_est, err_code);
      }
      else if (prefFunction == 1){
        UsualPrefKernelPlus f(vec, band, plus, vec(alt));
        res = integrate(f, lower, upper, err_est, err_code);
      }
      // else if (prefFunction == 2)
      // {
      //   UShapePrefKernel f(vec, band, plus, vec(i));
      //   const double res = integrate(f, lower, upper, err_est, err_code);
      //   sol(i) = sol(i) + res;
      //   std::cout << "res " << res << std::endl;
      // }
      //
      // else if (prefFunction == 3)
      // {
      //   VShapePrefKernel f(vec, band, plus, vec(i));
      //   const double res = integrate(f, lower, upper, err_est, err_code);
      //   sol(i) = sol(i) + res;
      //   std::cout << "res " << res << std::endl;
      // }
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
//    }


 // phi = phi.array() + sol.array() * weights(c);
 // phi = phi.array()/weights.sum();
 // double min = phi.minCoeff();
 // double max = phi.maxCoeff();
 // if(normalize == true) phi = (phi.array() - min)/(max - min);

  return(res);
}





/////////////////////////////////////////////////////
////////////////////////////////////////////////////

class UsualPrefKernelMinus: public Func {
private:
  Eigen::VectorXd vec;
  double band;
  bool plus;
  int alt;
public:
  UsualPrefKernelMinus(Eigen::VectorXd vec_, double band_, bool plus_, int alt_) : vec(vec_), band(band_), plus(plus_), alt(alt_) {}

  double operator()(const double& x) const
  {
    // Band para cada critério, phi menos, testes, pref functions e colunas

    //Get the number of rows
    int rows = vec.size();


    //Calculate kernel
    const double pinum = 3.14159265359;
    const double pi2 = 1.0/std::sqrt(2*pinum);
    // Kernel Function
    double K = 0;
    double res = 0;
    // Melhoria: tirar o cálculo do K de dentro da função
    for(int i = 0; i < rows; i++){
      K = (1/(vec.size()*band)) * pi2 * exp(-(1/2)*pow((((vec(alt)-vec(i)) - (vec(alt)-x))/band), 2));

      // Preference function using delta of alternative with x
      double delta = x - vec(alt);
      if(plus){
        if(delta <= 0){
          double qq = 0;
          res = res + K * qq;
        } else {
          double qq = 1;
          res = res + K * qq;
        }
      } else{
        if(delta >= 0){
          double qq = 0;
          res = res + K * qq;
        } else {
          double qq = 1;
          res = res + K * qq;
        }
      }
    }


    //Create the Flow vector
    // Eigen::VectorXd phiPlus  = Eigen::VectorXd::Zero(rows);
    // Eigen::VectorXd phiMinus = Eigen::VectorXd::Zero(rows);
    //
    // for(int row = 0; row < rows; row++){
    //   phiPlus(row) = matDelta.row(row).sum();
    //   phiMinus(row) = matDelta.col(row).sum();
    // }
    //    std::cout << "K = " << K << std::endl;
    //    std::cout << "vec(alt) " << alt << " = " << vec(alt) << std::endl;
    return(res);
  }
};




// [[Rcpp::export]]
double integrate_KernelPrometheeMinus(Eigen::VectorXd dados, int prefFunction, Eigen::VectorXd weights, Eigen::VectorXd parms, double band, bool normalize, int alt)
{
  //  int colunas = dados.cols();
  //  Eigen::VectorXd sol(dados.rows());
  //  double sol = 0;
  Eigen::VectorXd phi(dados.rows());

  //  int c = 0;
  //  for(c = 0; c < colunas; c++){

  //    double lower = dados.colwise().minCoeff()(c);
  //    double upper = dados.colwise().maxCoeff()(c);
  double lower = dados.minCoeff();
  double upper = dados.maxCoeff();
  double res = 0;
  //    Eigen::VectorXd vec = dados.col(c);
  Eigen::VectorXd vec = dados;
  bool plus = true;

  // Variables to be used by integral function
  double err_est;
  int err_code;
  if (prefFunction == 0){
    GaussianPrefKernel f(vec, band, parms(1), plus);
    res = integrate(f, lower, upper, err_est, err_code);
  }
  else if (prefFunction == 1){
    UsualPrefKernelMinus f(vec, band, plus, vec(alt));
    res = integrate(f, lower, upper, err_est, err_code);
  }
  // else if (prefFunction == 2)
  // {
  //   UShapePrefKernel f(vec, band, plus, vec(i));
  //   const double res = integrate(f, lower, upper, err_est, err_code);
  //   sol(i) = sol(i) + res;
  //   std::cout << "res " << res << std::endl;
  // }
  //
  // else if (prefFunction == 3)
  // {
  //   VShapePrefKernel f(vec, band, plus, vec(i));
  //   const double res = integrate(f, lower, upper, err_est, err_code);
  //   sol(i) = sol(i) + res;
  //   std::cout << "res " << res << std::endl;
  // }
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
  //    }


  // phi = phi.array() + sol.array() * weights(c);
  // phi = phi.array()/weights.sum();
  // double min = phi.minCoeff();
  // double max = phi.maxCoeff();
  // if(normalize == true) phi = (phi.array() - min)/(max - min);

  return(res);
}




























// [[Rcpp::export]]
double Ktest(Eigen::VectorXd vec, double band, bool plus, int alt)
{
  //Get the number of rows
  int rows = vec.size();

  //Calculate kernel
  const double pinum = 3.14159265359;
  const double pi2 = 1.0/std::sqrt(2*pinum);
  // Kernel Function
  double K = 0;
  double res = 0;
  double x = 5;
  for(int i = 0; i < rows; i++){
    K = (1/(vec.size()*band)) * pi2 * exp(-(1/2)*pow((((vec(alt)-vec(i)) - (vec(alt)-x))/band), 2));

    // Preference function using delta of alternative with x
    double delta = vec(alt) - x;
    if(plus){
      if(delta <= 0){
        double qq = 0;
        res = res + K * qq;
      } else {
        double qq = 1;
        res = res + K * qq;
      }
    } else{
      if(delta >= 0){
        double qq = 0;
        res = res + K * qq;
      } else {
        double qq = 1;
        res = res + K * qq;
      }
    }
  }


  //Create the Flow vector
  // Eigen::VectorXd phiPlus  = Eigen::VectorXd::Zero(rows);
  // Eigen::VectorXd phiMinus = Eigen::VectorXd::Zero(rows);
  //
  // for(int row = 0; row < rows; row++){
  //   phiPlus(row) = matDelta.row(row).sum();
  //   phiMinus(row) = matDelta.col(row).sum();
  // }

  return(res);

  return(res);
}









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


