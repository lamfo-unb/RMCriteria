#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// Finalizar fun??es de prefer?ncia, pensar no PhiMinus e encapsulamento no R (fun??o S4 que chame o Rcpp)

// [[Rcpp::export]]
double GaussianPrefKernel(double y,Eigen::VectorXd vec, double band, double sigma, bool plus) {
    //Calculate kernel
    const double pinum = 3.14159265359;
    const double pi2 = 1.0/std::sqrt(2*pinum);
    if(plus == false) vec = -vec;
// Kernel Function
    double K = (-(vec.array()-y).square().array()/(2*band*band)).exp().sum();
    K = K * pi2;
    K = (1/(vec.size()*band)) * K;

//  Preference function
    double res = 0;
    for(int j = 0; j < vec.size(); j++){
      for(int i = 0; i < vec.size(); i++){
        if(vec(j) <= vec(i)){
          double qq = 1 - std::exp(-(std::pow((vec(i)-vec(j)), 2))/(2*sigma));
          res = res + (qq * K);
        }
      }
    }

    return(res);
  }



// [[Rcpp::export]]
double UsualPrefKernel(double y,Eigen::VectorXd vec, double band, bool plus) {
  //Calculate kernel
  const double pinum = 3.14159265359;
  const double pi2 = 1.0/std::sqrt(2*pinum);
  if(plus == false) vec = -vec;
  // Kernel Function
  double K = (-(vec.array()-y).square().array()/(2*band*band)).exp().sum();
  K = K * pi2;
  K = (1/(vec.size()*band)) * K;

  //  Preference function
  double res = 0;
  for(int j = 0; j < vec.size(); j++){
    for(int i = 0; i < vec.size(); i++){
      if(vec(j) <= vec(i)){
        double qq = 1;
        res = res + (qq * K);
      }
    }
  }
  return(res);
}


// [[Rcpp::export]]
double UShapePrefKernel(double y,Eigen::VectorXd vec, double band, bool plus, double q) {
  //Calculate kernel
  const double pinum = 3.14159265359;
  const double pi2 = 1.0/std::sqrt(2*pinum);
  if(plus == false) vec = -vec;
  // Kernel Function
  double K = (-(vec.array()-y).square().array()/(2*band*band)).exp().sum();
  K = K * pi2;
  K = (1/(vec.size()*band)) * K;

  //  Preference function
  double res = 0;
  for(int j = 0; j < vec.size(); j++){
    for(int i = 0; i < vec.size(); i++){
      double deltaji = vec(j) - vec(i);
      if(deltaji >= q){
        double qq = 1;
        res = res + (qq * K);
      }
    }
  }
  return(res);
}


// [[Rcpp::export]]
double LevelPrefKernel(double y,Eigen::VectorXd vec, double band, bool plus, double q, double p) {
  //Calculate kernel
  const double pinum = 3.14159265359;
  const double pi2 = 1.0/std::sqrt(2*pinum);
  if(plus == false) vec = -vec;
  // Kernel Function
  double K = (-(vec.array()-y).square().array()/(2*band*band)).exp().sum();
  K = K * pi2;
  K = (1/(vec.size()*band)) * K;

  //  Preference function
  double res = 0;
  for(int j = 0; j < vec.size(); j++){
    for(int i = 0; i < vec.size(); i++){
      double deltaji = vec(j) - vec(i);
      if(deltaji >= p){
        double qq = 1;
        res = res + (qq * K);
      }
      else if(deltaji >= q){
        double qq = 0.5;
        res = res + (qq * K);
      }
    }
  }
  return(res);
}

// [[Rcpp::export]]
double VShapePrefKernel(double y,Eigen::VectorXd vec, double band, bool plus, double p) {
  //Calculate kernel
  const double pinum = 3.14159265359;
  const double pi2 = 1.0/std::sqrt(2*pinum);
  if(plus == false) vec = -vec;
  // Kernel Function
  double K = (-(vec.array()-y).square().array()/(2*band*band)).exp().sum();
  K = K * pi2;
  K = (1/(vec.size()*band)) * K;

  //  Preference function
  double res = 0.0;
  for(int j = 0.0; j < vec.size(); j++){
    for(int i = 0.0; i < vec.size(); i++){
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
    }
  }
  return(res);
}



// [[Rcpp::export]]
double VShapeIndPrefKernel(double y,Eigen::VectorXd vec, double band, bool plus, double q, double p) {
  //Calculate kernel
  const double pinum = 3.14159265359;
  const double pi2 = 1.0/std::sqrt(2*pinum);
  if(plus == false) vec = -vec;
  // Kernel Function
  double K = (-(vec.array()-y).square().array()/(2*band*band)).exp().sum();
  K = K * pi2;
  K = (1/(vec.size()*band)) * K;

  //  Preference function
  double res = 0.0;
  for(int j = 0.0; j < vec.size(); j++){
    for(int i = 0.0; i < vec.size(); i++){
      double deltaji = vec(j) - vec(i);
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
    }
  }
  return(res);
}




/*** R
library(RcppEigen)
library(RcppNumerical)
library(Rcpp)
vec <- c(5.2, 4.3, 6.7)
band <- 0.5^2
point <- 2
sigma <- 0.7^2
q <- 0.4
p <- 2


K <- (1/(length(vec)*band)) * (1/(sqrt(2*pi))) * sum(exp(-(vec-2)^2/(2*band*band)))

######################
# Gaussian Preference

res <- 0
for(j in 1:length(vec)){
  for(i in 1:length(vec)){
    if(vec[j] <= vec[i]){
      qq <- (1 - exp(-((vec[i]-vec[j])^2)/(2*sigma)))
      res <- res + (qq * K)
    }
  }
}

paste("Gaussian Rcpp: ", GaussianPrefKernel(point, vec, band, sigma, TRUE))
paste("Gaussian R", res)

######################
# Usual Preference

res <- 0
for(j in 1:length(vec)){
  for(i in 1:length(vec)){
    if(vec[j] <= vec[i]){
      qq <- 1
      res <- res + (qq * K)
    }
  }
}

paste("Usual Rcpp", UsualPrefKernel(point, vec, band, TRUE))
paste("Usual R", res)

######################
# U Shape Preference

res <- 0
for(j in 1:length(vec)){
  for(i in 1:length(vec)){
    deltaji <- vec[j] - vec[i]
    if(deltaji >= q){
      qq <- 1.0
      res <- res + (qq * K)
    }
  }
}

paste("UShape Rcpp: ", UShapePrefKernel(point, vec, band, TRUE, q))
paste("UShape R", res)

######################
# Level Preference

res <- 0
for(j in 1:length(vec)){
  for(i in 1:length(vec)){
    deltaji <- vec[j] - vec[i]
    if(deltaji >= p){
      qq <- 1.0
      res <- res + (qq * K)
    }
    else if(deltaji >= q){
      qq <- 0.5 * K
      res <- res + (qq * K)
    }
  }
}

paste("Level Rcpp: ", LevelPrefKernel(point, vec, band, TRUE, q, p))
paste("Level R", res)


######################
# V Shape Preference

res <- 0
for(j in 1:length(vec)){
  for(i in 1:length(vec)){
    deltaji <- vec[j] - vec[i]
    if(deltaji <= 0){
      qq <- 0.0
      res <- res + (qq * K)
    } else if(deltaji <= p){
      qq <- deltaji/p
      res <- res + (qq * K)
    } else {
      qq <- 1
      res <- res + (qq * K)
    }
  }
}


paste("V Shape Rcpp: ", VShapePrefKernel(point, vec, band, TRUE, p))
paste("V Shape R", res)



######################
# V Shape Indifference Preference

res <- 0
for(j in 1:length(vec)){
  for(i in 1:length(vec)){
    deltaji <- vec[j] - vec[i]
    if(deltaji <= q){
      qq <- 0.0
      res <- res + (qq * K)
    } else if(deltaji <= p){
      qq <- (deltaji-q)/(p-q)
      res <- res + (qq * K)
    } else {
      qq <- 1
      res <- res + (qq * K)
    }
  }
}


paste("V Shape Ind Rcpp: ", VShapeIndPrefKernel(point, vec, band, TRUE, q, p))
paste("V Shape Ind R", res)

*/
