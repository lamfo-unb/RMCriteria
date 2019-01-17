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
library(RMCriteria)
vec <- matrix(c(5.2, 4.3, 6.7), ncol = 1)
vec <- matrix(c(5.2, 4.3, 6.7,
                4.4, 1.3, 2.8), ncol = 2)
weights <- c(0.3, 0.7)
band <- 0.5^2
normalize  <-  FALSE
parms <- matrix(c(NA, NA), byrow = TRUE, nrow = 2)
# dados, prefFunction, weights, parms, band, normalize)
integrate_KernelPromethee(dados = vec, prefFunction = 1, weights = weights, parms = parms, band = band, normalize = normalize, alt = 0)
integrate_KernelPromethee(dados = vec, prefFunction = 1, weights = weights, parms = parms, band = band, normalize = normalize, alt = 1)
integrate_KernelPromethee(dados = vec, prefFunction = 1, weights = weights, parms = parms, band = band, normalize = normalize, alt = 2)

Ktest(vec, band, plus = TRUE, alt = 0)
Ktest(vec, band, plus = TRUE, alt = 1)
Ktest(vec, band, plus = TRUE, alt = 2)

point <- 2
sigma <- 0.7^2
q <- 0.4
p <- 2


K <- (1/(length(vec)*band)) * (1/(sqrt(2*pi))) * sum(exp(-(vec-2)^2/(2*band*band)))
#K = (-(vec.array()-x).square().array()/(2*band*band)).exp().sum();

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


res <- rep(NA, length(vec))
for(i in 1:length(vec)){
  lista <- integrate_UsualPref(vec[i])
  res[i] <- lista[[1]]
}

vec.df <- as.data.frame(vec)
res <- apply(vec, 1, function(x){
  lista <- integrate_UsualPref(x)
  lista[[1]]
})


paste("Usual Rcpp", UsualPrefKernel(point, vec, band, TRUE))
paste("Usual R", res)
paste("K ", K)

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
paste("UShape K", K)

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
