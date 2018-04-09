#include <RcppEigen.h>
#include "Preference.h"
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

class Usual: public Func {
private:
  double y;
  bool plus;
public:
  Usual(double y_, bool plus_): y(y_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=0){
        return(0.0);
      }
      else
      {
        return(1.0);
      }
    }
    else{
      if((x-y)<=0){
        return(0.0);
      }
      else
      {
        return(1.0);
      }
    }
  }
};


class UShape: public Func {
private:
  double y;
  double q;
  bool plus;
public:
  UShape(double y_, double q_, bool plus_) :  y(y_), q(q_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=q){
        return(0.0);
      }
      else
      {
        return(1.0);
      }
    }
    else{
      if((x-y)<=q){
        return(0.0);
      }
      else
      {
        return(1.0);
      }
    }
  }
};

class VShape: public Func {
private:
  double y;
  double p;
  bool plus;
public:
  VShape(double y_, double p_, bool plus_) : y(y_), p(p_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=0){
        return(0.0);
      }
      else if((y-x)<=p){
        return((y-x)/p);
      }
      else
      {
        return(1.0);
      }
    }
    else{
      if((x-y)<=0){
        return(0.0);
      }
      else if((x-y)<=p){
        return((x-y)/p);
      }
      else
      {
        return(1.0);
      }
    }
  }
};

class Level: public Func {
private:
  double y;
  double q;
  double p;
  bool plus;
public:
  Level(double y_, double q_, double p_, bool plus_) : y(y_), q(q_), p(p_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=q){
        return(0.0);
      }
      else if((y-x)<=p){
        return(0.5);
      }
      else
      {
        return(1.0);
      }
    }
    else{
      if((x-y)<=q){
        return(0.0);
      }
      else if((x-y)<=p){
        return(0.5);
      }
      else
      {
        return(1.0);
      }
    }
  }
};

class VShapeIndPref: public Func {
private:
  double y;
  double q;
  double p;
  bool plus;
public:
  VShapeIndPref(double y_, double q_, double p_, bool plus_) : y(y_), q(q_), p(p_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=q){
        return(0.0);
      }
      else if((y-x)<=p){
        return(((y-x)-q)/(p-q));
      }
      else
      {
        return(1.0);
      }
    }
    else{
      if((x-y)<=q){
        return(0.0);
      }
      else if((x-y)<=p){
        return(((x-y)-q)/(p-q));
      }
      else
      {
        return(1.0);
      }
    }
  }
};


class GaussianPref: public Func {
private:
  double y;
  double sigma;
  bool plus;
public:
  GaussianPref(double y_, double sigma_, bool plus_) : y(y_), sigma(sigma_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus==true){
      return(1.0-std::exp(-1.0*std::pow(y-x,2)/(2*std::pow(sigma,2))));
    }
    else{
      return(-1.0+std::exp(-1.0*std::pow(y-x,2)/(2*std::pow(sigma,2))));
    }
    return(0);
  }
};

//' Calculates PROMETHEE IV method.
//'
//' @param datMat A matrix containing the data from criterias and alternatives.
//' @param vecWeights A vector of weights for each criteria.
//' @param prefFunction A numerical vector to indicate the type of the
//' Preference Function:
//'   \itemize{
//'     \item \code{prefFunction = 0} Gaussian Preference Function
//'     \item \code{prefFunction = 1} Usual Preference Function
//'     \item \code{prefFunction = 2} U-Shape Preference Function
//'     \item \code{prefFunction = 3} V-Shape Preference Function
//'     \item \code{prefFunction = 4} Level Preference Function
//'     \item \code{prefFunction = 5} V-Shape Preference and Indiference Function
//'     }
//' @param parms A numerical matrix with parameters associated to the Preference
//'  Function. They're defined as a matrix of n columns and m rows. The maximum
//'  number of parameters is 3 and m is the number of criterias. The parameters
//'  are:
//'   \itemize{
//'   \item{Indifference Threshold (\code{q})}
//'   \item{Preference Threshold (\code{p})}
//'   \item{Gaussian Threshold (\code{s})}
//'   }
//' @param normalize A boolean to normalize the index.
//' @return Preference Matrix
//' @export
// [[Rcpp::export]]


Rcpp::List PrometheeIV(Eigen::MatrixXd datMat, Eigen::VectorXd vecWeights, Eigen::VectorXi prefFunction, Eigen::MatrixXd parms, bool normalize){
  //Get the number of rows
  int rows=datMat.rows();

  //Get the number of columns
  int cols=datMat.cols();

  //Create the Flow vector
  Eigen::VectorXd phiPlusVec  = Eigen::VectorXd::Zero(rows);
  Eigen::VectorXd phiMinusVec = Eigen::VectorXd::Zero(rows);

  //Get maximum
  Eigen::VectorXd maxVec = datMat.colwise().maxCoeff();
  Eigen::VectorXd minVec = datMat.colwise().minCoeff();


  //Foreach alternative
  for(int r=0;r<rows;r++){
    //Phi+
    double phiPlus=0.0;
    double phiMinus=0.0;
     //Foreach criterion
    for(int c=0;c<cols;c++){
      //Gaussian
      if(prefFunction(c)==0){
        Eigen::VectorXd sigmaVec = parms.row(c);
        GaussianPref fplus(datMat(r,c),sigmaVec(0),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        GaussianPref fminus(datMat(r,c),sigmaVec(0),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==1){
        Usual fplus(datMat(r,c),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        Usual fminus(datMat(r,c),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==2){
        Eigen::VectorXd parmsVec = parms.row(c);
        UShape fplus(datMat(r,c),parmsVec(0),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        UShape fminus(datMat(r,c),parmsVec(0),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==3){
        Eigen::VectorXd parmsVec = parms.row(c);
        VShape fplus(datMat(r,c),parmsVec(0),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        VShape fminus(datMat(r,c),parmsVec(0),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==4){
        Eigen::VectorXd parmsVec = parms.row(c);
        Level fplus(datMat(r,c),parmsVec(0),parmsVec(1),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        Level fminus(datMat(r,c),parmsVec(0),parmsVec(1),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==5){
        Eigen::VectorXd parmsVec = parms.row(c);
        VShapeIndPref fplus(datMat(r,c),parmsVec(0),parmsVec(1),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        VShapeIndPref fminus(datMat(r,c),parmsVec(0),parmsVec(1),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
    }

    phiPlusVec(r) = phiPlus;
    phiMinusVec(r) =phiMinus;
  }

  //Index
  Eigen::VectorXd phi = phiPlusVec-phiMinusVec;

  //Normalize
  double min = phi.minCoeff();
  double max = phi.maxCoeff();
  if(normalize==true) phi = (phi.array() - min)/(max-min);

  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("Phi+") = phiPlusVec,
                                          Rcpp::Named("Phi-")  = phiMinusVec,
                                          Rcpp::Named("Index")  = phi);

  return (resTemp);
}

