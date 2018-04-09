#include <RcppEigen.h>
#include "Preference.h"
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

class UsualKernel: public Func {
private:
  double y;
  Eigen::VectorXd vec;
  double band;
  bool plus;
public:
  UsualKernel(double y_, Eigen::VectorXd vec_, double band_, bool plus_): y(y_), vec(vec_), band(band_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=0){
        return(0.0);
      }
      else
      {
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*1.0);
      }
    }
    else{
      if((x-y)<=0){
        return(0.0);
      }
      else
      {
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*1.0);
      }
    }
  }
};


class UShapeKernel: public Func {
private:
  double y;
  Eigen::VectorXd vec;
  double band;
  double q;
  bool plus;
public:
  UShapeKernel(double y_, Eigen::VectorXd vec_, double band_, double q_, bool plus_) : y(y_), vec(vec_), band(band_), q(q_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=q){
        return(0.0);
      }
      else
      {
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*1.0);
      }
    }
    else{
      if((x-y)<=q){
        return(0.0);
      }
      else
      {
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*1.0);
      }
    }
  }
};

class VShapeKernel: public Func {
private:
  double y;
  Eigen::VectorXd vec;
  double band;
  double p;
  bool plus;
public:
  VShapeKernel(double y_, Eigen::VectorXd vec_, double band_, double p_, bool plus_) : y(y_), vec(vec_), band(band_), p(p_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=0){
        return(0.0);
      }
      else if((y-x)<=p){
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*((y-x)/p));
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
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*((x-y)/p));
      }
      else
      {
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*1.0);
      }
    }
  }
};

class LevelKernel: public Func {
private:
  double y;
  Eigen::VectorXd vec;
  double band;
  double q;
  double p;
  bool plus;
public:
  LevelKernel(double y_, Eigen::VectorXd vec_, double band_, double q_, double p_, bool plus_) : y(y_),vec(vec_), band(band_), q(q_), p(p_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=q){
        return(0.0);
      }
      else if((y-x)<=p){
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*0.5);
      }
      else
      {
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*1.0);
      }
    }
    else{
      if((x-y)<=q){
        return(0.0);
      }
      else if((x-y)<=p){
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*0.5);
      }
      else
      {
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*1.0);
      }
    }
  }
};

class VShapeIndPrefKernel: public Func {
private:
  double y;
  Eigen::VectorXd vec;
  double band;
  double q;
  double p;
  bool plus;
  public:
  VShapeIndPrefKernel(double y_,Eigen::VectorXd vec_, double band_, double q_, double p_, bool plus_) : y(y_), vec(vec_), band(band_), q(q_), p(p_), plus(plus_) {}

  double operator()(const double& x) const
  {
    if(plus){
      if((y-x)<=q){
        return(0.0);
      }
      else if((y-x)<=p){
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*(((y-x)-q)/(p-q)));
      }
      else
      {
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*1.0);
      }
    }
    else{
      if((x-y)<=q){
        return(0.0);
      }
      else if((x-y)<=p){
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*(((x-y)-q)/(p-q)));
      }
      else
      {
        //Calculate kernel
        double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
        return(K*1.0);
      }
    }
  }
};


class GaussianPrefKernel: public Func {
private:
  double y;
  Eigen::VectorXd vec;
  double band;
  double sigma;
  bool plus;
public:
  GaussianPrefKernel(double y_,Eigen::VectorXd vec_, double band_, double sigma_, bool plus_) : y(y_),vec(vec_), band(band_), sigma(sigma_), plus(plus_) {}

  double operator()(const double& x) const
  {
    //Calculate kernel
    double K = ((-(vec.array()-x)/band).square().array()/2.0).exp().sum();
    if(plus==true){
      return(K*(1.0-std::exp(-1.0*std::pow(y-x,2)/(2*std::pow(sigma,2)))));
    }
    else{
      return(K*(-1.0+std::exp(-1.0*std::pow(y-x,2)/(2*std::pow(sigma,2)))));
    }
    return(0);
  }
};

//' Calculates PROMETHEE IV KERNEL method.
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
//' @param parms a numerical matrix with parameters associated to the Preference
//'  Function. They're defined as a matrix of n columns and m rows. The maximum
//'  number of parameters is 3 and m is the number of criterias. The parameters
//'  are:
//'   \itemize{
//'   \item{Indifference Threshold (\code{q})}
//'   \item{Preference Threshold (\code{p})}
//'   \item{Gaussian Threshold (\code{s})}
//'   }
//' @param band A numerical matrix with m rows corresponding to each criteria
//' and one column corresponding to the bandwitch estimated for that criteria.
//' This bandwitch is used for Kernel Density Estimation in Promethee IV Kernel.
//'  By default, it is calculated using bw.nrd0.
//' @param normalize A boolean to normalize the index.
//' @return Preference Matrix
//' @export
// [[Rcpp::export]]

Rcpp::List PrometheeIVKernel(Eigen::MatrixXd datMat, Eigen::VectorXd vecWeights, Eigen::VectorXi prefFunction, Eigen::MatrixXd parms, Eigen::MatrixXd band, bool normalize){
  //Get the number of rows
  int rows=datMat.rows();

  //Get the number of columns
  int cols=datMat.cols();

  //Create the Flow vector
  Eigen::VectorXd phiPlusVec  = Eigen::VectorXd::Zero(rows);
  Eigen::VectorXd phiMinusVec = Eigen::VectorXd::Zero(rows);

  //Get maximum
  Eigen::RowVectorXd mean = datMat.colwise().mean();
  Eigen::RowVectorXd std = ((datMat.rowwise() - mean).array().square().colwise().sum() / (datMat.rows() - 1)).sqrt();
  Eigen::RowVectorXd minVec = mean.array() - 6.0*std.array();
  Eigen::RowVectorXd maxVec = mean.array() + 6.0*std.array();

  //Foreach alternative
  for(int r=0;r<rows;r++){
    //Phi+
    double phiPlus=0.0;
    double phiMinus=0.0;
     //Foreach criterion
    for(int c=0;c<cols;c++){
      //Get the vector
      Eigen::VectorXd vec = datMat.col(c);

      //Gaussian
      if(prefFunction(c)==0){
        Eigen::VectorXd sigmaVec = parms.row(c);
        GaussianPrefKernel fplus(datMat(r,c), vec, band(c), sigmaVec(0),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        GaussianPrefKernel fminus(datMat(r,c), vec, band(c),sigmaVec(0),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==1){
        UsualKernel fplus(datMat(r,c), vec, band(c), true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        UsualKernel fminus(datMat(r,c), vec, band(c), false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==2){
        Eigen::VectorXd parmsVec = parms.row(c);
        UShapeKernel fplus(datMat(r,c), vec, band(c), parmsVec(0),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        UShapeKernel fminus(datMat(r,c), vec, band(c), parmsVec(0),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==3){
        Eigen::VectorXd parmsVec = parms.row(c);
        VShapeKernel fplus(datMat(r,c), vec, band(c), parmsVec(0),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        VShapeKernel fminus(datMat(r,c), vec, band(c), parmsVec(0),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==4){
        Eigen::VectorXd parmsVec = parms.row(c);
        LevelKernel fplus(datMat(r,c), vec, band(c), parmsVec(0),parmsVec(1),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        LevelKernel fminus(datMat(r,c), vec, band(c), parmsVec(0),parmsVec(1),false);
        double err_est_minus;
        int err_code_minus;
        const double resMinus = Numer::integrate(fminus, minVec(c), maxVec(c), err_est_minus, err_code_minus);
        phiMinus=phiMinus+resMinus*vecWeights(c);
      }
      else if(prefFunction(c)==5){
        Eigen::VectorXd parmsVec = parms.row(c);
        VShapeIndPrefKernel fplus(datMat(r,c), vec, band(c), parmsVec(0),parmsVec(1),true);
        double err_est_plus;
        int err_code_plus;
        const double resPlus = Numer::integrate(fplus, minVec(c), maxVec(c), err_est_plus, err_code_plus);
        phiPlus=phiPlus+resPlus*vecWeights(c);

        VShapeIndPrefKernel fminus(datMat(r,c), vec, band(c), parmsVec(0),parmsVec(1),false);
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

