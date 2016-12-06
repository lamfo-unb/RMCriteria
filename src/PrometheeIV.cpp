#include <RcppEigen.h>
#include "Preference.h"
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// P(0.3 < X < 0.8), X ~ Beta(a, b)
class BetaPDF: public Func {
private:
  double a;
  double b;
public:
  BetaPDF(double a_, double b_) : a(a_), b(b_) {}

  double operator()(const double& x) const
  {
    return R::dbeta(x, a, b, 0);
  }
};

// [[Rcpp::export]]
Rcpp::List integrate_test()
{
  const double a = 3, b = 10;
  const double lower = 0.3, upper = 0.8;
  const double true_val = R::pbeta(upper, a, b, 1, 0) -
    R::pbeta(lower, a, b, 1, 0);

  BetaPDF f(a, b);
  double err_est;
  int err_code;
  const double res = Numer::integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("true") = true_val,
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}


// Create the Kernel matrix
// @param datVec  Column of the dataset
// @param int Type od preference function
// @return Preference Matrix
Eigen::MatrixXd matPrometheeIV(Eigen::VectorXd datVec,int prefFunction, Eigen::VectorXd parms){
  //Get the number of rows
  int rows=datVec.size();

  //Initialize the matPromethee
  Eigen::MatrixXd matPromethee = Eigen::MatrixXd::Zero(rows, rows);

  //Variance
  double sum2 = 0;
  double cont = 0;
  //Initialize the computation
  for(int i=0;i<rows;i++){
    for(int j=i+1;j<rows;j++){
      //Compute the difference
      double delta = datVec(i)-datVec(j);
      matPromethee(i,j)=delta;
      matPromethee(j,i)=(-1.0)*delta;
      sum2=sum2+2.0*std::pow(matPromethee(j,i),2);
      cont=cont+1;
    }
  }

  //Calculate the sigma
  double sigma = (sum2/(2*cont));

  return(matPromethee);
}


// Create the Kernel matrix
// @param datVec  Column of the dataset
// @param int Type od preference function
// @return Preference Matrix
// [[Rcpp::export]]
Rcpp::List PrometheeIV(Eigen::MatrixXd datMat, Eigen::VectorXd vecWeights, Eigen::VectorXi prefFunction,Eigen::VectorXi alphaVector, Eigen::MatrixXd parms){
  //Get the number of rows
  int rows=datMat.rows();

  //Get the number of columns
  int cols=datMat.cols();

  //Initialize the matPromethee
  Eigen::MatrixXd matPromethee = Eigen::MatrixXd::Zero(rows, rows);





  //Create the Flow vector
  Eigen::VectorXd limInf  = Eigen::VectorXd::Zero(rows);
  Eigen::VectorXd limSup = Eigen::VectorXd::Zero(rows);

  //'Store the results
  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("limInf") = limInf,
                                    Rcpp::Named("limSup")  = limSup);

  return (resTemp);
}

