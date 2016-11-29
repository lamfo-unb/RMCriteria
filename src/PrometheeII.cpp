#include <RcppEigen.h>
#include "Preference.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;


//Error Measure
double GaussianPreference(double delta, Eigen::VectorXd parms);


// Create the Kernel matrix
// @param datMat  Matrix with the data
// @param function Kernel Function
// @param parms vector of parameters fot the kernel
// @return Kernel Matrix
Eigen::MatrixXd KernelMatrix(Eigen::MatrixXd datMat,const std::function<double (arma::vec, arma::vec, arma::vec)> kernel, arma::vec parms){
  //Get the number of rows
  int rows=datMat.rows();
  //Typecasting
  arma::mat datMat2 = convertEigenToArma(datMat);
  //Initialize the matriz
  Eigen::MatrixXd matKernel = Eigen::MatrixXd::Zero(rows,rows);
  for(unsigned int c1=0;c1<rows;c1++){
    for(unsigned int c2=c1;c2<rows;c2++){
      //First column with variables
      arma::vec vec1 = datMat2.row(c1).t();
      //Second column with variables
      arma::vec vec2 = datMat2.row(c2).t();
      //Calculate the kernel value
      double val= kernel(vec1,vec2,parms);
      //Store the kernel value
      matKernel(c1,c2)=matKernel(c2,c1)=val;
    }
  }

  return(matKernel);
}
