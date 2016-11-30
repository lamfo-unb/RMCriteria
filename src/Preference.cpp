#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

//Mean Square Error
void GaussianPreference(Eigen::MatrixXd &matDelta, double sigma){
  int rows = matDelta.rows();
  for(int i=0;i<rows;i++){
    for(int j=0;j<rows;j++){
      matDelta(i,j) = 1-std::exp(-1.0*std::pow(matDelta(i,j),2.0))
    }
  }
}
