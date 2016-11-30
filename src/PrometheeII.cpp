#include <RcppEigen.h>
#include "Preference.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;


//Error Measure
void GaussianPreference(Eigen::MatrixXd &matDelta, double sigma);


// Create the Kernel matrix
// @param datVec  Column of the dataset
// @param int Type od preference function
// @return Preference Matrix
Eigen::MatrixXd matPrometheeII(Eigen::VectorXd datVec,int prefFunction){
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

  //Gaussian Preference
  if(prefFunction==0){
    GaussianPreference(matPromethee,sigma);
  }
  else if(prefFunction==1){
  }

  return(matPromethee);
}


// Create the Kernel matrix
// @param datVec  Column of the dataset
// @param int Type od preference function
// @return Preference Matrix
// [[Rcpp::export]]
Eigen::VectorXd PrometheeII(Eigen::MatrixXd datMat, Eigen::VectorXd vecWeights, Eigen::VectorXi prefFunction, bool normalize){
  //Get the number of rows
  int rows=datMat.rows();

  //Get the number of columns
  int cols=datMat.cols();

  //Initialize the matPromethee
  Eigen::MatrixXd matPromethee = Eigen::MatrixXd::Zero(rows, rows);

  //For each column
  for(int col=0;col<cols;col++){
    if(vecWeights(col)>0){
      //Create the delta matrix
      Eigen::MatrixXd matTemp  = matPrometheeII(datMat.col(col), prefFunction(col));
      //Multiply the weight
      matTemp = matTemp*vecWeights(col);
      //Accumulate the matrix
      matPromethee = matPromethee+matTemp;
    }
  }

  //Normalize the matrix
  matPromethee = matPromethee/vecWeights.sum();

  //Create the Flow vector
  Eigen::VectorXd phiPlus  = Eigen::VectorXd::Zero(rows);
  Eigen::VectorXd phiMinus = Eigen::VectorXd::Zero(rows);

  for(int row=0;row<rows;row++){
    phiPlus(row) = matPromethee.row(row).sum();
    phiMinus(row) = matPromethee.col(row).sum();
  }

    //Index
  Eigen::VectorXd phi = phiPlus-phiMinus;

  //Normalize
  double min = phi.minCoeff();
  double max = phi.maxCoeff();
  if(normalize==true) phi = (phi.array() - min)/(max-min);
  return (phi);
}

