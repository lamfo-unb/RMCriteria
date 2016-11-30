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
  double sum  = 0;
  double cont = 0;
  //Initialize the computation
  for(int i=0;i<rows;i++){
    for(int j=i+1;j<rows-1;j++){
      //Compute the difference
      double delta = datVec(i)-datVec(j);
      if(delta>0){
        matPromethee(i,j)=delta;
      }
      else{
        matPromethee(j,i)=(1.0)*delta;
      }
      sum2=sum2+std::pow(matPromethee(j,i),2);
      sum=sum+delta;
      cont=cont+1;
    }
  }

  //Calculate the sigma
  double sigma = (sum2/cont)-(std::pow(sum/cont,2));

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
Eigen::VectorXd PrometheeII(Eigen::MatrixXd datMat, Eigen::VectorXd vecWeights, Eigen::VectorXi prefFunction){
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

  //Create the Flow matrix
  Eigen::MatrixXd matFlow = Eigen::MatrixXd::Zero(cols, cols);
  for(int col=0;col<cols;col++){
    matFlow.row(col) = matPromethee.row(col).sum();
  }

}

