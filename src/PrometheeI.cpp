#include <RcppEigen.h>
#include "Preference.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;


//Error Measure
void GaussianPreference(Eigen::MatrixXd &matDelta, double sigma);
//Usual Preference Function
void UsualPreference(Eigen::MatrixXd &matDelta);
//UShape Preference Function
void UShapePreference(Eigen::MatrixXd &matDelta, double q);
//VShape Preference Function
void VShapePreference(Eigen::MatrixXd &matDelta, double p);
//Level Preference Function
void LevelPreference(Eigen::MatrixXd &matDelta,double q, double p);
//Level Preference Function
void VShapeIndPreference(Eigen::MatrixXd &matDelta, double q, double p);


// Create the Kernel matrix
// @param datVec  Column of the dataset
// @param int Type od preference function
// @return Preference Matrix
Eigen::MatrixXd matPrometheeI(Eigen::VectorXd datVec,int prefFunction, Eigen::VectorXd parms){
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
    if(std::isnan(parms(0))){
      GaussianPreference(matPromethee,sigma);
    }
    else{
      sigma=parms(0);
      GaussianPreference(matPromethee,sigma);
    }

  }
  else if(prefFunction==1){
    UsualPreference(matPromethee);
  }
  else if(prefFunction==2){
    double q = parms(0);
    UShapePreference(matPromethee, q);
  }
  else if(prefFunction==3){
    double p=parms(0);
    VShapePreference(matPromethee, p);
  }
  else if(prefFunction==4){
    double q=parms(0);
    double p=parms(1);
    LevelPreference(matPromethee, q,  p);
  }
  else if(prefFunction==5){
    double q=parms(0);
    double p=parms(1);
    VShapeIndPreference(matPromethee, q,  p);
  }
  return(matPromethee);
}


// Create the Kernel matrix
// @param datVec  Column of the dataset
// @param int Type od preference function
// @return Preference Matrix
// [[Rcpp::export]]
List PrometheeI(Eigen::MatrixXd datMat, Eigen::VectorXd vecWeights, Eigen::VectorXi prefFunction, Eigen::MatrixXd parms, bool normalize){
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
      Eigen::MatrixXd matTemp  = matPrometheeI(datMat.col(col), prefFunction(col), parms.row(col));
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

  //Normalize
  double minPlus = phiPlus.minCoeff();
  double maxPlus = phiPlus.maxCoeff();
  double minMinus = phiMinus.minCoeff();
  double maxMinus = phiMinus.maxCoeff();
  if(normalize==true){
    phiPlus = (phiPlus.array() - minPlus)/(maxPlus-minPlus);
    phiMinus = (phiMinus.array() - minMinus)/(maxMinus-minMinus);
  }

  //'Store the results
  List resTemp = Rcpp::List::create(Rcpp::Named("Phi+") = phiPlus,
                                    Rcpp::Named("Phi-")  = phiMinus);

  return (resTemp);
}

