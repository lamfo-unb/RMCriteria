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


// [[Rcpp::export]]
Eigen::MatrixXd matPrometheeIII(Eigen::VectorXd datVec,int prefFunction, Eigen::VectorXd parms){
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
  double sigma2 = (sum2/(2*cont));

  //Gaussian Preference
  if(prefFunction==0){
    if(std::isnan(parms(0))){
      GaussianPreference(matPromethee,std::sqrt(sigma2));
    }
    else{
      double sigma=parms(0);
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


//' Calculates PROMETHEE III method.
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
//' @param alphaVector A numerical vector to indicate the size of the interval
//' for each alternative in Promethee III ranking.
//' @param parms a numerical matrix with parameters associated to the Preference
//'  Function. They're defined as a matrix of n columns and m rows. The maximum
//'  number of parameters is 3 and m is the number of criterias. The parameters
//'  are:
//'   \itemize{
//'   \item{Indifference Threshold (\code{q})}
//'   \item{Preference Threshold (\code{p})}
//'   \item{Gaussian Threshold (\code{s})}
//'   }
//' @return Preference Matrix
//' @export
// [[Rcpp::export]]


List PrometheeIII(Eigen::MatrixXd datMat, Eigen::VectorXd vecWeights, Eigen::VectorXi prefFunction,Eigen::VectorXi alphaVector, Eigen::MatrixXd parms){
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
      Eigen::MatrixXd matTemp  = matPrometheeIII(datMat.col(col), prefFunction(col), parms.row(col));
      //Multiply the weight
      matTemp = matTemp*vecWeights(col);
      //Accumulate the matrix
      matPromethee = matPromethee+matTemp;
    }
  }

  //Normalize the matrix
  matPromethee = matPromethee/vecWeights.sum();

  //Create the Flow vector
  Eigen::VectorXd limInf  = Eigen::VectorXd::Zero(rows);
  Eigen::VectorXd limSup = Eigen::VectorXd::Zero(rows);

  for(int row=0;row<rows;row++){
    Eigen::VectorXd rowMat =  matPromethee.row(row);
    Eigen::VectorXd colMat = matPromethee.col(row);
    Eigen::VectorXd diff = rowMat-colMat;
    int nElements = matPromethee.cols();

    //Compute mean and variance
    double phiBar =diff.mean();
    double phiVar = diff.squaredNorm()/nElements - (phiBar*phiBar);

    limInf(row) = phiBar-alphaVector(row)*std::sqrt(phiVar);
    limSup(row) = phiBar+alphaVector(row)*std::sqrt(phiVar);
  }

  //'Store the results
  List resTemp = Rcpp::List::create(Rcpp::Named("limInf") = limInf,
                                    Rcpp::Named("limSup")  = limSup);

  return (resTemp);
}

