#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;


//Usual Preference Function
void UsualPreference(Eigen::MatrixXd &matDelta){
  int rows = matDelta.rows();
  for(int i=0;i<rows;i++){
    for(int j=0;j<rows;j++){
      if(matDelta(i,j)<=0){
        matDelta(i,j) = 0.0;
      }
      else{
        matDelta(i,j) = 1.0;
      }

    }
  }
}


//UShape Preference Function
void UShapePreference(Eigen::MatrixXd &matDelta, double q){
  int rows = matDelta.rows();
  for(int i=0;i<rows;i++){
    for(int j=0;j<rows;j++){
      if(matDelta(i,j)<=q){
        matDelta(i,j) = 0.0;
      }
      else{
        matDelta(i,j) = 1.0;
      }

    }
  }
}

//VShape Preference Function
void VShapePreference(Eigen::MatrixXd &matDelta, double p){
  int rows = matDelta.rows();
  for(int i=0;i<rows;i++){
    for(int j=0;j<rows;j++){
      if(matDelta(i,j)<=0){
        matDelta(i,j) = 0.0;
      }
      else if(matDelta(i,j)<=p){
        matDelta(i,j) = matDelta(i,j)/p;
      }
      else{
        matDelta(i,j) = 1.0;
      }
    }
  }
}

//Level Preference Function
void LevelPreference(Eigen::MatrixXd &matDelta,double q, double p){
  int rows = matDelta.rows();
  for(int i=0;i<rows;i++){
    for(int j=0;j<rows;j++){
      if(matDelta(i,j)<=q){
        matDelta(i,j) = 0.0;
      }
      else if(matDelta(i,j)<=p){
        matDelta(i,j) = 0.5;
      }
      else{
        matDelta(i,j) = 1.0;
      }
    }
  }
}

//Level Preference Function
void VShapeIndPreference(Eigen::MatrixXd &matDelta, double q, double p){
  int rows = matDelta.rows();
  for(int i=0;i<rows;i++){
    for(int j=0;j<rows;j++){
      if(matDelta(i,j)<=q){
        matDelta(i,j) = 0.0;
      }
      else if(matDelta(i,j)<=p){
        matDelta(i,j) = (matDelta(i,j)-q)/(p-q);
      }
      else{
        matDelta(i,j) = 1.0;
      }
    }
  }
}

//Gaussian Preference Function
void GaussianPreference(Eigen::MatrixXd &matDelta, double sigma){
  int rows = matDelta.rows();
  for(int i=0;i<rows;i++){
    for(int j=0;j<rows;j++){
      if(matDelta(i,j)<=0){
        matDelta(i,j) = 0.0;
      }
      else{
        matDelta(i,j) = 1-std::exp((-1.0*std::pow(matDelta(i,j),2.0))/(2.0*sigma*sigma));
      }

    }
  }
}
