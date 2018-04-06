// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef PREFERENCE_H
#define PREFERENCE_H

//Gaussian Preference Function
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



// This is the end of the header guard
#endif


