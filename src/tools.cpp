#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /** Calculate root mean squared error. 
   * Straight from the lessons.
   */
  assert(estimations.size() > 0); 
  assert(estimations.size() == ground_truth.size()); 

  VectorXd rmse = VectorXd(4);
  rmse <<  0, 0, 0, 0;
  
  // collect squared error, compute the mean
  for(int i=0; i < estimations.size(); ++i){
    
    VectorXd rmse_working = estimations[i] - ground_truth[i]; 
    rmse_working = rmse_working.array() * rmse_working.array();
    rmse += rmse_working;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();

  return rmse; 
}

float Tools::WrapAnglePi(const float angle)  {
  // scale module 2pi to  wrap the angle between [-pi, pi]
  float angle_result = angle;
  if ( angle > 0 ) { 
    angle_result = fmod(angle + M_PI, 2*M_PI) - M_PI;
  } else {
    angle_result = fmod(angle - M_PI, 2*M_PI) + M_PI;
  }
  
  return angle_result;
}
