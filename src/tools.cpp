#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  //accumulate squared residuals
  //need to compute vx and vy
  for(unsigned int i=0; i < ground_truth.size(); ++i){
    VectorXd tmp = estimations[i] - ground_truth[i];
    tmp = tmp.array()*tmp.array();
    rmse += tmp;
  }

  //calculate the mean
  rmse = rmse / ground_truth.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}