#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

/**
* A helper method to calculate RMSE.
*/
VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

/**
* A helper method to calculate Jacobians.
*/
MatrixXd CalculateRadarJacobian(const VectorXd& x_state);

VectorXd PolarToCartesian(const VectorXd& observation_vector);
VectorXd CartesianToPolar(const VectorXd& state_vector);

#endif /* TOOLS_H_ */
