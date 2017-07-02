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

/**
 * A pair of methods for converting between polar (radar)
 * and carthesian coordinates.
 * @param observation_vector
 * @return
 */
VectorXd PolarToCartesian(const VectorXd& observation_vector);
VectorXd CartesianToPolar(const VectorXd& state_vector);

#endif /* TOOLS_H_ */
