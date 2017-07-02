#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size() || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for (unsigned int i = 0; i < estimations.size(); ++i) {

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array() * residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;

}

VectorXd PolarToCartesian(const VectorXd& observation_vector) {
    VectorXd cartesian(4);

    double r_measured = observation_vector(0);
    double phi_measured = observation_vector(1);
    double r_dot_measured = observation_vector(2);

    double x = r_measured * cos(phi_measured);
    double y = r_measured * sin(phi_measured);

    double x_dot = r_dot_measured * cos(phi_measured);
    double y_dot = r_dot_measured * sin(phi_measured);

    cartesian << x, y, x_dot, y_dot;

    return cartesian;
};

VectorXd CartesianToPolar(const VectorXd& state_vector) {
    VectorXd polar(3);

    double x = state_vector[0];
    double y = state_vector[1];
    double x_dot = state_vector[2];
    double y_dot = state_vector[3];

    double rho = sqrt(x * x + y * y);

    if (fabs(rho) < 1e-5) {
        std::cout << "Division by zero in CartesianToPolarConverter "
                     "- returning all zeros state" << std::endl;
        return VectorXd(3);
    }

    double phi = atan2(y, x);
    double rho_dot = (x * x_dot + y * y_dot) / rho;

    polar << rho, phi, rho_dot;

    return polar;
};

MatrixXd CalculateRadarJacobian(const VectorXd& x_state) {
    // Calculate Jacobian - a linear approximation of radar measurment function

    MatrixXd H(3, 4);

    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    double norm_squared = px * px + py * py;
    if (fabs(norm_squared) < 1e-4) {
        cout << "Division by zero in CalculateRadarJacobian "
                "- returning zero H matrix" << endl;
        return H;
    }
    double norm = sqrt(norm_squared);
    double norm_3 = norm_squared * norm;

    H << px / norm,                         py / norm,                          0,         0,
         -py / norm_squared,                px / norm_squared,                  0,         0,
         py * (vx * py - vy * px) / norm_3, px * ( vy * px - vx * py) / norm_3, px / norm, py / norm;

    return H;
}
