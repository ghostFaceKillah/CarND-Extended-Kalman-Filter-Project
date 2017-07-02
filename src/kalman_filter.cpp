#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::UpdateLinear(const Eigen::VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  Update(z - z_pred);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_pred = CartesianToPolar(x_);
  VectorXd innovation = z - z_pred;

  // bring the yaw innovation back into reasonable range
  while (innovation(1) >= M_PI_2) {
    innovation(1) -= M_PI;
  }
  while (innovation(1) <= -M_PI_2) {
    innovation(1) += M_PI;
  }

  Update(innovation);
}

void KalmanFilter::Update(const VectorXd &innovation) {
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  MatrixXd S = R_ + (H_ * P_ * H_.transpose());    // compute innovation covariance
  MatrixXd K = P_ * H_.transpose() * S.inverse();  // compute optimal Kalman gain
  x_ = x_ + K * innovation;                        // update the state
  P_ = (I - K * H_) * P_;                          // update the state covariance
}
