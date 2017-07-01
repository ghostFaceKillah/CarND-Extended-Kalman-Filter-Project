#include "kalman_filter.h"
#include "tools.h"

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
  // predict the state

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::UpdateLinear(const Eigen::VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  Update(z, z_pred);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_pred = CartesianToPolar(x_);
  Update(z, z_pred);
}

void KalmanFilter::Update(const VectorXd &z, const VectorXd &z_pred) {
  // update the state by using Kalman Filter equations
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  VectorXd y = z - z_pred;                         // compute innovation
  MatrixXd S = R_ + H_ * P_ * H_.transpose();      // compute innovation covariance
  MatrixXd K = P_ * H_.transpose() * S.inverse();  // compute optimal Kalman gain
  x_ = x_ + K * y;                                 // update the state
  P_ = (I - K * H_) * P_;                          // update the state covariance
}
