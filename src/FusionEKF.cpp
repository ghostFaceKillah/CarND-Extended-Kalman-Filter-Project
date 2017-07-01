#include "FusionEKF.h"
#include <math.h>
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>


/**
 * Data structure note copied from the course website
   For a row containing radar data, the columns are:
    * sensor_type,
    * rho_measured,
    * phi_measured,
    * rhodot_measured,
    * timestamp,
    * x_groundtruth,
    * y_groundtruth,
    * vx_groundtruth,
    * vy_groundtruth,
    * yaw_groundtruth,
    * yawrate_groundtruth.

    For a row containing lidar data, the columns are:
    *  sensor_type,
    *  x_measured,
    *  y_measured,
    *  timestamp,
    *  x_groundtruth,
    *  y_groundtruth,
    *  vx_groundtruth,
    *  vy_groundtruth,
    *  yaw_groundtruth,
    *  yawrate_groundtruth.

  Whereas radar has three measurements (rho, phi, rhodot),
  lidar has two measurements (x, y).
 */

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);


  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  // observation model
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // More intelligent value will be initialized on first measurmenet
  ekf_.x_ = VectorXd(4);

  // estimated state covariance
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0,   0,
             0, 1, 0,   0,
             0, 0, 1e3, 0,
             0, 0, 0,   1e3;

  // state transition model - parts will be modified based on time elapsed
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // Initialize the state ekf_.x_ with the first measurement.
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      ekf_.x_ = PolarToCartesian(measurement_pack.raw_measurements_)
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      double x = measurement_pack.raw_measurements_[0];
      double y = measurement_pack.raw_measurements_[1];
      ekf_.x_ << x, y, 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // Update the state transition matrix F according to the new elapsed time
  // but first convert elapsed time from Unix timestamp to seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1e6;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  double noise_ax = 9;
  double noise_ay = 9;

  // Update the process covariance matrix Q - for details see Udacity 5.9
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0,               dt_3/2*noise_ax, 0,
              0,               dt_4/4*noise_ay, 0,               dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0,               dt_2*noise_ax,   0,
              0,               dt_3/2*noise_ay, 0,               dt_2*noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = CalculateRadarJacobian(ekf_.x_);    // observation matrix
    ekf_.R_ = R_radar_;                           // observation covariance
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;  // observation matrix
    ekf_.R_ = R_laser_;  // observation covariance
    ekf_.UpdateLinear(measurement_pack.raw_measurements_);
  }
}

