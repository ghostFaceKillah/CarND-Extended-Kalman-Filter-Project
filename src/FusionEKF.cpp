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

  // what is this guy used for ??
  Hj_ = MatrixXd(3, 4);

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



  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  // estimated state mean will be initialized elsewhere
  // estimated state covariance
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0,   0,
             0, 1, 0,   0,
             0, 0, 1e3, 0,
             0, 0, 0,   1e3;

  // state transition model - parts wil lmodified
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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.

      double r_measured = measurement_pack.raw_measurements_[0];
      double phi_measured = measurement_pack.raw_measurements_[1];
      double x = r_measured * cos(phi_measured);
      double y = r_measured * sin(phi_measured);

      ekf_.x_ = VectorXd(4);
      ekf_.x_ << x, y, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      double x = measurement_pack.raw_measurements_[0];
      double y = measurement_pack.raw_measurements_[1];

      ekf_.x_ << x, y, 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     [X] Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     [X] Update the process noise Q covariance matrix.
   */

  // where does the funny Q_ formulation come from ??

  // compute the time elapsed between the current and previous measurments
  // it is possible that these time measurements in seconds (insteaad of Unix timestamps)
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

  //set the process covariance matrix Q
  // the form is nicely explained in 5.9 on the Udacity if you want to study it
  // in more detail
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0,               dt_3/2*noise_ax, 0,
              0,               dt_4/4*noise_ay, 0,               dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0,               dt_2*noise_ax,   0,
              0,               dt_3/2*noise_ay, 0,               dt_2*noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */


  // aha - we need to set up good matrix H for given measurment type
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    // TODO: is it correct to use ekf_.x_ here?
    ekf_.H_ = CalculateRadarJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

MatrixXd CalculateRadarJacobian(const VectorXd& x_state) {
   MatrixXd H(3, 4);

  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  double norm_squared = px * px + py * py;
  double norm = sqrt(norm_squared);
  double norm_3_2 = pow(norm, 1.5);

  H << px / norm,                            py / norm,                            0,         0,
        -py / norm_squared,                  px / norm_squared,                    0,         0,
        py * (vx * py - vy * px) / norm_3_2, px * ( vy * px - vx * py) / norm_3_2, px / norm, py / norm;

  return H;
}

