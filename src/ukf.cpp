#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 7.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;
  time_us_ = 0;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  P_ << 0.0043,  -0.0013,   0.0030,  -0.0022,  -0.0020,
       -0.0013,   0.0077,   0.0011,   0.0071,   0.0060,
        0.0030,   0.0011,   0.0054,   0.0007,   0.0008,
       -0.0022,   0.0071,   0.0007,   0.0098,   0.0100,
       -0.0020,   0.0060,   0.0008,   0.0100,   0.0123;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // set weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i < 2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage &meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    x_ << 0, 0, 0, 0, 0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      float py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
      x_ << px, py, 0, 0, 0;
    }

    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }

  // Generated Sigma Points and augmentation
  sigma_points_matrix_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  VectorXd x_aug_ = VectorXd(n_aug_);

  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  MatrixXd sq_root_P_aug_ = P_aug_.llt().matrixL();

  sigma_points_matrix_aug_.col(0) = x_aug_;
  for (int i=0; i < n_aug_; ++i) {
    sigma_points_matrix_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_) * sq_root_P_aug_.col(i);
    sigma_points_matrix_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * sq_root_P_aug_.col(i);
  }

  // Predict
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);
  std::cout << meas_package.sensor_type_ << std::endl;

  // Update
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }

  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Predict x_(k+1) elements
  for (int i=0; i< 2*n_aug_+1; ++i) {
    //extract values for better readability
    double p_x = sigma_points_matrix_aug_(0,i);
    double p_y = sigma_points_matrix_aug_(1,i);
    double v = sigma_points_matrix_aug_(2,i);
    double yaw = sigma_points_matrix_aug_(3,i);
    double yawd = sigma_points_matrix_aug_(4,i);
    double nu_a = sigma_points_matrix_aug_(5,i);
    double nu_yawdd = sigma_points_matrix_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // Predicted Mean and Covariance
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage &meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  Zsig.fill(0.0);

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  std::cout << "Lidar" << std::endl;

  //Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  //transform sigma points into measurement space
  for(int i=0; i<2*n_aug_+1; i++) {
      float px = Xsig_pred_(0,i);
      float py = Xsig_pred_(1,i);

      Zsig(0,i) = px;
      Zsig(1,i) = py;
  }

  //calculate mean predicted measurement
  for(int i=0; i<2*n_aug_+1; i++){
      z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  for(int i=0; i<2*n_aug_+1; i++){
      VectorXd diff = VectorXd(n_z);
      diff = Zsig.col(i) - z_pred;

      S = S + weights_(i)*diff*diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;

  S = S + R;

  MatrixXd Tc = MatrixXd(n_x_, n_z); //5x2
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for(int i=0; i<2 * n_aug_ + 1; i++){
      VectorXd z_diff = Zsig.col(i) - z_pred;
      VectorXd x_diff = Xsig_pred_.col(i) - x_;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse(); //5x2 * 2*2

  //update state mean and covariance matrix
  x_ = x_ + K*(meas_package.raw_measurements_ - z_pred);
  P_ = P_ - K*S*K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage &meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // set measurement dimension, ro, phi, and rodot.
  std::cout << "Radar" << std::endl;

  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  Zsig.fill(0.0);

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  //transform sigma points into measurement space
  for(int i=0; i<2*n_aug_+1; i++) {
      float px = Xsig_pred_(0,i);
      float py = Xsig_pred_(1,i);
      float v = Xsig_pred_(2,i);
      float psi = Xsig_pred_(3,i);

      float ro = sqrt(px*px + py*py);
      float phi = atan2(py, px);
      float ro_dot = (px*cos(psi)*v + py*sin(psi)*v)/ro;

      Zsig(0,i) = ro;
      Zsig(1,i) = phi;
      Zsig(2,i) = ro_dot;
  }

  //calculate mean predicted measurement
  for(int i=0; i<2*n_aug_+1; i++){
      z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  for(int i=0; i<2*n_aug_+1; i++){
      VectorXd diff = Zsig.col(i) - z_pred;

      if(diff[1] > M_PI){
          diff[1] -= 2*M_PI;
      } else if (diff[1] < -M_PI){
          diff[1] += 2*M_PI;
      }

      S = S + weights_(i)*diff*diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;

  S = S + R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for(int i=0; i<2 * n_aug_ + 1; i++){
      VectorXd z_diff = Zsig.col(i) - z_pred;
      VectorXd x_diff = Xsig_pred_.col(i) - x_;

      if(z_diff(1) > M_PI){
          z_diff(1) = z_diff(1) - 2*M_PI;
      } else if (z_diff(1) < -M_PI){
          z_diff(1) = z_diff(1) + 2*M_PI;
      }

      if(x_diff(3) > M_PI){
          x_diff(3) = x_diff(3) - 2*M_PI;
      } else if (z_diff(1) < -M_PI){
          x_diff(3) = x_diff(3) + 2*M_PI;
      }

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //update state mean and covariance matrix
  x_ = x_ + K*(meas_package.raw_measurements_ - z_pred);
  P_ = P_ - K*S*K.transpose();
}
