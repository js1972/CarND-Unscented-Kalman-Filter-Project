#include <iostream>
#include "ukf.h"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "IncompatibleTypes"
using namespace std;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5; //30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55; //30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.0175; //0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.1; //0.3;

  previous_timestamp_ = 0;

  n_x_ = 5;
  n_aug_ = 7;

  lambda_ = 3 - n_aug_;
  n_sigma_ = 2 * n_aug_ + 1;

  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

  // calculate weights
  weights_ = VectorXd(n_sigma_);
  weights_(0)= lambda_/(lambda_ + n_aug_);
  for (int i=1; i<n_sigma_; i++) {
    weights_(i) = 1/(2*(lambda_ + n_aug_));
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    previous_timestamp_  = meas_package.timestamp_;
    cout << "Measurement package raw measurements: " << endl << meas_package.raw_measurements_ << endl;

    // initialize the state vector
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR and use_radar_) {
      float ro = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      x_ << ro * cos(phi), ro * sin(phi), 0, 0, 0;
      cout << "\nRADAR" << endl;
      DebugPrintState();

      is_initialized_ = true;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER and use_laser_) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      cout << "\nLASER" << endl;
      DebugPrintState();

      is_initialized_ = true;

    } else {
      std::cout << "Invalid sensor type" << std::endl;
    }
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // Compute the time elapsed between the current and previous measurements, in seconds
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;

  // skip prediction, if measurement has (about) the same timestamp as previous
  // timestamps are microseconds
  if (previous_timestamp_ < meas_package.timestamp_ - 100) {
    Prediction(dt);
    // update prediction timestamp only if we did prediction
    previous_timestamp_ = meas_package.timestamp_;
  } else {
    cout << "measurement skipped : dt < 100msec" << endl;
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR and use_radar_) {
    cout << "New RADAR measurement" << endl;
    UpdateRadar(meas_package);

  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER and use_laser_) {
    cout << "New LASER measurement" << endl;
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  cout << "PREDICTION()" << endl;

  // Create the sigma point matrix
  MatrixXd Xsig_aug(n_aug_, n_sigma_);
  //MatrixXd Xsig_aug_(n_aug_, n_sigma_);

  GenerateAugmentedSigmaPoints(x_, P_, Xsig_aug);
  //GenerateAugmentedSigmaPoints();
  SigmaPointPrediction(delta_t, Xsig_aug, Xsig_pred_);
  PredictMeanAndCovariance(Xsig_pred_, x_, P_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  cout << "UPDATE LIDAR" << endl;

  int n_y = meas_package.raw_measurements_.rows();

  VectorXd y_pred(n_y);
  MatrixXd S(n_y, n_y);
  MatrixXd Ysig(n_y, n_sigma_);

  PredictLidarMeasurement(Xsig_pred_, Ysig, y_pred, S);
  UpdateState(meas_package.raw_measurements_, y_pred, S, Xsig_pred_, Ysig, x_, P_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = meas_package.raw_measurements_.rows();

  VectorXd z_pred(n_z);
  MatrixXd S(n_z, n_z);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig(n_z, n_sigma_);

  PredictRadarMeasurement(Xsig_pred_, Zsig, z_pred, S);
  UpdateState(meas_package.raw_measurements_, z_pred, S, Xsig_pred_, Zsig, x_, P_);
}

void UKF::GenerateAugmentedSigmaPoints(const VectorXd &x, const MatrixXd &P, MatrixXd &Xsig_aug) {

  //set state dimension
  int n_x = x.rows();

  //create augmented mean vector
  VectorXd x_aug(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug(n_aug_, n_aug_);

  //create augmented mean state - adding zero for the extra noise values on the end
  x_aug << x, 0, 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x, n_x) = P;
  P_aug(n_aug_-2,n_aug_-2) = std_a_ * std_a_;
  P_aug(n_aug_-1,n_aug_-1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
}

void UKF::SigmaPointPrediction(double delta_t, const MatrixXd &Xsig_aug, MatrixXd &Xsig_pred) {

  //predict sigma points
  for (int i = 0; i < n_sigma_; i++) {

    //extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance(const MatrixXd &Xsig_pred, VectorXd &x, MatrixXd &P) {
  cout << "PREDICT MEAN AND COVARIANCE" << endl;

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred.col(i);
  }

  DebugPrintState();

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::PredictLidarMeasurement(const MatrixXd &Xsig_pred, MatrixXd &Ysig, VectorXd &y_pred, MatrixXd &S) {

  int n_y = y_pred.rows();

  // transform sigma points into measurement space
  for (int i = 0; i < n_sigma_; i++) {
    Ysig(0, i) = Xsig_pred(0, i);
    Ysig(1, i) = Xsig_pred(1, i);
  }

  // mean predicted measurement
  y_pred.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    y_pred = y_pred + weights_(i) * Ysig.col(i);
  }

  // measurement covariance matrix S
  S.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    // residual
    VectorXd y_diff = Ysig.col(i) - y_pred;
    S = S + weights_(i) * y_diff * y_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_y, n_y);
  R << std_laspx_*std_laspx_, 0,
       0,                     std_laspy_*std_laspy_;

  S = S + R;
}

void UKF::PredictRadarMeasurement(const MatrixXd &Xsig_pred, MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S) {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = z_pred.rows();

  // transform sigma points into measurement space

  for (int i = 0; i < n_sigma_; i++) {

    // extract values for better readability
    double px  = Xsig_pred(0,i);
    double py  = Xsig_pred(1,i);
    double v   = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(px*px + py*py);              // r
    if (fabs(Zsig(0, i)) > 0.001) {
      Zsig(1, i) = atan2(py,px);                   // phi
      Zsig(2, i) = (px*v1 + py*v2 ) / Zsig(0, i);  // r_dot
    } else {
      Zsig(1, i) = 0.0;  // phi
      Zsig(2, i) = 0.0;  // r_dot
    }
  }

  // mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix S
  S.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix to the measurement covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R <<    std_radr_*std_radr_, 0,                       0,
          0,                   std_radphi_*std_radphi_, 0,
          0,                   0,                       std_radrd_*std_radrd_;

  S = S + R;
}

void UKF::UpdateState(const VectorXd &z, const VectorXd &z_pred, const MatrixXd &S,
                      const MatrixXd &Xsig_pred, const MatrixXd &Zsig, VectorXd &x, MatrixXd &P) {
  cout << "UPDATE STATE" << endl;

  //set measurement dimension, radar can measure r, phi, and r_dot and lidar px and py
  int n_z = z_pred.rows();
  int n_x = x.rows();

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;
  DebugPrintState();

  P = P - K*S*K.transpose();
}

void UKF::DebugPrintState() {
  cout << "Current State Vector:" << endl;
  cout << "p_x\t\t\t" << x_(0) << endl;
  cout << "p_y\t\t\t" << x_(1) << endl;
  cout << "v\t\t\t" << x_(2) << endl;
  cout << "psi\t\t\t" << x_(3) << endl;
  cout << "psi_dot\t\t" << x_(4) << endl << endl;
}

#pragma clang diagnostic pop