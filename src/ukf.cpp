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

  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;

  n_aug_ = n_x_ + 2;

  lambda_ = 3 - n_aug_;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.1;

  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for ( int i = 1; i < 2*n_aug_+1; ++i)
    weights_(i) = 0.5/(lambda_+n_aug_);

  NIS_radar_ = 0.0;

  NIS_laser_ = 0.0;
}

bool UKF::processMeasurement(Eigen::VectorXd &result, double &nis, const MeasurementRadar *measurement)
{
  std::cout << "processing radar measurement\n" << measurement->getCartesianMeasurement() << std::endl;
  if (!initialMeasurementProcessing(measurement))
    return false;
  UpdateRadar(measurement);
  result = x_;
  std::cout << "\tresult: " << result << std::endl;
  return true;
}

bool UKF::processMeasurement(Eigen::VectorXd &result, double &nis, const MeasurementLaser *measurement)
{
  std::cout << "processing laser measurement\n" << measurement->getCartesianMeasurement() << std::endl;
  if (!initialMeasurementProcessing(measurement))
    return false;
  UpdateLidar(measurement);
  result = x_;
  std::cout << "\tresult: " << result << std::endl;
  return true;
}

bool UKF::initialMeasurementProcessing(const Measurement *measurement)
{
  if (!measurement->valid())
    return false;

  if (!is_initialized_) {
    init(measurement);
    return false;
  }

  double delta_t = getDeltaTime(measurement);
  Prediction(delta_t);

  return true;
}

double UKF::getDeltaTime(const Measurement *measurement)
{
  double dt = (measurement->getTimestamp() - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = measurement->getTimestamp();
  return dt;
}

void UKF::init(const Measurement *measurement) {
  x_ = measurement->getCartesianMeasurement();
  time_us_ = measurement->getTimestamp();
  is_initialized_ = true;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //1) Create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  std::cout << "Xsig_aug:\n" << Xsig_aug << std::endl;

  //2) Predict sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    double px = Xsig_aug.col(i)(0);
    double py = Xsig_aug.col(i)(1);
    double v = Xsig_aug.col(i)(2);
    double psi = Xsig_aug.col(i)(3);
    double dpsi = Xsig_aug.col(i)(4);
    double nu_a = Xsig_aug.col(i)(5);
    double nu_d2psi = Xsig_aug.col(i)(6);

    VectorXd xkp1 = VectorXd(5);

    double p2 = v + delta_t*nu_a;
    double p3 = psi + dpsi * delta_t + 0.5*delta_t*delta_t*nu_d2psi;
    double p4 = dpsi + delta_t * nu_d2psi;
    if (fabs(dpsi) > 0.001 )
    {
      double p0 = px + (v/dpsi)*(sin(psi + dpsi*delta_t) - sin(psi)) + 0.5*delta_t*delta_t*cos(psi)*nu_a;
      double p1 = py + (v/dpsi)*(-cos(psi + dpsi*delta_t) + cos(psi)) + 0.5*delta_t*delta_t*sin(psi)*nu_a;
      xkp1 << p0, p1, p2, p3, p4;
    }
    else
    {
      double p0 = v*cos(psi)*delta_t + 0.5*delta_t*delta_t*cos(psi)*nu_a;
      double p1 = v*sin(psi)*delta_t + 0.5*delta_t*delta_t*sin(psi)*nu_a;
      xkp1 << p0, p1, p2, p3, p4;
    }
    Xsig_pred_.col(i) = xkp1;
  }

  std::cout << "XSIG_PRED:\n" << Xsig_pred_ << std::endl;

  //3) Predict state and covariance matrix
  x_ = VectorXd(n_x_);

  P_ = MatrixXd(n_x_, n_x_);

  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  std::cout << "Predicted x: " << x_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementLaser} lidar measurement
 */
void UKF::UpdateLidar(const MeasurementLaser *measurement) {
  int n_z = 2; //dimension of radar measurement space

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //TODO: Move R to member variable
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_*std_laspx_, 0,
      0, std_laspy_*std_laspy_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = measurement->getMeasurement() - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  /**
  TODO: calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementRadar} radar measurement
 */
void UKF::UpdateRadar(const MeasurementRadar *measurement) {
  int n_z = 3; //dimension of radar measurement space

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  std::cout << "ZSIG:\n" << Zsig << std::endl;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  std::cout << "z_pred: " << z_pred << std::endl;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //TODO: Move R to member variable
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = measurement->getMeasurement() - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  /**
  TODO: calculate the radar NIS.
  */
}
