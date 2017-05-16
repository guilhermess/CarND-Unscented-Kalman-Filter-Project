
#include "ukf.h"
#include "Math.hpp"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::set;

/**
 * Unscented Kalman Filter Constructor.
 * @param use_laser use laser measurements.
 * @param use_radar use radar measurements.
 * @param std_a acceleration standard deviation.
 * @param std_yawdd yaw standard deviation.
 * @param std_laspx laser px standard deviation.
 * @param std_laspy laser py standard deviation.
 * @param std_radr radar r standard deviation.
 * @param std_radphi radar phi standard deviation.
 * @param std_radrd radar r' standard deviation.
 */
UKF::UKF(bool use_laser, bool use_radar,
         double std_a, double std_yawdd,
         double std_laspx, double std_laspy,
         double std_radr, double std_radphi, double std_radrd) :
    is_initialized_{false},
    use_laser_{use_laser}, use_radar_{use_radar},
    n_x_{5}, n_aug_{n_x_+2}, lambda_{static_cast<double>(3 - n_aug_)},
    x_{VectorXd(n_x_)},
    P_{MatrixXd::Identity(5,5)},
    var_a_{std_a*std_a},
    var_yawdd_{std_yawdd*std_yawdd},
    var_laspx_{std_laspx*std_laspx},
    var_laspy_{std_laspy*std_laspy},
    var_radr_{std_radr*std_radr},
    var_radphi_{std_radphi*std_radphi},
    var_radrd_{std_radrd*std_radrd},
    weights_{VectorXd(2*n_aug_+1)},
    R_radar_{MatrixXd(3,3)}, R_lidar_{MatrixXd(2,2)}
{
  weights_(0) = lambda_/(lambda_ + n_aug_);
  weights_(1) = 0.5/(lambda_ + n_aug_);
  for ( int i = 2; i < 2*n_aug_+1; ++i)
    weights_(i) = weights_(1);

  R_lidar_ << var_laspx_, 0,
      0, var_laspy_;

  R_radar_ << var_radr_, 0, 0,
      0, var_radphi_, 0,
      0, 0, var_radrd_;
}

/**
 * Process Radar Measurement
 * @param result state vector result.
 * @param nis NIS reuslt.
 * @param measurement Radar Measurement
 * @return true if processing is successful, otherwise returns false.
 */
bool UKF::processMeasurement(Eigen::VectorXd &result, double &nis, const MeasurementRadar *measurement)
{
  if ( !use_radar_ )
    return false;

  if (!initialMeasurementProcessing(measurement))
    return false;
  UpdateRadar(nis, measurement);
  result = x_;
  return true;
}

/**
 * Process Lidar Measurement
 * @param result state vector result.
 * @param nis NIS result.
 * @param measurement Lidar measurement.
 * @return true if processing is successful, otherwise returns false.
 */
bool UKF::processMeasurement(Eigen::VectorXd &result, double &nis, const MeasurementLaser *measurement)
{
  if ( !use_laser_ )
    return false;
  if (!initialMeasurementProcessing(measurement))
    return false;
  UpdateLidar(nis, measurement);
  result = x_;
  return true;
}


/**
 * Initial measurement processing, common to both Radar and Lidar measurements.
 * @param measurement Radar or Lidar measurement
 * @return true if initial processing is successful, otherwise returns false.
 */
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

/**
 * Get the delta time between last measurement and current measurement. It also updates the
 * current time stamp variable time_us_;
 * @param measurement new measurement.
 * @return delta time (in seconds) from the last measurement to current measurement.
 */
double UKF::getDeltaTime(const Measurement *measurement)
{
  double dt = (measurement->getTimestamp() - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = measurement->getTimestamp();
  return dt;
}

/**
 * Initialize Unscented Kalman Filter with a first valid measurement.
 * @param measurement first measurement.
 */
void UKF::init(const Measurement *measurement) {
  x_ = measurement->getCartesianMeasurement();
  time_us_ = measurement->getTimestamp();
  is_initialized_ = true;
}

/**
 * Prediction step of Unscented Kalman Filter.
 * @param delta_t time step between last and current measurement.
 */
void UKF::Prediction(double delta_t) {
  auto Xsig_aug = getAugmentedSigmaPoints(x_, P_);
  Xsig_pred_ = predictSigmaPoints(delta_t, Xsig_aug);
  math::mean(x_, Xsig_pred_, weights_, 2 * n_aug_ + 1);
  math::covariance(P_, x_, Xsig_pred_, weights_, 2 * n_aug_ + 1, {3});
}

/**
 * Get augmented sigma points based on state vector x and state covariance matrix P
 * @param x current state vector.
 * @param P current covariance matrix
 * @return augmented sigma points.
 */
MatrixXd UKF::getAugmentedSigmaPoints(const VectorXd &x, const MatrixXd &P)
{
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P;
  P_aug(5,5) = var_a_;
  P_aug(6,6) = var_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //1) Create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  return Xsig_aug;
}


/**
 * Predict Sigma Points based on augmented sigma points and time step delta_t.
 * @param delta_t time step between last and current measurement.
 * @param Xsig_aug augmented sigma points.
 * @return predicted sigma points.
 */
MatrixXd UKF::predictSigmaPoints(double delta_t, const MatrixXd &Xsig_aug)
{
  auto Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    double px = Xsig_aug.col(i)(0);
    double py = Xsig_aug.col(i)(1);
    double v = Xsig_aug.col(i)(2);
    double yaw = Xsig_aug.col(i)(3);
    double yawd = Xsig_aug.col(i)(4);
    double nu_a = Xsig_aug.col(i)(5);
    double nu_yawdd = Xsig_aug.col(i)(6);

    VectorXd xkp1 = VectorXd(5);

    double p2 = v + delta_t*nu_a;
    double p3 = yaw + yawd * delta_t + 0.5*delta_t*delta_t*nu_yawdd;
    double p4 = yawd + delta_t * nu_yawdd;
    if (fabs(yawd) > 0.001 )
    {
      double p0 = px + (v/yawd)*(sin(yaw + yawd*delta_t) - sin(yaw)) + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
      double p1 = py + (v/yawd)*(-cos(yaw + yawd*delta_t) + cos(yaw)) + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
      xkp1 << p0, p1, p2, p3, p4;
    }
    else
    {
      double p0 = px + v*delta_t*cos(yaw);
      double p1 = py + v*delta_t*sin(yaw);
      xkp1 << p0, p1, p2, p3, p4;
    }
    Xsig_pred.col(i) = xkp1;
  }
  return Xsig_pred;
}



/**
 * Update state based on Radar Measurement.
 * @param nis NIS result variable.
 * @param measurement radar measurement.
 */
void UKF::UpdateRadar(double &nis, const MeasurementRadar *measurement) {
  int n_z = 3;

  auto Zsig = transformToRadarMeasurementSpace(Xsig_pred_);
  VectorXd z_pred(n_z);
  math::mean(z_pred, Zsig, weights_, 2 * n_aug_ + 1);

  MatrixXd S(n_z, n_z);
  math::covariance(S, z_pred, Zsig, weights_, 2 * n_aug_ + 1, {1});
  S = S + R_radar_;
  auto Tc = getCrossCorrelationMatrix(x_, Xsig_pred_, z_pred, Zsig, n_z, {1}, {3});

  //Kalman gain K;
  auto K = getKalmanGain(S, Tc);

  //residual
  VectorXd z_diff = measurement->getMeasurement() - z_pred;

  //angle normalization
  z_diff(1) = math::normalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  nis = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Update state based on Lidar Measurement.
 * @param nis NIS result variable.
 * @param measurement lidar measurement.
 */
void UKF::UpdateLidar(double &nis, const MeasurementLaser *measurement) {
  int n_z = 2;
  auto Zsig = transformToLidarMeasurementSpace(Xsig_pred_);

  VectorXd z_pred(n_z);
  math::mean(z_pred, Zsig, weights_, 2 * n_aug_ + 1);

  MatrixXd S(n_z, n_z);
  math::covariance(S, z_pred, Zsig, weights_, 2 * n_aug_ + 1, {});
  S = S + R_lidar_;
  auto Tc = getCrossCorrelationMatrix(x_, Xsig_pred_, z_pred, Zsig, n_z, {}, {});

  //Kalman gain K;
  auto K = getKalmanGain(S, Tc);

  //residual
  VectorXd z_diff = measurement->getMeasurement() - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  nis = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Get Kalman gain.
 * @param S measurement covariance matrix.
 * @param Tc Cross correlation matrix.
 * @return Kalman Gain matrix.
 */
MatrixXd UKF::getKalmanGain(const MatrixXd &S, const MatrixXd &Tc)
{
  MatrixXd K = Tc * S.inverse();
  return K;
}

/**
 * Get cross correlation matrix.
 * @param x process model mean.
 * @param Xsig_pred process model predicted sigma points.
 * @param z_pred measurement mean.
 * @param Zsig measurement space sigma points.
 * @param n_z measurement space state dimension.
 * @param z_angle_indexes dimensions in measurement space that represent angles.
 * @param x_angle_indexes dimensions in the process space that represent angles.
 * @return cross correlation matrix.
 */
MatrixXd UKF::getCrossCorrelationMatrix(const VectorXd &x,
                                        const MatrixXd &Xsig_pred,
                                        const VectorXd &z_pred,
                                        const MatrixXd &Zsig,
                                        int n_z,
                                        set<int> z_angle_indexes,
                                        set<int> x_angle_indexes)
{
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    for (auto angle_index : z_angle_indexes)
      z_diff(angle_index) = math::normalizeAngle(z_diff(angle_index));

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    for (auto angle_index : x_angle_indexes)
      x_diff(angle_index) = math::normalizeAngle(x_diff(angle_index));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  return Tc;
}

/**
 * Transform sigma points from process model space to lidar measurement space.
 * @param Xsig_pred sigma points in process model space.
 * @return sigma points in measurement space.
 */
MatrixXd UKF::transformToLidarMeasurementSpace(const MatrixXd &Xsig_pred)
{
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }
  return Zsig;
}

/**
 * Transform sigma points from process model space to radar measurement space.
 * @param Xsig_pred sigma points in process model space.
 * @return sigma points in measurement space.
 */
MatrixXd UKF::transformToRadarMeasurementSpace(const MatrixXd &Xsig_pred)
{
  int n_z = 3; //dimension of radar measurement space

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
  return Zsig;
}
