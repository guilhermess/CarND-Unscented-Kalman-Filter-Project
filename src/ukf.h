#ifndef UKF_H
#define UKF_H

#include "MeasurementRadar.hpp"
#include "MeasurementLaser.hpp"
#include "Math.hpp"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <set>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class TestUKF;

class UKF {
public:

  UKF(bool use_laser, bool use_radar,
      double std_a, double std_yawdd,
      double std_laspx, double std_laspy,
      double std_radr, double std_radphi, double std_radrd);

  bool processMeasurement(VectorXd &result, double &nis, const MeasurementRadar *measurement);

  bool processMeasurement(VectorXd &result, double &nis, const MeasurementLaser *measurement);

  //Predict Methods
  MatrixXd getAugmentedSigmaPoints(const VectorXd &x, const MatrixXd &P);
  MatrixXd predictSigmaPoints(double delta_t, const MatrixXd &Xsig_aug);

  //Measurement Update Methods
  MatrixXd transformToLidarMeasurementSpace(const MatrixXd &Xsig_pred);
  MatrixXd transformToRadarMeasurementSpace(const MatrixXd &Xsig_pred);
  MatrixXd getCrossCorrelationMatrix(const VectorXd &x,
                                     const MatrixXd &Xsig_pred,
                                     const VectorXd &z_pred,
                                     const MatrixXd &Zsig,
                                     int n_z,
                                     std::set<int> z_angle_indexes,
                                     std::set<int> x_angle_indexes);
  MatrixXd getKalmanGain(const MatrixXd &S, const MatrixXd &Tc);


private:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise variance longitudinal acceleration
  double var_a_;

  ///* Process noise variance yaw acceleration
  double var_yawdd_;

  ///* Laser measurement noise variance position1
  double var_laspx_;

  ///* Laser measurement noise variance position2
  double var_laspy_;

  ///* Radar measurement noise variance radius
  double var_radr_;

  ///* Radar measurement noise variance angle in rad
  double var_radphi_;

  ///* Radar measurement noise variance radius change
  double var_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  MatrixXd R_radar_;

  MatrixXd R_lidar_;

  bool initialMeasurementProcessing(const Measurement *measurement);

  double getDeltaTime(const Measurement *measurement);

  void init(const Measurement *measurement);

  void Prediction(double delta_t);

  void UpdateLidar(double &nis, const MeasurementLaser *measurement);

  void UpdateRadar(double &nis, const MeasurementRadar *measurement);
};

#endif /* UKF_H */
