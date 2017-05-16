
#include "Math.hpp"
#include <iostream>

namespace math {

using namespace Eigen;

VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                       const std::vector<Eigen::VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if ((estimations.size() == 0) || (estimations.size() != ground_truth.size()))
  {
    std::cerr << "Could not calculate RMSE: estimations empty or not the same size as ground truth vector"
              << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  for (int i = 0; i < estimations.size(); ++i)
  {
    VectorXd estimated = estimations[i];
    VectorXd truth = ground_truth[i];
    VectorXd residual = estimated - truth;
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

double normalizeAngle(double angleRadian) {
  while ( angleRadian < -M_PI ) {
    angleRadian += 2 * M_PI;
  }
  while ( angleRadian > M_PI ) {
    angleRadian -= 2*M_PI;
  }
  return angleRadian;
}

Eigen::VectorXd stateVectorToPositionAndVelocity(const Eigen::VectorXd &state)
{
  VectorXd result{4};
  result << state(0), state(1),
      state(2) * cos(state(3)),
      state(2) * sin(state(3));
  return result;
}

void mean(VectorXd &result, const MatrixXd &matrix, const VectorXd &weights, int size)
{
  //predict state mean
  result.fill(0.0);
  for (int i = 0; i < size; i++)
    result = result + weights(i) * matrix.col(i);
}

void covariance(MatrixXd &result,
                const VectorXd &mean_value,
                const MatrixXd &matrix,
                const VectorXd &weights,
                int size,
                std::set<int> angle_indexes)
{
  //predicted state covariance matrix
  result.fill(0.0);
  for (int i = 0; i < size; i++) {
    VectorXd diff = matrix.col(i) - mean_value;

    for (auto angle_index : angle_indexes)
      diff(angle_index) = math::normalizeAngle(diff(angle_index));

    result = result + weights(i) * diff * diff.transpose() ;
  }
}

unsigned int countInRange(const std::vector<double> data, double lower_bound, double upper_bound)
{
  return std::count_if(data.begin(), data.end(),
                       [=](double val) { return (val > lower_bound && val < upper_bound); });
}


}
