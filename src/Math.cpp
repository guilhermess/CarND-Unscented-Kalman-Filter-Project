
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

  float normalizeAngle(float angleRadian) {
    while ( angleRadian < -M_PI )
      angleRadian += 2*M_PI;
    while ( angleRadian > M_PI )
      angleRadian -= 2*M_PI;
    return angleRadian;
  }

}
