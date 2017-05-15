#ifndef MATH_HPP_
#define MATH_HPP_

#include <vector>
#include "Eigen/Dense"

namespace math {

  /**
   * Calculate RMSE.
   */
 Eigen::VectorXd
 CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /*
   * Normalize angle in radians between -PI and PI
   */
  float normalizeAngle(float angleRadian);

}

#endif