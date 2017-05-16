#ifndef MATH_HPP_
#define MATH_HPP_

#include <vector>
#include <set>
#include "Eigen/Dense"

namespace math {

using Eigen::VectorXd;
using Eigen::MatrixXd;

VectorXd CalculateRMSE(const std::vector<VectorXd> &estimations,
                       const std::vector<VectorXd> &ground_truth);

double normalizeAngle(double angleRadian);

Eigen::VectorXd stateVectorToPositionAndVelocity(const VectorXd &state);

void mean(VectorXd &result, const MatrixXd &matrix, const VectorXd &weights, int size);

void covariance(MatrixXd &result,
                const VectorXd &mean_value,
                const MatrixXd &matrix,
                const VectorXd &weights,
                int size,
                std::set<int> angle_indexes);

unsigned int countInRange(const std::vector<double> data, double lower_bound, double upper_bound);

}

#endif