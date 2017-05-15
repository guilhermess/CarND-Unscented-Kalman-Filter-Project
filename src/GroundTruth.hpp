#ifndef GROUND_TRUTH_H_
#define GROUND_TRUTH_H_

#include "Eigen/Dense"

class GroundTruth {
public:
  GroundTruth(const Eigen::VectorXd &value) :
      value_{value}
  {}

  inline const Eigen::VectorXd &getValue() const {
    return value_;
  }

private:
  Eigen::VectorXd value_;

};

#endif