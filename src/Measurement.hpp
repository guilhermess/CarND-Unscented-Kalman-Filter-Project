#ifndef MEASUREMENT_HPP_
#define MEASUREMENT_HPP_

#include "Eigen/Dense"

class Measurement {
public:
  inline Measurement( const Eigen::VectorXd &measurement, long long timestamp )
      : measurement_{measurement}, timestamp_{timestamp}
  {}

  virtual ~Measurement() {}

  inline const Eigen::VectorXd &getMeasurement() const { return measurement_; }

  inline long long getTimestamp() const { return timestamp_; }

  inline virtual Eigen::VectorXd getCartesianMeasurement() const
  {
    Eigen::VectorXd measurement(5);
    measurement << measurement_(0), measurement_(1), 0, 0, 0;
    return measurement;
  }

  virtual bool valid() const  = 0;

protected:
  /*
   * Min value considered a measurement too small (i.e. too close to zero).
   */
  static const float eps;

private:
  Eigen::VectorXd measurement_;
  long long timestamp_;

};

#endif