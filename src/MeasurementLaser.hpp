#ifndef MEASUREMENT_LASER_HPP_
#define MEASUREMENT_LASER_HPP_

#include "Measurement.hpp"

class MeasurementLaser : public Measurement {
public:
  inline MeasurementLaser( const Eigen::VectorXd &measurement, long long timestamp )
      : Measurement{measurement, timestamp}
  {}

  inline bool valid() const override {
    if ( fabs(getMeasurement()[0]) < Measurement::eps &&
         fabs(getMeasurement()[1]) < Measurement::eps )
      return false;
    return true;
  }

};

#endif