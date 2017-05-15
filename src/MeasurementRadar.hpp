#ifndef MEASUREMENT_RADAR_HPP_
#define MEASUREMENT_RADAR_HPP_

#include "Measurement.hpp"

class MeasurementRadar : public Measurement {
public:
  inline MeasurementRadar( const Eigen::VectorXd &measurement, long long timestamp )
      : Measurement{measurement, timestamp}
  {}

  inline virtual Eigen::VectorXd getCartesianMeasurement() const override
  {
    Eigen::VectorXd measurement(5);
    float rho = getMeasurement()[0];
    float phi = getMeasurement()[1];
    float x = rho * cos(phi);
    float y = rho * sin(phi);
    measurement << x, y, 0, 0, 0;
    return measurement;
  }

  inline bool valid() const override
  {
    if ( fabs(getMeasurement()[0]) < Measurement::eps )
      return false;
    return true;
  }




};

#endif