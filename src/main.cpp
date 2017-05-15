
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "MeasurementLaser.hpp"
#include "MeasurementRadar.hpp"
#include "GroundTruth.hpp"
#include "Math.hpp"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[]) {

  check_arguments(argc, argv);

  string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  // column names for output file
  out_file_ << "time_stamp" << "\t";
  out_file_ << "px_state" << "\t";
  out_file_ << "py_state" << "\t";
  out_file_ << "v_state" << "\t";
  out_file_ << "yaw_angle_state" << "\t";
  out_file_ << "yaw_rate_state" << "\t";
  out_file_ << "sensor_type" << "\t";
  out_file_ << "NIS" << "\t";
  out_file_ << "px_measured" << "\t";
  out_file_ << "py_measured" << "\t";
  out_file_ << "px_ground_truth" << "\t";
  out_file_ << "py_ground_truth" << "\t";
  out_file_ << "vx_ground_truth" << "\t";
  out_file_ << "vy_ground_truth" << "\n";

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  // Create a UKF instance
  UKF ukf;

  // used to compute the RMSE later
  vector<VectorXd> estimation_vector;
  vector<VectorXd> ground_truth_vector;

  string line;

  while (getline(in_file_, line)) {
    string sensor_type;
    istringstream iss(line);
    long long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    VectorXd measurement_result(5);
    VectorXd cartesian_measurement;
    string measurement_type;
    double nis = 0.0;
    if (sensor_type.compare("L") == 0) {
      // LASER MEASUREMENT
      VectorXd measurement{2};
      float x;
      float y;
      iss >> x;
      iss >> y;
      measurement << x, y;
      iss >> timestamp;
      auto laser_measurement = MeasurementLaser(measurement, timestamp);
      if ( !ukf.processMeasurement(measurement_result, nis, &laser_measurement) )
        continue;
      cartesian_measurement = laser_measurement.getCartesianMeasurement();
      measurement_type = "lidar";
    } else if (sensor_type.compare("R") == 0) {
      // RADAR MEASUREMENT
      VectorXd measurement{3};
      float ro;
      float phi;
      float ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      measurement << ro, phi, ro_dot;
      iss >> timestamp;
      auto radar_measurement = MeasurementRadar(measurement, timestamp);
      if ( !ukf.processMeasurement(measurement_result, nis, &radar_measurement) )
        continue;
      cartesian_measurement = radar_measurement.getCartesianMeasurement();
      measurement_type = "radar";
    }
    // read ground truth data to compare later
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    VectorXd value{4};
    value << x_gt, y_gt, vx_gt, vy_gt;
    auto ground_truth = GroundTruth(value);

    // output the estimation
    out_file_ << measurement_result(0) << "\t";
    out_file_ << measurement_result(1) << "\t";
    out_file_ << measurement_result(2) << "\t";
    out_file_ << measurement_result(3) << "\t";
    out_file_ << measurement_result(4) << "\t";

    //output type
    out_file_ << measurement_type << "\t";

    //output NIS value
    out_file_ << nis << "\t";

    out_file_ << cartesian_measurement(0) << "\t";
    out_file_ << cartesian_measurement(1) << "\t";

    // output the ground truth packages
    out_file_ << ground_truth.getValue()(0) << "\t";
    out_file_ << ground_truth.getValue()(1) << "\t";
    out_file_ << ground_truth.getValue()(2) << "\t";
    out_file_ << ground_truth.getValue()(3) << "\n";

    estimation_vector.push_back(measurement_result);
    ground_truth_vector.push_back(ground_truth.getValue());
  }

  // compute the accuracy (RMSE)
  cout << "Accuracy - RMSE:" << endl << math::CalculateRMSE(estimation_vector, ground_truth_vector) << endl;

  // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  cout << "Done!" << endl;
  return 0;
}
