/**
 * @file TimeScaling.hpp
 * @brief
 * @version 0.13.1
 * @date 2022-08-12
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _TIME_SCALING_HPP_
#define _TIME_SCALING_HPP_


namespace MoFEM {

struct ScalingMethod {

  /**
   * @brief Get scaling at given time
   *
   * @param time
   * @return double
   */
  virtual double getScale(const double time);

  ScalingMethod() = default;
  virtual ~ScalingMethod() = default;
};

/** \brief Force scale operator for reading two columns
 */
struct TimeScale : public ScalingMethod {

  TimeScale(std::string file_name = "-time_scalar_file",
            bool error_if_file_not_given = false);
  TimeScale(std::string file_name,
            char delimiter,
            bool error_if_file_not_given = false);
   double getScale(const double time);

private:
  MoFEMErrorCode timeData();
  double getScaleFromData(const double time);
  double getLinearScale(const double time);
  double lerp(const double a, const double b, double t);
  std::map<double, double> tSeries;
  int debug;
  std::string fileName;
  const char defaultDelimiter = ',';
  char delimiter = ',';
  bool errorIfFileNotGiven;
  std::function<double(double)> scalingMethod = [](double time){ return time; };
};
}

#endif //_TIME_SCALING_HPP_