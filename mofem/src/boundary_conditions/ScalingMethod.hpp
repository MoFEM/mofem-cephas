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

  /**
   * @brief TimeScale constructor
   *
   * @param file_name Path to input CSV data file
   * @param error_if_file_not_given If file name is not provided, the
   * constructor will throw an error if this flag is set to true or throw a
   * warning and use linear scaling if this flag is set to false
   */
  TimeScale(std::string file_name = "", bool error_if_file_not_given = false);

  /**
   * @brief TimeScale constructor
   *
   * @param file_name Path to input CSV data file
   * @param delimiter Character which is used to separate the data in a csv row,
   * by default it is ','
   * @param error_if_file_not_given If file name is not provided, the
   * constructor will throw an error if this flag is set to true or throw a
   * warning and use linear scaling if this flag is set to false
   */
  TimeScale(std::string file_name, char delimiter,
            bool error_if_file_not_given = false);

  /**
   * @brief Get scaling at a given time
   *
   * @param time
   * @return double
   */
  double getScale(const double time);

private:
  MoFEMErrorCode timeData();

  /**
   * @brief Get scaling at a given time when the scalar values have been
   * provided. Uses linear interpolation on the nearest time range to calculate
   * scaling if the provided time is not present in the data.
   * @return double
   */
  double getScaleFromData(const double time);

  /**
   * @brief Returns the value of time.
   * @return double
   */
  double getLinearScale(const double time);

  std::map<double, double> tSeries;
  std::string fileName = "";
  std::string fileNameFlag = "-time_scalar_file";
  static const char defaultDelimiter = ',';
  char delimiter = ',';
  bool errorIfFileNotGiven;
  std::function<double(double)> scalingMethod = [](double time) {
    return time;
  };
};
} // namespace MoFEM

#endif //_TIME_SCALING_HPP_