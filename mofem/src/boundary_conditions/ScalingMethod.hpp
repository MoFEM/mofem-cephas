/**
 * @file ScalingMethod.hpp
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

  using ScalingFun = std::function<double(double)>;

  /**
   * @brief TimeScale constructor
   *
   * @param file_name Path to input CSV data file
   * @param error_if_file_not_given If file name is not provided, the
   * constructor will throw an error if this flag is set to true or throw a
   * warning and use linear scaling if this flag is set to false
   */
  TimeScale(
      std::string file_name = "", bool error_if_file_not_given = false,
      ScalingFun def_scaling_fun = [](double time) { return time; });

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
  TimeScale(
      std::string file_name, std::string delimiter,
        bool error_if_file_not_given = false,
      ScalingFun def_scaling_fun = [](double time) { return time; });

  /**
   * @brief Get scaling at a given time
   *
   * @param time
   * @return double
   */
  double getScale(const double time);

  std::string fileName = "";            //< file CSV data file
  PetscBool argFileScale = PETSC_FALSE; //< get file name from command line

private:
  MoFEMErrorCode timeData(std::string fileName, std::string delimiter);

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
  std::string fileNameFlag = "-time_scalar_file";

  static const std::string
      defaultDelimiter; // "(\\s*,\\s*|\\s+)"; 

  bool errorIfFileNotGiven;
  ScalingFun defScalingMethod = [](double time) { return time; };
  ScalingFun scalingMethod = defScalingMethod;
};

/** \brief Force scale operator for reading four columns (time and vector)
 */
template <int SPACE_DIM> struct TimeScaleVector : public ScalingMethod {

  TimeScaleVector(std::string name = "-time_vector_file",
                  bool error_if_file_not_given = false);

  TimeScaleVector(std::string name, int ms_id,
                  bool error_if_file_not_given = false);

  virtual FTensor::Tensor1<double, SPACE_DIM> getVector(const double time);
  virtual FTensor::Tensor1<double, SPACE_DIM>
  getVectorFromData(const double time);

private:
  MoFEMErrorCode timeData();

  std::map<double, FTensor::Tensor1<double, 3>> tSeries;
  int readFile, debug;
  string nAme;
  bool errorIfFileNotGiven;

  PetscBool fLg;
  std::function<FTensor::Tensor1<double, SPACE_DIM>(double)> scalingMethod =
      [this](double time) {
        FTensor::Index<'i', SPACE_DIM> i;
        FTensor::Tensor1<double, SPACE_DIM> s;
        s(i) = time;
        return s;
      };
};
using TimeScaleVector3 = TimeScaleVector<3>;
using TimeScaleVector2 = TimeScaleVector<3>;

} // namespace MoFEM

#endif //_TIME_SCALING_HPP_