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

  /**
   * @brief  Scale vector
   *
   * Function is virtual, expected to be implemented at every derived class
   *
   * @param fe
   * @param Nf
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode scaleNf(const FEMethod *fe, VectorDouble &Nf) = 0;

  /**
   * @brief Apply scaling for series of scaling methods
   *
   * @param fe
   * @param methods_op
   * @param nf
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode
  applyScalingSeries(const FEMethod *fe,
                     boost::ptr_vector<ScalingMethod> &methods_op,
                     VectorDouble &nf);

  ScalingMethod() = default;
  virtual ~ScalingMethod() = default;
};

/** \brief Force scale operator for reading two columns
 */
struct TimeForceScale : public ScalingMethod {

  TimeForceScale(std::string name = "-time_data_file",
                 bool error_if_file_not_given = false);

   double getScale(const double time);

  /**
   * @brief Scale force the right hand vector
   *
   * @param fe
   * @param Nf
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode scaleNf(const FEMethod *fe, VectorDouble &Nf);

private:
  MoFEMErrorCode timeData();

  std::map<double, double> tSeries;
  int readFile, debug;
  string nAme;
  bool errorIfFileNotGiven;

  PetscBool fLg;
};

/**
 * @brief Read acceleration data
 *
 * Read data file with four columns, first column is time, and remaining columns
 * are acceleration in ax, ay, az direction.
 *
 */
struct TimeAccelerogram : public ScalingMethod {

  TimeAccelerogram(std::string name = "-accelerogram_data_file");

	/**
	 * @brief Add acceleration to vector Nf
	 * 
	 * @note Is assumed that Nf.size() == 3
	 * 
	 * @param fe 
	 * @param Nf vector with three elements
	 * @return MoFEMErrorCode 
	 */
  MoFEMErrorCode scaleNf(const FEMethod *fe, VectorDouble &Nf);

private:
  MoFEMErrorCode timeData();

  std::map<double, VectorDouble> tSeries;
  int readFile, debug;
  string nAme;
};
}

#endif //_TIME_SCALING_HPP_