/** \file NormsOperators.hpp
  * \brief User data Operators for calculating norms

*/

#ifndef __NORM_OPERATORS_HPP___
#define __NORM_OPERATORS_HPP___

namespace MoFEM {

/** \brief Get norm of input VectorDouble for Tesnosr0
 *
 */

struct OpCalcNormL2Tensor0 : public ForcesAndSourcesCore::UserDataOperator {

  OpCalcNormL2Tensor0(const std::string field_name,
                      boost::shared_ptr<VectorDouble> data_ptr,
                      SmartPetscObj<Vec> data_vec, const int index);

  /**
   * \brief calculate values of scalar field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<VectorDouble> dataPtr;
  SmartPetscObj<Vec> dataVec;
  const int iNdex;
};

/** \brief Get norm of input MatrixDouble for Tesnosr1
 *
 */
template <int DIM>
struct OpCalcNormL2Tensor1 : public ForcesAndSourcesCore::UserDataOperator {

  OpCalcNormL2Tensor1(const std::string field_name,
                      boost::shared_ptr<MatrixDouble> data_ptr,
                      SmartPetscObj<Vec> data_vec, const int index);

  /**
   * \brief calculate values of scalar field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<MatrixDouble> dataPtr;
  SmartPetscObj<Vec> dataVec;
  const int iNdex;
};

/** \brief Get norm of input MatrixDouble for Tesnosr2
 *
 */
template <int DIM_1, int DIM_2>
struct OpCalcNormL2Tensor2 : public ForcesAndSourcesCore::UserDataOperator {

  OpCalcNormL2Tensor2(const std::string field_name,
                      boost::shared_ptr<MatrixDouble> data_ptr,
                      SmartPetscObj<Vec> data_vec, const int index);

  /**
   * \brief calculate values of scalar field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<MatrixDouble> dataPtr;
  SmartPetscObj<Vec> dataVec;
  const int iNdex;
};

/** \brief Get difference between input VectorDouble and VectorDouble
 * intended for Tensor0 fields
 *
 */

struct OpCalcDifferenceTensor0 : public ForcesAndSourcesCore::UserDataOperator {

  OpCalcDifferenceTensor0(const std::string field_name,
                          boost::shared_ptr<VectorDouble> data_1_ptr,
                          boost::shared_ptr<VectorDouble> data_2_ptr,
                          boost::shared_ptr<VectorDouble> result_ptr);

  /**
   * \brief calculate values of scalar field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<VectorDouble> data1Ptr;
  boost::shared_ptr<VectorDouble> data2Ptr;
  boost::shared_ptr<VectorDouble> resultPtr;
};

/** \brief Get difference between input MatrixDouble and MatrixDouble
 * intended for Tensor1 fields
 *
 */

template <int DIM>
struct OpCalcDifferenceTensor1 : public ForcesAndSourcesCore::UserDataOperator {

  OpCalcDifferenceTensor1(const std::string field_name,
                          boost::shared_ptr<MatrixDouble> data_1_ptr,
                          boost::shared_ptr<MatrixDouble> data_2_ptr,
                          boost::shared_ptr<MatrixDouble> result_ptr);

  /**
   * \brief calculate values of scalar field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<MatrixDouble> data1Ptr;
  boost::shared_ptr<MatrixDouble> data2Ptr;
  boost::shared_ptr<MatrixDouble> resultPtr;
};

/** \brief Get difference between input MatrixDouble and MatrixDouble
 * intended for Tensor1 fields
 *
 */

template <int DIM_1, int DIM_2>
struct OpCalcDifferenceTensor2 : public ForcesAndSourcesCore::UserDataOperator {

  OpCalcDifferenceTensor2(const std::string field_name,
                          boost::shared_ptr<MatrixDouble> data_1_ptr,
                          boost::shared_ptr<MatrixDouble> data_2_ptr,
                          boost::shared_ptr<MatrixDouble> result_ptr);

  /**
   * \brief calculate values of scalar field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<MatrixDouble> data1Ptr;
  boost::shared_ptr<MatrixDouble> data2Ptr;
  boost::shared_ptr<MatrixDouble> resultPtr;
};

} // namespace MoFEM

#endif // __NORM_OPERATORS_HPP___
