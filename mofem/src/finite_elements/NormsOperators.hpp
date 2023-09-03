/** \file NormsOperators.hpp
  * \brief User data Operators for calculating norms

*/

#ifndef __NORM_OPERATORS_HPP__
#define __NORM_OPERATORS_HPP__

namespace MoFEM {

/** \brief Get norm of input VectorDouble for Tensor0
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

/** \brief Get norm of input MatrixDouble for Tensor1
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

/** \brief Get norm of input MatrixDouble for Tensor2
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

} // namespace MoFEM

#endif // __NORM_OPERATORS_HPP__
