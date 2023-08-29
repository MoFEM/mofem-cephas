/** \file NormsOperators.hpp
  * \brief User data Operators for calculating norms

*/

#ifndef __NORM_OPERATORS_HPP__
#define __NORM_OPERATORS_HPP__

namespace MoFEM {

/** \brief Get norm of input MatrixDouble
 *
 */
template <int FIELD_DIM>
struct OpCalcNormL2Tesnosr1 : public ForcesAndSourcesCore::UserDataOperator {

  OpCalcNormL2Tesnosr1(const std::string field_name,
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
