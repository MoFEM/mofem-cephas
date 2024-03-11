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

  OpCalcNormL2Tensor0(boost::shared_ptr<VectorDouble> data_ptr,
                      SmartPetscObj<Vec> data_vec, const int index,
                      boost::shared_ptr<VectorDouble> diff_data_ptr = nullptr);

  /**
   * \brief calculate norm of scalar values at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<VectorDouble> dataPtr;
  boost::shared_ptr<VectorDouble> diffDataPtr;
  SmartPetscObj<Vec> dataVec;
  const int iNdex;
};

/** \brief Get norm of input MatrixDouble for Tensor1
 *
 */
template <int DIM>
struct OpCalcNormL2Tensor1 : public ForcesAndSourcesCore::UserDataOperator {

  OpCalcNormL2Tensor1(boost::shared_ptr<MatrixDouble> data_ptr,
                      SmartPetscObj<Vec> data_vec, const int index,
                      boost::shared_ptr<MatrixDouble> diff_data_ptr = nullptr);

  /**
   * \brief calculate norm of vector values at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<MatrixDouble> dataPtr;
  boost::shared_ptr<MatrixDouble> diffDataPtr;
  SmartPetscObj<Vec> dataVec;
  const int iNdex;
};

/** \brief Get norm of input MatrixDouble for Tensor2
 *
 */
template <int DIM_1, int DIM_2>
struct OpCalcNormL2Tensor2 : public ForcesAndSourcesCore::UserDataOperator {

  OpCalcNormL2Tensor2(boost::shared_ptr<MatrixDouble> data_ptr,
                      SmartPetscObj<Vec> data_vec, const int index,
                      boost::shared_ptr<MatrixDouble> diff_data_ptr = nullptr);

  /**
   * \brief calculate norm of tensor values at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<MatrixDouble> dataPtr;
  boost::shared_ptr<MatrixDouble> diffDataPtr;
  SmartPetscObj<Vec> dataVec;
  const int iNdex;
};

typedef boost::function<double(const double, const double, const double)>
    ScalarFunc;

typedef boost::function<VectorDouble(const double, const double, const double)>
    VectorFunc;

/** \brief Get values from scalar function at integration points and save them
 * to VectorDouble for Tensor0
 *
 */

struct OpGetTensor0fromFunc : public ForcesAndSourcesCore::UserDataOperator {

  OpGetTensor0fromFunc(boost::shared_ptr<VectorDouble> data_ptr,
                       ScalarFunc scalar_function);

  /**
   * \brief calculate values of scalar function at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  ScalarFunc sFunc;
  boost::shared_ptr<VectorDouble> dataPtr;
};

template <int SPACE_DIM, int BASE_DIM>
struct OpGetTensor1fromFunc : public ForcesAndSourcesCore::UserDataOperator {

  OpGetTensor1fromFunc(boost::shared_ptr<MatrixDouble> data_ptr,
                       VectorFunc vector_function);

  /**
   * \brief calculate values of vector function at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  VectorFunc vFunc;
  boost::shared_ptr<MatrixDouble> dataPtr;
};

} // namespace MoFEM

#endif // __NORM_OPERATORS_HPP__
