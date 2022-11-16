/** \file FormsIntegrators.hpp
  * \brief Forms inteegrators
  * \ingroup mofem_form

*/



#ifndef __FORMS_INTEGRATORS_HPP__
#define __FORMS_INTEGRATORS_HPP__

namespace MoFEM {

//! [Storage and set boundary conditions]

struct EssentialBcStorage : public EntityStorage {
  EssentialBcStorage(VectorInt &indices) : entityIndices(indices) {}
  VectorInt entityIndices;
  using HashVectorStorage =
      map<std::string, std::vector<boost::shared_ptr<EssentialBcStorage>>>;
  static HashVectorStorage feStorage;
};

/**
 * @brief Set indices on entities on finite element
 * @ingroup mofem_forms
 *
 * If index is marked, its value is set to -1. DOF with such index is not
 * assembled into the system.
 *
 * Indices are stored on the entity.
 *
 */
struct OpSetBc : public ForcesAndSourcesCore::UserDataOperator {
  OpSetBc(std::string field_name, bool yes_set,
          boost::shared_ptr<std::vector<unsigned char>> boundary_marker);
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

public:
  bool yesSet;
  boost::shared_ptr<std::vector<unsigned char>> boundaryMarker;
};

struct OpUnSetBc : public ForcesAndSourcesCore::UserDataOperator {
  OpUnSetBc(std::string field_name);
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * @brief Set values to vector in operator
 * @ingroup mofem_forms
 *
 * @param V
 * @param data
 * @param ptr
 * @param iora
 * @return MoFEMErrorCode
 */
template <>
MoFEMErrorCode
VecSetValues<EssentialBcStorage>(Vec V,
                                 const EntitiesFieldData::EntData &data,
                                 const double *ptr, InsertMode iora);

/**
 * @brief Set values to matrix in operator
 *
 * @param M
 * @param row_data
 * @param col_data
 * @param ptr
 * @param iora
 * @return MoFEMErrorCode
 */
template <>
MoFEMErrorCode MatSetValues<EssentialBcStorage>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const double *ptr,
    InsertMode iora);

//! [Storage and set boundary conditions]

/**
 * @brief Form integrator assembly types
 * @ingroup mofem_forms
 *
 */
enum AssemblyType { PETSC, USER_ASSEMBLE, LAST_ASSEMBLE };

/**
 * @brief Form integrator integration types
 * @ingroup mofem_forms
 *
 */
enum IntegrationType { GAUSS, USER_INTEGRATION, LAST_INTEGRATION };

/**
 * @brief Scalar function type
 * @ingroup mofem_forms
 *
 */
using ScalarFun =
    boost::function<double(const double, const double, const double)>;

inline double scalar_fun_one(const double, const double, const double) {
  return 1;
}

/**
 * @brief Lambda function used to scale with time
 * 
 */
using TimeFun = boost::function<double(double)>;

/**
 * @brief Lambda function used to scale with time
 * 
 */
using FEFun = boost::function<double(const FEMethod *fe_ptr)>;

/**
 * @brief Constant function type
 *
 */
using ConstantFun = boost::function<double()>;

/**
 * @brief Vector function type
 * @ingroup mofem_forms
 *
 * @tparam DIM dimension of the return
 */
template <int DIM>
using VectorFun = boost::function<FTensor::Tensor1<double, DIM>(
    const double, const double, const double)>;

template <AssemblyType A, typename EleOp> struct OpBaseImpl : public EleOp {
  using OpType = typename EleOp::OpType;
  using EntData = EntitiesFieldData::EntData;

  OpBaseImpl(const std::string row_field_name, const std::string col_field_name,
             const OpType type, boost::shared_ptr<Range> ents_ptr = nullptr)
      : EleOp(row_field_name, col_field_name, type, false),
        assembleTranspose(false), onlyTranspose(false), entsPtr(ents_ptr) {}

  /**
   * \brief Do calculations for the left hand side
   * @param  row_side row side number (local number) of entity on element
   * @param  col_side column side number (local number) of entity on element
   * @param  row_type type of row entity 
   * @param  col_type type of column entity 
   * @param  row_data data for row
   * @param  col_data data for column
   * @return          error code
   */
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type, EntData &row_data,
                        EntData &col_data);

  /**
   * @brief Do calculations for the right hand side
   *
   * @param row_side
   * @param row_type
   * @param row_data
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode doWork(int row_side, EntityType row_type, EntData &row_data);

  TimeFun timeScalingFun; ///< assumes that time variable is set
  FEFun feScalingFun;     ///< assumes that time variable is set
  boost::shared_ptr<Range> entsPtr; ///< Entities on which element is run

protected:
  template <int DIM>
  inline FTensor::Tensor1<FTensor::PackPtr<double *, DIM>, DIM> getNf() {
    return getFTensor1FromArray<DIM, DIM>(locF);
  }

  template <int DIM>
  inline FTensor::Tensor2<FTensor::PackPtr<double *, DIM>, DIM, DIM>
  getLocMat(const int rr) {
    return getFTensor2FromArray<DIM, DIM, DIM>(locMat, rr);
  }

  int nbRows;             ///< number of dofs on rows
  int nbCols;             ///< number if dof on column
  int nbIntegrationPts;   ///< number of integration points
  int nbRowBaseFunctions; ///< number or row base functions

  int rowSide;           ///< row side number
  int colSide;           ///< column side number
  EntityType rowType;    ///< row type
  EntityType colType;    ///< column type

  bool assembleTranspose;
  bool onlyTranspose;

  MatrixDouble locMat;          ///< local entity block matrix
  MatrixDouble locMatTranspose; ///< local entity block matrix
  VectorDouble locF;            ///< local entity vector

  /**
   * \brief Integrate grad-grad operator
   * @param  row_data row data (consist base functions on row entity)
   * @param  col_data column data (consist base functions on column entity)
   * @return          error code
   */
  virtual MoFEMErrorCode iNtegrate(EntData &row_data, EntData &col_data) {
    return MOFEM_NOT_IMPLEMENTED;
  }

  virtual MoFEMErrorCode aSsemble(EntData &row_data, EntData &col_data,
                                  const bool trans) = 0;

  /**
   * \brief Class dedicated to integrate operator
   * @param  data entity data on element row
   * @return      error code
   */
  virtual MoFEMErrorCode iNtegrate(EntData &data) {
    return MOFEM_NOT_IMPLEMENTED;
  }

  virtual MoFEMErrorCode aSsemble(EntData &data) = 0;
};

/**
 * @brief Integrator forms
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 */
template <typename EleOp> struct FormsIntegrators {

  using EntData = EntitiesFieldData::EntData;
  using OpType = typename EleOp::OpType;

  /**
   * @brief Assembly methods
   * @ingroup mofem_forms
   *
   * @tparam A
   */
  template <AssemblyType A> struct Assembly {

    using OpBase = OpBaseImpl<A, EleOp>;

    /**
     * @brief Linear form
     * @ingroup mofem_forms
     *
     * @tparam I
     */
    template <IntegrationType I> struct LinearForm;

    /**
     * @brief Bi linear form
     * @ingroup mofem_forms
     *
     * @tparam I
     */
    template <IntegrationType I> struct BiLinearForm;

  }; // Assembly
};   // namespace MoFEM

template <AssemblyType A, typename EleOp>
MoFEMErrorCode
OpBaseImpl<A, EleOp>::doWork(int row_side, int col_side, EntityType row_type,
                             EntityType col_type,
                             EntitiesFieldData::EntData &row_data,
                             EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  if (entsPtr) {
    if (entsPtr->find(this->getFEEntityHandle()) == entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }

  // get number of dofs on row
  nbRows = row_data.getIndices().size();
  // if no dofs on row, exit that work, nothing to do here
  if (!nbRows)
    MoFEMFunctionReturnHot(0);
  rowSide = row_side;
  rowType = row_type;
  // get number of dofs on column
  nbCols = col_data.getIndices().size();
  // if no dofs on column, exit nothing to do here
  if (!nbCols)
    MoFEMFunctionReturnHot(0);
  colSide = col_side;
  colType = col_type;
  // get number of integration points
  nbIntegrationPts = EleOp::getGaussPts().size2();
  // get row base functions
  nbRowBaseFunctions = row_data.getN().size2();
  // set size of local entity bock
  locMat.resize(nbRows, nbCols, false);
  // clear matrix
  locMat.clear();
  // integrate local matrix for entity block
  CHKERR this->iNtegrate(row_data, col_data);

  // assemble local matrix
  auto check_if_assemble_transpose = [&] {
    if (this->sYmm) {
      if (row_side != col_side || row_type != col_type)
        return true;
      else 
        return false;
    } else if (assembleTranspose) {
      return true;
    }

    return false;
  };
  CHKERR aSsemble(row_data, col_data, check_if_assemble_transpose());
  MoFEMFunctionReturn(0);
}

template <AssemblyType A, typename EleOp>
MoFEMErrorCode OpBaseImpl<A, EleOp>::doWork(int row_side, EntityType row_type,
                                            EntData &row_data) {
  MoFEMFunctionBegin;

  if (entsPtr) {
    if (entsPtr->find(this->getFEEntityHandle()) == entsPtr->end())
      MoFEMFunctionReturnHot(0);
  }

  // get number of dofs on row
  nbRows = row_data.getIndices().size();
  rowSide = row_side;
  rowType = row_type;

  if (!nbRows)
    MoFEMFunctionReturnHot(0);
  // get number of integration points
  nbIntegrationPts = EleOp::getGaussPts().size2();
  // get row base functions
  nbRowBaseFunctions = row_data.getN().size2();
  // resize and clear the right hand side vector
  locF.resize(nbRows);
  locF.clear();
  // integrate local vector
  CHKERR this->iNtegrate(row_data);
  // assemble local vector
  CHKERR this->aSsemble(row_data);
  MoFEMFunctionReturn(0);
}

template <typename EleOp>
struct OpBaseImpl<PETSC, EleOp> : public OpBaseImpl<LAST_ASSEMBLE, EleOp> {
  using OpBaseImpl<LAST_ASSEMBLE, EleOp>::OpBaseImpl;

protected:
  MoFEMErrorCode aSsemble(EntitiesFieldData::EntData &row_data,
                          EntitiesFieldData::EntData &col_data,
                          const bool transpose) {
    MoFEMFunctionBegin;

    if (!this->timeScalingFun.empty())
      this->locMat *= this->timeScalingFun(this->getFEMethod()->ts_t);
    if (!this->feScalingFun.empty())
      this->locMat *= this->feScalingFun(this->getFEMethod());

    // Assemble transpose
    if (transpose) {
      this->locMatTranspose.resize(this->locMat.size2(), this->locMat.size1(),
                                   false);
      noalias(this->locMatTranspose) = trans(this->locMat);
      CHKERR MatSetValues<EssentialBcStorage>(
          this->getKSPB(), col_data, row_data,
          &*this->locMatTranspose.data().begin(), ADD_VALUES);
    }

    if (!this->onlyTranspose) {
      // assemble local matrix
      CHKERR MatSetValues<EssentialBcStorage>(
          this->getKSPB(), row_data, col_data, &*this->locMat.data().begin(),
          ADD_VALUES);
    }

    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode aSsemble(EntitiesFieldData::EntData &data) {

    if (!this->timeScalingFun.empty())
      this->locF *= this->timeScalingFun(this->getFEMethod()->ts_t);
    if (!this->feScalingFun.empty())
      this->locF *= this->feScalingFun(this->getFEMethod());

    return VecSetValues<EssentialBcStorage>(
        this->getKSPf(), data, &*this->locF.data().begin(), ADD_VALUES);
  }
};

} // namespace MoFEM

/**
 * \defgroup mofem_forms Forms Integrators
 *
 * \brief Classes and functions used to evaluate fields at integration pts,
 *jacobians, etc..
 *
 * \ingroup mofem_forces_and_sources
 **/

#endif //__FORMS_INTEGRATORS_HPP__