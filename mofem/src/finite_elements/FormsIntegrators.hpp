/** \file FormsIntegrators.hpp
  * \brief Forms inteegrators

*/

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __FORMS_INTEGRATORS_HPP__
#define __FORMS_INTEGRATORS_HPP__

namespace MoFEM {

//! [Storage and set boundary conditions]

struct EssentialBcStorage;

/**
 * @brief Set indices on entities on finite element
 *
 * If indices is marked, set its value to -1. DOF which such indice is not
 * assembled into system.
 *
 * Indices are strored on on entity.
 *
 */
struct OpSetBc : public ForcesAndSourcesCore::UserDataOperator {
  OpSetBc(std::string field_name, bool yes_set,
          boost::shared_ptr<std::vector<bool>> boundary_marker);
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

public:
  bool yesSet;
  boost::shared_ptr<std::vector<bool>> boundaryMarker;
};

struct OpUnSetBc : public ForcesAndSourcesCore::UserDataOperator {
  OpUnSetBc(std::string field_name);
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * @brief Set values to vector in operator
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
                                 const DataForcesAndSourcesCore::EntData &data,
                                 const double *ptr, InsertMode iora);

/**
 * @brief Set valyes to matrix in operator
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
    Mat M, const DataForcesAndSourcesCore::EntData &row_data,
    const DataForcesAndSourcesCore::EntData &col_data, const double *ptr,
    InsertMode iora);

//! [Storage and set boundary conditions]

enum AssemblyType { PETSC };
enum IntegrationType { GAUSS };

/**
 * @brief Integrator forms
 * @ingroup mofem_form
 * 
 * @tparam EleOp 
 */
template <typename EleOp> struct FormsIntegrators {

  using EntData = DataForcesAndSourcesCore::EntData;
  using OpType = typename EleOp::OpType;

  typedef boost::function<double(const double, const double, const double)>
      VectorFun;

  typedef boost::function<double(const double, const double, const double)>
      ScalarFun;

  template <AssemblyType A> struct OpBase;

  /**
   * @brief Assembly methods
   * @ingroup mofem_form
   * 
   * @tparam A 
   */
  template <AssemblyType A> struct Assembly {

    using OpBase = OpBase<A>;

    template <IntegrationType I> struct LinearForm;
    template <IntegrationType I> struct BiLinearForm;

  }; // Assembly
};   // namespace MoFEM

template <typename EleOp>
template <AssemblyType A>
struct FormsIntegrators<EleOp>::OpBase : public EleOp {

  OpBase(const std::string row_field_name, const std::string col_field_name,
         const OpType type)
      : EleOp(row_field_name, col_field_name, type, false) {}

  /**
   * \brief Do calculations for the left hand side
   * @param  row_side row side number (local number) of entity on element
   * @param  col_side column side number (local number) of entity on element
   * @param  row_type type of row entity MBVERTEX, MBEDGE, MBTRI or MBTET
   * @param  col_type type of column entity MBVERTEX, MBEDGE, MBTRI or MBTET
   * @param  row_data data for row
   * @param  col_data data for column
   * @return          error code
   */
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type, EntData &row_data,
                        EntData &col_data) {
    MoFEMFunctionBegin;
    // get number of dofs on row
    OpBase::nbRows = row_data.getIndices().size();
    // if no dofs on row, exit that work, nothing to do here
    if (!OpBase::nbRows)
      MoFEMFunctionReturnHot(0);
    // get number of dofs on column
    OpBase::nbCols = col_data.getIndices().size();
    // if no dofs on Columbia, exit nothing to do here
    if (!OpBase::nbCols)
      MoFEMFunctionReturnHot(0);
    // get number of integration points
    OpBase::nbIntegrationPts = OpBase::getGaussPts().size2();
    // get row base functions
    OpBase::nbRowBaseFunctions = row_data.getN().size2();
    // set size of local entity bock
    OpBase::locMat.resize(OpBase::nbRows, OpBase::nbCols, false);
    // clear matrix
    OpBase::locMat.clear();
    // integrate local matrix for entity block
    CHKERR this->iNtegrate(row_data, col_data);
    // assemble local matrix
    CHKERR this->aSsemble(row_data, col_data);
    MoFEMFunctionReturn(0);
  }

  /**
   * @brief Do calculations for the right hand side
   *
   * @param row_side
   * @param row_type
   * @param row_data
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode doWork(int row_side, EntityType row_type, EntData &row_data) {
    MoFEMFunctionBegin;
    // get number of dofs on row
    OpBase::nbRows = row_data.getIndices().size();
    if (!OpBase::nbRows)
      MoFEMFunctionReturnHot(0);
    // get number of integration points
    OpBase::nbIntegrationPts = OpBase::getGaussPts().size2();
    // get row base functions
    OpBase::nbRowBaseFunctions = row_data.getN().size2();
    // resize and clear the right hand side vector
    OpBase::locF.resize(nbRows);
    OpBase::locF.clear();
    // integrate local vector
    CHKERR this->iNtegrate(row_data);
    // assemble local vector
    CHKERR this->aSsemble(row_data);
    MoFEMFunctionReturn(0);
  }

protected:
  template <int DIM>
  inline FTensor::Tensor1<FTensor::PackPtr<double *, DIM>, DIM> getNf() {
    static_assert(DIM != DIM, "Not Implemented");
    return FTensor::Tensor1<FTensor::PackPtr<double *, DIM>, DIM>();
  }

  template <>
  inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 1> getNf() {
    return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 1>{&OpBase::locF[0]};
  }

  template <>
  inline FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> getNf() {
    return FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>{&OpBase::locF[0],
                                                              &OpBase::locF[1]};
  }

  template <>
  inline FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> getNf() {
    return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>{
        &OpBase::locF[0], &OpBase::locF[1], &OpBase::locF[2]};
  }

  template <int DIM>
  inline FTensor::Tensor2<FTensor::PackPtr<double *, DIM>, DIM, DIM>
  getLocMat(const int rr) {
    static_assert(DIM != DIM, "Not Implemented");
    return FTensor::Tensor2<FTensor::PackPtr<double *, DIM>, DIM, DIM>();
  }

  template <>
  inline FTensor::Tensor2<FTensor::PackPtr<double *, 2>, 2, 2>
  getLocMat(const int rr) {
    return FTensor::Tensor2<FTensor::PackPtr<double *, 2>, 2, 2>{
        &OpBase::locMat(rr + 0, 0), &OpBase::locMat(rr + 0, 1),
        &OpBase::locMat(rr + 1, 0), &OpBase::locMat(rr + 1, 2)};
  }

  template <>
  inline FTensor::Tensor2<FTensor::PackPtr<double *, 3>, 3, 3>
  getLocMat(const int rr) {
    return FTensor::Tensor2<FTensor::PackPtr<double *, 3>, 3, 3>{
        &OpBase::locMat(rr + 0, 0), &OpBase::locMat(rr + 0, 1),
        &OpBase::locMat(rr + 0, 2), &OpBase::locMat(rr + 1, 0),
        &OpBase::locMat(rr + 1, 1), &OpBase::locMat(rr + 1, 2),
        &OpBase::locMat(rr + 2, 0), &OpBase::locMat(rr + 2, 1),
        &OpBase::locMat(rr + 2, 2)};
  }

  int nbRows;             ///< number of dofs on rows
  int nbCols;             ///< number if dof on column
  int nbIntegrationPts;   ///< number of integration points
  int nbRowBaseFunctions; ///< number or row base functions

  MatrixDouble locMat; ///< local entity block matrix
  VectorDouble locF;   ///< local entity vector

  /**
   * \brief Integrate grad-grad operator
   * @param  row_data row data (consist base functions on row entity)
   * @param  col_data column data (consist base functions on column entity)
   * @return          error code
   */
  virtual MoFEMErrorCode iNtegrate(EntData &row_data, EntData &col_data) {
    return MOFEM_NOT_IMPLEMENTED;
  }

  template <AssemblyType T>
  MoFEMErrorCode aSsembleImpl(EntData &row_data, EntData &col_data) {
    return MOFEM_NOT_IMPLEMENTED;
  }

  template <>
  inline MoFEMErrorCode aSsembleImpl<PETSC>(EntData &row_data,
                                            EntData &col_data) {
    MoFEMFunctionBegin;
    // assemble local matrix
    CHKERR MatSetValues<EssentialBcStorage>(this->getKSPB(), row_data, col_data,
                                            &*locMat.data().begin(),
                                            ADD_VALUES);
    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode aSsemble(EntData &row_data, EntData &col_data) {
    return aSsembleImpl<A>(row_data, col_data);
  }

  /**
   * \brief Class dedicated to integrate operator
   * @param  data entity data on element row
   * @return      error code
   */
  virtual MoFEMErrorCode iNtegrate(EntData &data) {
    return MOFEM_NOT_IMPLEMENTED;
  }

  /**
   * \brief Class dedicated to assemble operator to global system vector
   * @param  data entity data (indices, base functions, etc. ) on element row
   * @return      error code
   */
  template <AssemblyType T> MoFEMErrorCode aSsembleImpl(EntData &data);

  template <>
  MoFEMErrorCode aSsembleImpl<PETSC>(FormsIntegrators<EleOp>::EntData &data) {
    MoFEMFunctionBegin;
    // get values from local vector
    const double *vals = &*locF.data().begin();
    // assemble vector
    CHKERR VecSetValues<EssentialBcStorage>(this->getKSPf(), data, vals,
                                            ADD_VALUES);
    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode aSsemble(EntData &data) { return aSsembleImpl<A>(data); }
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