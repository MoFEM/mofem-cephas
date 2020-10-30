/** \file FormsIntegrators.hpp
  * \brief Forms inteegrators
  * \ingroup mofem_form

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
 * @ingroup mofem_forms
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

/**
 * @brief Form integrator assembly tpes
 * @ingroup mofem_forms
 * 
 */
enum AssemblyType { PETSC, USER_ASSEMBLE, LAST_ASSEMBLE };

/**
 * @brief Fom integrayors inegrator types
 * @ingroup mofem_forms
 * 
 */
enum IntegrationType { GAUSS, USER_INTEGRATION, LAST_INTEGRATION };

/**
 * @brief Sacalr function type
 * @ingroup mofem_forms
 * 
 */
using ScalarFun =
    boost::function<double(const double, const double, const double)>;

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
  using EntData = DataForcesAndSourcesCore::EntData;

  OpBaseImpl(const std::string row_field_name, const std::string col_field_name,
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

  virtual MoFEMErrorCode aSsemble(EntData &row_data, EntData &col_data) = 0;

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

  using EntData = DataForcesAndSourcesCore::EntData;
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
                             DataForcesAndSourcesCore::EntData &row_data,
                             DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;
  // get number of dofs on row
  nbRows = row_data.getIndices().size();
  // if no dofs on row, exit that work, nothing to do here
  if (!nbRows)
    MoFEMFunctionReturnHot(0);
  // get number of dofs on column
  nbCols = col_data.getIndices().size();
  // if no dofs on Columbia, exit nothing to do here
  if (!nbCols)
    MoFEMFunctionReturnHot(0);
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
  CHKERR this->aSsemble(row_data, col_data);
  MoFEMFunctionReturn(0);
}

template <AssemblyType A, typename EleOp>
MoFEMErrorCode OpBaseImpl<A, EleOp>::doWork(int row_side, EntityType row_type,
                                            EntData &row_data) {
  MoFEMFunctionBegin;
  // get number of dofs on row
  nbRows = row_data.getIndices().size();
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
  MoFEMErrorCode aSsemble(DataForcesAndSourcesCore::EntData &row_data,
                          DataForcesAndSourcesCore::EntData &col_data) {
    // assemble local matrix
    return MatSetValues<EssentialBcStorage>(
        this->getKSPB(), row_data, col_data, &*this->locMat.data().begin(),
        ADD_VALUES);
  }

  MoFEMErrorCode aSsemble(DataForcesAndSourcesCore::EntData &data) {
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