/** \file DataOperators.hpp

  Number of structures preforming operations on integration points.

  For example:
  - calculate Jacobian
  - calculate Piola-Transform on shape functions
  - calculate normals and tangent vectors at integration points on triangle

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

#ifndef __DATAOPERATORS_HPP
#define __DATAOPERATORS_HPP

using namespace boost::numeric;

namespace MoFEM {

/** \brief base operator to do operations at Gauss Pt. level
 * \ingroup mofem_forces_and_sources
 */
struct DataOperator {

  DataOperator(const bool symm = true);

  virtual ~DataOperator() = default;

  /** \brief Operator for bi-linear form, usually to calculate values on left
   * hand side
   */
  virtual MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                                EntityType col_type,
                                DataForcesAndSourcesCore::EntData &row_data,
                                DataForcesAndSourcesCore::EntData &col_data) {
    MoFEMFunctionBeginHot;
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "doWork function is not implemented for this operator");
    MoFEMFunctionReturnHot(0);
  }

  virtual MoFEMErrorCode opLhs(DataForcesAndSourcesCore &row_data,
                               DataForcesAndSourcesCore &col_data);

  /** \brief Operator for linear form, usually to calculate values on right hand
   * side
   */
  virtual MoFEMErrorCode doWork(int side, EntityType type,
                                DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBeginHot;
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "doWork function is not implemented for this operator");
    MoFEMFunctionReturnHot(0);
  }

  virtual MoFEMErrorCode opRhs(DataForcesAndSourcesCore &data,
                               const bool error_if_no_base = false);

  bool sYmm; ///< If true assume that matrix is symmetric structure

  std::array<bool, MBMAXTYPE>
      doEntities; ///< If true operator is executed for entity.

  // Deprecated variables. Use doEntities instead. I keep them for back
  // compatibility with some older modules. It will be removed in some future.

  bool &doVertices; ///< \deprectaed If false skip vertices
  bool &doEdges;    ///< \deprectaed If false skip edges
  bool &doQuads;    ///< \deprectaed
  bool &doTris;     ///< \deprectaed
  bool &doTets;     ///< \deprectaed
  bool &doPrisms;   ///< \deprectaed

  /**
   * \brief Get if operator uses symmetry of DOFs or not
   *
   * If symmetry is used, only not repeating combinations of entities are
   * looped.  For an example pair of (Vertex, Edge_0) and (Edge_0, Vertex) will
   * calculate the same matrices only transposed. Implementing that this can be
   * exploited by integrating only one pair.
   *
   * @return true if symmetry
   */
  inline bool getSymm() const { return sYmm; }

  /// set if operator is executed taking in account symmetry
  inline void setSymm() { sYmm = true; }

  /// unset if operator is executed for  non symmetric problem
  inline void unSetSymm() { sYmm = false; }

private:
  template <bool Symm>
  inline MoFEMErrorCode opLhs(DataForcesAndSourcesCore &row_data,
                              DataForcesAndSourcesCore &col_data);

  template <bool ErrorIfNoBase>
  inline MoFEMErrorCode opRhs(DataForcesAndSourcesCore &data,
                              const std::array<bool, MBMAXTYPE> &do_entities);
};

/**
 * \brief Transform local reference derivatives of shape function to global
 derivatives

 * \ingroup mofem_forces_and_sources
 */
struct OpSetInvJacH1 : public DataOperator {

  FTensor::Tensor2<double *, 3, 3> tInvJac;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  OpSetInvJacH1(MatrixDouble3by3 &inv_jac)
      : tInvJac(&inv_jac(0, 0), &inv_jac(0, 1), &inv_jac(0, 2), &inv_jac(1, 0),
                &inv_jac(1, 1), &inv_jac(1, 2), &inv_jac(2, 0), &inv_jac(2, 1),
                &inv_jac(2, 2)) {}

  MatrixDouble diffNinvJac;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief brief Transform local reference derivatives of shape function to
 global derivatives

 * \ingroup mofem_forces_and_sources
 */
struct OpSetInvJacHdivAndHcurl : public DataOperator {

  FTensor::Tensor2<double *, 3, 3> tInvJac;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  OpSetInvJacHdivAndHcurl(MatrixDouble3by3 &inv_jac)
      : tInvJac(&inv_jac(0, 0), &inv_jac(0, 1), &inv_jac(0, 2), &inv_jac(1, 0),
                &inv_jac(1, 1), &inv_jac(1, 2), &inv_jac(2, 0), &inv_jac(2, 1),
                &inv_jac(2, 2)) {}

  MatrixDouble diffHdivInvJac;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief transform local reference derivatives of shape function to global
 derivatives if higher order geometry is given

 * \ingroup mofem_forces_and_sources
*/
struct OpSetHoInvJacH1 : public DataOperator {

  MatrixDouble &invHoJac;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  OpSetHoInvJacH1(MatrixDouble &inv_ho_jac) : invHoJac(inv_ho_jac) {}

  MatrixDouble diffNinvJac;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief transform local reference derivatives of shape function to global
 derivatives if higher order geometry is given
 *

 * \ingroup mofem_forces_and_sources
*/
struct OpSetHoInvJacHdivAndHcurl : public DataOperator {

  MatrixDouble &invHoJac;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  OpSetHoInvJacHdivAndHcurl(MatrixDouble &inv_ho_jac) : invHoJac(inv_ho_jac) {}

  MatrixDouble diffHdivInvJac;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief apply contravariant (Piola) transfer to Hdiv space

Contravariant Piola transformation
\f[
\psi_i|_t = \frac{1}{\textrm{det}(J)}J_{ij}\hat{\psi}_j\\
\left.\frac{\partial \psi_i}{\partial \xi_j}\right|_t
=
\frac{1}{\textrm{det}(J)}J_{ik}\frac{\partial \hat{\psi}_k}{\partial \xi_j}
\f]

* \ingroup mofem_forces_and_sources

*/
struct OpSetContravariantPiolaTransform : public DataOperator {

  double &vOlume;

  FTensor::Tensor2<double *, 3, 3> tJac;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  OpSetContravariantPiolaTransform(double &volume, MatrixDouble3by3 &jac)
      : vOlume(volume),
        // jAc(jac),
        tJac(&jac(0, 0), &jac(0, 1), &jac(0, 2), &jac(1, 0), &jac(1, 1),
             &jac(1, 2), &jac(2, 0), &jac(2, 1), &jac(2, 2)) {}

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief Apply contravariant (Piola) transfer to Hdiv space for HO geometr

* \ingroup mofem_forces_and_sources
*/
struct OpSetHoContravariantPiolaTransform : public DataOperator {

  VectorDouble &detHoJac;
  MatrixDouble &hoJac;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  OpSetHoContravariantPiolaTransform(VectorDouble &det_jac, MatrixDouble &jac)
      : detHoJac(det_jac), hoJac(jac) {}

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief Apply covariant (Piola) transfer to Hcurl space for HO geometry
 */
struct OpSetHoCovariantPiolaTransform : public DataOperator {

  MatrixDouble &hoInvJac;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  OpSetHoCovariantPiolaTransform(MatrixDouble &inv_jac) : hoInvJac(inv_jac) {}

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief apply covariant transfer to Hcurl space

Contravariant Piola transformation
\f[
\psi_i|_t = \frac{1}{\textrm{det}(J)}J_{ij}\hat{\psi}_j\\
\left.\frac{\partial \psi_i}{\partial \xi_j}\right|_t
=
\frac{1}{\textrm{det}(J)}J_{ik}\frac{\partial \hat{\psi}_k}{\partial \xi_j}
\f]


* \ingroup mofem_forces_and_sources
*/
struct OpSetCovariantPiolaTransform : public DataOperator {

  FTensor::Tensor2<double *, 3, 3> tInvJac;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  OpSetCovariantPiolaTransform(MatrixDouble3by3 &inv_jac)
      : tInvJac(&inv_jac(0, 0), &inv_jac(0, 1), &inv_jac(0, 2), &inv_jac(1, 0),
                &inv_jac(1, 1), &inv_jac(1, 2), &inv_jac(2, 0), &inv_jac(2, 1),
                &inv_jac(2, 2)) {}

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief Get field values and gradients at Gauss points
 * \ingroup mofem_forces_and_sources
 */
template <int RANK, int DIM> struct OpGetDataAndGradient : public DataOperator {

  MatrixDouble &dataAtGaussPts;
  MatrixDouble &dataGradAtGaussPts;

  OpGetDataAndGradient(MatrixDouble &data_at_gauss_pt,
                       MatrixDouble &data_grad_at_gauss_pt)
      : dataAtGaussPts(data_at_gauss_pt),
        dataGradAtGaussPts(data_grad_at_gauss_pt) {}

  /**
   * Return tensor associated with matrix storing values
   */
  template <int R>
  FTensor::Tensor1<double *, R> getValAtGaussPtsTensor(MatrixDouble &data) {
    THROW_MESSAGE("Not implemented");
  }

  /**
   * Return tensor associated with matrix storing gradient values
   */
  template <int R, int D>
  FTensor::Tensor2<double *, R, D> getGradAtGaussPtsTensor(MatrixDouble &data) {
    THROW_MESSAGE("Not implemented");
  }

  /**
   * \brief Calculate gradient and values at integration points
   * @param  side side of entity on element
   * @param  type type of entity
   * @param  data data stored on entity (dofs values, dofs indices, etc.)
   * @return      error code
   */
  MoFEMErrorCode calculateValAndGrad(int side, EntityType type,
                                     DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBeginHot;
    const int nb_base_functions = data.getN().size2();
    bool constant_diff = false;
    if (type == MBVERTEX && data.getDiffN().size1() * data.getDiffN().size2() ==
                                DIM * nb_base_functions) {
      constant_diff = true;
    }
    const int nb_dofs = data.getFieldData().size();
    for (unsigned int gg = 0; gg < data.getN().size1(); gg++) {
      double *data_ptr, *n_ptr, *diff_n_ptr;
      n_ptr = &data.getN()(gg, 0);
      if (constant_diff) {
        diff_n_ptr = &data.getDiffN()(0, 0);
      } else {
        diff_n_ptr = &data.getDiffN()(gg, 0);
      }
      data_ptr = &*data.getFieldData().data().begin();
      for (int rr = 0; rr < RANK; rr++, data_ptr++) {
        dataAtGaussPts(gg, rr) +=
            cblas_ddot(nb_dofs / RANK, n_ptr, 1, data_ptr, RANK);
        double *diff_n_ptr2 = diff_n_ptr;
        for (unsigned int dd = 0; dd < DIM; dd++, diff_n_ptr2++) {
          dataGradAtGaussPts(gg, DIM * rr + dd) +=
              cblas_ddot(nb_dofs / RANK, diff_n_ptr2, DIM, data_ptr, RANK);
        }
      }
    }
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

    MoFEMFunctionBegin;

    if (data.getFieldData().size() == 0) {
      MoFEMFunctionReturnHot(0);
    }

    unsigned int nb_dofs = data.getFieldData().size();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);

    if (nb_dofs % RANK != 0) {
      SETERRQ4(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "data inconsistency, type %d, side %d, nb_dofs %d, rank %d",
               type, side, nb_dofs, RANK);
    }
    if (nb_dofs / RANK > data.getN().size2()) {
      std::cerr << side << " " << type << " "
                << ApproximationBaseNames[data.getBase()] << std::endl;
      std::cerr << data.getN() << std::endl;
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "data inconsistency nb_dofs >= data.N.size2(), i.e. %u >= %u",
               nb_dofs, data.getN().size2());
    }

    if (type == MBVERTEX) {
      dataAtGaussPts.resize(data.getN().size1(), RANK, false);
      dataGradAtGaussPts.resize(data.getN().size1(), RANK * DIM, false);
      dataAtGaussPts.clear();
      dataGradAtGaussPts.clear();
    }

    CHKERR calculateValAndGrad(side, type, data);

    MoFEMFunctionReturn(0);
  }
};

/**
 * \brief Specialization for field with 3 coefficients in 3 dimension
 */
template <>
template <>
FTensor::Tensor1<double *, 3>
OpGetDataAndGradient<3, 3>::getValAtGaussPtsTensor<3>(MatrixDouble &data);

/**
 * \brief Specialization for field with 3 coefficients in 3 dimension
 */
template <>
template <>
FTensor::Tensor2<double *, 3, 3>
OpGetDataAndGradient<3, 3>::getGradAtGaussPtsTensor<3, 3>(MatrixDouble &data);

/**
 * \brief Specialization for field with 3 coefficients in 3 dimension
 */
template <>
MoFEMErrorCode OpGetDataAndGradient<3, 3>::calculateValAndGrad(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data);

/**
 * \brief Specialization for field with for scalar field in 3 dimension
 */
template <>
MoFEMErrorCode OpGetDataAndGradient<1, 3>::calculateValAndGrad(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data);

/** \brief Calculate normals at Gauss points of triangle element
 * \ingroup mofem_forces_and_source
 */
struct OpGetCoordsAndNormalsOnFace : public DataOperator {

  MatrixDouble &cOords_at_GaussPt;
  MatrixDouble &nOrmals_at_GaussPt;
  MatrixDouble &tAngent1_at_GaussPt;
  MatrixDouble &tAngent2_at_GaussPt;

  OpGetCoordsAndNormalsOnFace(MatrixDouble &coords_at_gausspt,
                              MatrixDouble &normals_at_gausspt,
                              MatrixDouble &tangent1_at_gausspt,
                              MatrixDouble &tangent2_at_gausspt)
      : cOords_at_GaussPt(coords_at_gausspt),
        nOrmals_at_GaussPt(normals_at_gausspt),
        tAngent1_at_GaussPt(tangent1_at_gausspt),
        tAngent2_at_GaussPt(tangent2_at_gausspt) {}

  MatrixDouble sPin;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

  MoFEMErrorCode calculateNormals();
};

/** \brief calculate normals at Gauss points of triangle element
 * \ingroup mofem_forces_and_sources
 */
struct OpGetCoordsAndNormalsOnPrism : public DataOperator {

  MatrixDouble &cOords_at_GaussPtF3;
  MatrixDouble &nOrmals_at_GaussPtF3;
  MatrixDouble &tAngent1_at_GaussPtF3;
  MatrixDouble &tAngent2_at_GaussPtF3;
  MatrixDouble &cOords_at_GaussPtF4;
  MatrixDouble &nOrmals_at_GaussPtF4;
  MatrixDouble &tAngent1_at_GaussPtF4;
  MatrixDouble &tAngent2_at_GaussPtF4;

  OpGetCoordsAndNormalsOnPrism(
      MatrixDouble &coords_at_gaussptf3, MatrixDouble &normals_at_gaussptf3,
      MatrixDouble &tangent1_at_gaussptf3, MatrixDouble &tangent2_at_gaussptf3,
      MatrixDouble &coords_at_gaussptf4, MatrixDouble &normals_at_gaussptf4,
      MatrixDouble &tangent1_at_gaussptf4, MatrixDouble &tangent2_at_gaussptf4)
      : cOords_at_GaussPtF3(coords_at_gaussptf3),
        nOrmals_at_GaussPtF3(normals_at_gaussptf3),
        tAngent1_at_GaussPtF3(tangent1_at_gaussptf3),
        tAngent2_at_GaussPtF3(tangent2_at_gaussptf3),
        cOords_at_GaussPtF4(coords_at_gaussptf4),
        nOrmals_at_GaussPtF4(normals_at_gaussptf4),
        tAngent1_at_GaussPtF4(tangent1_at_gaussptf4),
        tAngent2_at_GaussPtF4(tangent2_at_gaussptf4) {}

  MatrixDouble sPin;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

  MoFEMErrorCode calculateNormals();
};

/** \brief transform Hdiv base fluxes from reference element to physical
 * triangle \ingroup mofem_forces_and_sources
 */
struct OpSetContravariantPiolaTransformOnFace : public DataOperator {

  const VectorDouble &nOrmal;
  const MatrixDouble &normalsAtGaussPts;

  OpSetContravariantPiolaTransformOnFace(const VectorDouble &normal,
                                         const MatrixDouble &normals_at_pts)
      : nOrmal(normal), normalsAtGaussPts(normals_at_pts) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief transform Hcurl base fluxes from reference element to physical
 * triangle \ingroup mofem_forces_and_sources
 */
struct OpSetCovariantPiolaTransformOnFace : public DataOperator {

  const VectorDouble &nOrmal;
  const MatrixDouble &normalsAtGaussPts;
  const VectorDouble &tAngent0;
  const MatrixDouble &tangent0AtGaussPt;
  const VectorDouble &tAngent1;
  const MatrixDouble &tangent1AtGaussPt;

  OpSetCovariantPiolaTransformOnFace(const VectorDouble &normal,
                                     const MatrixDouble &normals_at_pts,
                                     const VectorDouble &tangent0,
                                     const MatrixDouble &tangent0_at_pts,
                                     const VectorDouble &tangent1,
                                     const MatrixDouble &tangent1_at_pts)
      : nOrmal(normal), normalsAtGaussPts(normals_at_pts), tAngent0(tangent0),
        tangent0AtGaussPt(tangent0_at_pts), tAngent1(tangent1),
        tangent1AtGaussPt(tangent1_at_pts) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief Calculate tangent vector on edge form HO geometry approximation
 * \ingroup mofem_forces_and_sources
 */
struct OpGetHoTangentOnEdge : public DataOperator {

  MatrixDouble &tAngent;

  OpGetHoTangentOnEdge(MatrixDouble &tangent) : tAngent(tangent) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief transform Hcurl base fluxes from reference element to physical edge
 * \ingroup mofem_forces_and_sources
 */
struct OpSetCovariantPiolaTransformOnEdge : public DataOperator {

  const VectorDouble &tAngent;
  const MatrixDouble &tangentAtGaussPt;

  OpSetCovariantPiolaTransformOnEdge(const VectorDouble &tangent,
                                     const MatrixDouble &tangent_at_pts)
      : tAngent(tangent), tangentAtGaussPt(tangent_at_pts) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

} // namespace MoFEM

#endif //__DATAOPERATORS_HPP

/**
 * \defgroup mofem_forces_and_sources Forces and sources
 **/
