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

  virtual ~DataOperator() {}

  /** \brief Operator for bi-linear form, usually to calculate values on right hand side
    */
  virtual PetscErrorCode doWork(
    int row_side,int col_side,
    EntityType row_type,EntityType col_type,
    DataForcesAndSurcesCore::EntData &row_data,
    DataForcesAndSurcesCore::EntData &col_data) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode opLhs(
    DataForcesAndSurcesCore &row_data,
    DataForcesAndSurcesCore &col_data,
    bool symm = true
  );

  /** \brief Operator for linear form, usually to calculate values on left hand side
    */
  virtual PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode opRhs(
    DataForcesAndSurcesCore &data,
    const bool do_vertices = true,
    const bool do_edges = true,
    const bool do_quads = true,
    const bool do_tris = true,
    const bool do_tets = true,
    const bool do_prisms = true,
    const bool er_ror_if_no_base = true
  );

};

/**
 * \brief Calculate inverse of tensor rank 2 at integration points

 * \ingroup mofem_forces_and_sources
 */
template<int Tensor_Dim,class T,class L,class A>
inline PetscErrorCode invertTensor3by3(
  ublas::matrix<T,L,A> &jac_data,
  ublas::vector<T,A> &det_data,
  ublas::matrix<T,L,A> &inv_jac_data
) {
  PetscFunctionBegin;
  SETERRQ(
    PETSC_COMM_SELF,
    MOFEM_NOT_IMPLEMENTED,
    "Specialization for this template not yet implemented"
  );
  PetscFunctionReturn(0);
}

template<>
inline PetscErrorCode invertTensor3by3<3,double,ublas::row_major,ublas::unbounded_array<double> >(
  MatrixDouble &jac_data,
  VectorDouble &det_data,
  MatrixDouble &inv_jac_data
);

/**
 * \brief Calculate determinant

 * \ingroup mofem_forces_and_sources
 */
template<class T1,class T2>
inline PetscErrorCode determinantTensor3by3(
  T1 &t,T2 &det
) {
  PetscFunctionBegin;
  det =
    +t(0,0)*t(1,1)*t(2,2) + t(1,0)*t(2,1)*t(0,2)
    +t(2,0)*t(0,1)*t(1,2) - t(0,0)*t(2,1)*t(1,2)
    -t(2,0)*t(1,1)*t(0,2) - t(1,0)*t(0,1)*t(2,2);
  PetscFunctionReturn(0);
}

/**
 * \brief Calculate matrix inverse

 * \ingroup mofem_forces_and_sources
 */
template<class T1,class T2,class T3>
inline PetscErrorCode invertTensor3by3(
  T1 &t,T2 &det,T3 &inv_t
) {
  PetscFunctionBegin;
  inv_t(0,0) = (t(1,1)*t(2,2)-t(1,2)*t(2,1))/det;
  inv_t(0,1) = (t(0,2)*t(2,1)-t(0,1)*t(2,2))/det;
  inv_t(0,2) = (t(0,1)*t(1,2)-t(0,2)*t(1,1))/det;
  inv_t(1,0) = (t(1,2)*t(2,0)-t(1,0)*t(2,2))/det;
  inv_t(1,1) = (t(0,0)*t(2,2)-t(0,2)*t(2,0))/det;
  inv_t(1,2) = (t(0,2)*t(1,0)-t(0,0)*t(1,2))/det;
  inv_t(2,0) = (t(1,0)*t(2,1)-t(1,1)*t(2,0))/det;
  inv_t(2,1) = (t(0,1)*t(2,0)-t(0,0)*t(2,1))/det;
  inv_t(2,2) = (t(0,0)*t(1,1)-t(0,1)*t(1,0))/det;
  PetscFunctionReturn(0);
}

/**
 * \brief Specialization for symmetric tensor

 * \ingroup mofem_forces_and_sources
 */
template<>
inline PetscErrorCode invertTensor3by3<
FTensor::Tensor2_symmetric<double,3>,double,FTensor::Tensor2_symmetric<double,3>
>(
  FTensor::Tensor2_symmetric<double,3> &t,
  double &det,
  FTensor::Tensor2_symmetric<double,3> &inv_t
) {
  PetscFunctionBegin;
  inv_t(0,0) = (t(1,1)*t(2,2)-t(1,2)*t(2,1))/det;
  inv_t(0,1) = (t(0,2)*t(2,1)-t(0,1)*t(2,2))/det;
  inv_t(0,2) = (t(0,1)*t(1,2)-t(0,2)*t(1,1))/det;
  inv_t(1,1) = (t(0,0)*t(2,2)-t(0,2)*t(2,0))/det;
  inv_t(1,2) = (t(0,2)*t(1,0)-t(0,0)*t(1,2))/det;
  inv_t(2,2) = (t(0,0)*t(1,1)-t(0,1)*t(1,0))/det;
  PetscFunctionReturn(0);
}

/**
 * \brief Specialization for symmetric (pointer) tensor

 * \ingroup mofem_forces_and_sources
 */
template<>
inline PetscErrorCode invertTensor3by3<
FTensor::Tensor2_symmetric<double,3>,double,FTensor::Tensor2_symmetric<double*,3>
>(
  FTensor::Tensor2_symmetric<double,3> &t,
  double &det,
  FTensor::Tensor2_symmetric<double*,3> &inv_t
) {
  PetscFunctionBegin;
  inv_t(0,0) = (t(1,1)*t(2,2)-t(1,2)*t(2,1))/det;
  inv_t(0,1) = (t(0,2)*t(2,1)-t(0,1)*t(2,2))/det;
  inv_t(0,2) = (t(0,1)*t(1,2)-t(0,2)*t(1,1))/det;
  inv_t(1,1) = (t(0,0)*t(2,2)-t(0,2)*t(2,0))/det;
  inv_t(1,2) = (t(0,2)*t(1,0)-t(0,0)*t(1,2))/det;
  inv_t(2,2) = (t(0,0)*t(1,1)-t(0,1)*t(1,0))/det;
  PetscFunctionReturn(0);
}

/**
 * \brief Transform local reference derivatives of shape function to global derivatives

 * \ingroup mofem_forces_and_sources
 */
struct OpSetInvJacH1: public DataOperator {

  FTensor::Tensor2<double*,3,3> tInvJac;
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;

  OpSetInvJacH1(MatrixDouble3by3 &inv_jac):
  tInvJac(
    &inv_jac(0,0),&inv_jac(0,1),&inv_jac(0,2),
    &inv_jac(1,0),&inv_jac(1,1),&inv_jac(1,2),
    &inv_jac(2,0),&inv_jac(2,1),&inv_jac(2,2)
  ) {
  }

  MatrixDouble diffNinvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data
  );

};

/**
 * \brief brief Transform local reference derivatives of shape function to global derivatives

 * \ingroup mofem_forces_and_sources
 */
struct OpSetInvJacHdivAndHcurl: public DataOperator {

  FTensor::Tensor2<double*,3,3> tInvJac;
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',3> k;

  OpSetInvJacHdivAndHcurl(MatrixDouble3by3 &inv_jac):
  tInvJac(
    &inv_jac(0,0),&inv_jac(0,1),&inv_jac(0,2),
    &inv_jac(1,0),&inv_jac(1,1),&inv_jac(1,2),
    &inv_jac(2,0),&inv_jac(2,1),&inv_jac(2,2)
  ) {
  }

  MatrixDouble diffHdivInvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/**
 * \brief transform local reference derivatives of shape function to global derivatives if higher order geometry is given

 * \ingroup mofem_forces_and_sources
*/
struct OpSetHoInvJacH1: public DataOperator {

  MatrixDouble &invHoJac;
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  OpSetHoInvJacH1(MatrixDouble &inv_ho_jac): invHoJac(inv_ho_jac) {}

  MatrixDouble diffNinvJac;
  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/**
 * \brief transform local reference derivatives of shape function to global derivatives if higher order geometry is given
 *

 * \ingroup mofem_forces_and_sources
*/
struct OpSetHoInvJacHdivAndHcurl: public DataOperator {

  MatrixDouble &invHoJac;
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',3> k;

  OpSetHoInvJacHdivAndHcurl(MatrixDouble &inv_ho_jac): invHoJac(inv_ho_jac) {}

  MatrixDouble diffHdivInvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data
  );

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
struct OpSetContravariantPiolaTransform: public DataOperator {

  double &vOlume;
  // MatrixDouble3by3 &jAc;

  FTensor::Tensor2<double*,3,3> tJac;
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',3> k;

  OpSetContravariantPiolaTransform(double &volume,MatrixDouble3by3 &jac):
  vOlume(volume),
  // jAc(jac),
  tJac(
    &jac(0,0),&jac(0,1),&jac(0,2),
    &jac(1,0),&jac(1,1),&jac(1,2),
    &jac(2,0),&jac(2,1),&jac(2,2)
  ) {
  }

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;

  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief Apply contravariant (Piola) transfer to Hdiv space for HO geometr

* \ingroup mofem_forces_and_sources
*/
struct OpSetHoContravariantPiolaTransform: public DataOperator {

  VectorDouble &detHoJac;
  MatrixDouble &hoJac;
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',3> k;

  OpSetHoContravariantPiolaTransform(VectorDouble &det_jac,MatrixDouble &jac):
  detHoJac(det_jac),hoJac(jac) {}

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;
  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief Apply covariant (Piola) transfer to Hcurl space for HO geometry
*/
struct OpSetHoCovariantPiolaTransform: public DataOperator {

  MatrixDouble &hoInvJac;
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',3> k;

  OpSetHoCovariantPiolaTransform(MatrixDouble &inv_jac):
  hoInvJac(inv_jac) {}

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;
  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

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
struct OpSetCovariantPiolaTransform: public DataOperator {


  FTensor::Tensor2<double*,3,3> tInvJac;
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',3> k;

  OpSetCovariantPiolaTransform(MatrixDouble3by3 &inv_jac):
  tInvJac(
    &inv_jac(0,0),&inv_jac(0,1),&inv_jac(0,2),
    &inv_jac(1,0),&inv_jac(1,1),&inv_jac(1,2),
    &inv_jac(2,0),&inv_jac(2,1),&inv_jac(2,2)
  ) {
  }

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;

  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};


/** \brief Get field values and gradients at Gauss points
  * \ingroup mofem_forces_and_sources
  */
struct OpGetDataAndGradient: public DataOperator {

  MatrixDouble &data_at_GaussPt;
  MatrixDouble &dataGrad_at_GaussPt;

  const unsigned int dIm;
  const FieldCoefficientsNumber rAnk;

  OpGetDataAndGradient(
    MatrixDouble &data_at_gauss_pt,
    MatrixDouble &data_grad_at_gauss_pt,
    FieldCoefficientsNumber rank,
    int dim = 3
  ):
  data_at_GaussPt(data_at_gauss_pt),
  dataGrad_at_GaussPt(data_grad_at_gauss_pt),
  dIm(dim),
  rAnk(rank) {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
  );

};

/** \brief Calculate normals at Gauss points of triangle element
  * \ingroup mofem_forces_and_source
  */
struct OpGetCoordsAndNormalsOnFace: public DataOperator {

  MatrixDouble &cOords_at_GaussPt;
  MatrixDouble &nOrmals_at_GaussPt;
  MatrixDouble &tAngent1_at_GaussPt;
  MatrixDouble &tAngent2_at_GaussPt;

  OpGetCoordsAndNormalsOnFace(
    MatrixDouble &coords_at_gausspt,
    MatrixDouble &normals_at_gausspt,
    MatrixDouble &tangent1_at_gausspt,
    MatrixDouble &tangent2_at_gausspt):
    cOords_at_GaussPt(coords_at_gausspt),
    nOrmals_at_GaussPt(normals_at_gausspt),
    tAngent1_at_GaussPt(tangent1_at_gausspt),
    tAngent2_at_GaussPt(tangent2_at_gausspt) {}

  MatrixDouble sPin;
  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
  );

  PetscErrorCode calculateNormals();

};

/** \brief calculate normals at Gauss points of triangle element
  * \ingroup mofem_forces_and_sources
  */
struct OpGetCoordsAndNormalsOnPrism: public DataOperator {

  MatrixDouble &cOords_at_GaussPtF3;
  MatrixDouble &nOrmals_at_GaussPtF3;
  MatrixDouble &tAngent1_at_GaussPtF3;
  MatrixDouble &tAngent2_at_GaussPtF3;
  MatrixDouble &cOords_at_GaussPtF4;
  MatrixDouble &nOrmals_at_GaussPtF4;
  MatrixDouble &tAngent1_at_GaussPtF4;
  MatrixDouble &tAngent2_at_GaussPtF4;

  OpGetCoordsAndNormalsOnPrism(
    MatrixDouble &coords_at_gaussptf3,
    MatrixDouble &normals_at_gaussptf3,
    MatrixDouble &tangent1_at_gaussptf3,
    MatrixDouble &tangent2_at_gaussptf3,
    MatrixDouble &coords_at_gaussptf4,
    MatrixDouble &normals_at_gaussptf4,
    MatrixDouble &tangent1_at_gaussptf4,
    MatrixDouble &tangent2_at_gaussptf4):
    cOords_at_GaussPtF3(coords_at_gaussptf3),
    nOrmals_at_GaussPtF3(normals_at_gaussptf3),
    tAngent1_at_GaussPtF3(tangent1_at_gaussptf3),
    tAngent2_at_GaussPtF3(tangent2_at_gaussptf3),
    cOords_at_GaussPtF4(coords_at_gaussptf4),
    nOrmals_at_GaussPtF4(normals_at_gaussptf4),
    tAngent1_at_GaussPtF4(tangent1_at_gaussptf4),
    tAngent2_at_GaussPtF4(tangent2_at_gaussptf4) {}

  MatrixDouble sPin;
  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
  );

  PetscErrorCode calculateNormals();

};

/** \brief transform Hdiv base fluxes from reference element to physical triangle
 * \ingroup mofem_forces_and_sources
 */
struct OpSetContravariantPiolaTransoformOnTriangle: public DataOperator {

  const VectorDouble &nOrmal;
  const MatrixDouble &normalsAtGaussPt;

  OpSetContravariantPiolaTransoformOnTriangle(
    const VectorDouble &normal,
    const MatrixDouble &normals_at_pts
  ):
  nOrmal(normal),
  normalsAtGaussPt(normals_at_pts) {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
  );

};

/** \brief transform Hcurl base fluxes from reference element to physical triangle
 * \ingroup mofem_forces_and_sources
 */
struct OpSetCovariantPiolaTransoformOnTriangle: public DataOperator {

  const VectorDouble &nOrmal;
  const MatrixDouble &normalsAtGaussPt;
  const VectorDouble &tAngent0;
  const MatrixDouble &tangent0AtGaussPt;
  const VectorDouble &tAngent1;
  const MatrixDouble &tangent1AtGaussPt;

  OpSetCovariantPiolaTransoformOnTriangle(
    const VectorDouble &normal,
    const MatrixDouble &normals_at_pts,
    const VectorDouble &tangent0,
    const MatrixDouble &tangent0_at_pts,
    const VectorDouble &tangent1,
    const MatrixDouble &tangent1_at_pts
  ):
  nOrmal(normal),
  normalsAtGaussPt(normals_at_pts),
  tAngent0(tangent0),
  tangent0AtGaussPt(tangent0_at_pts),
  tAngent1(tangent1),
  tangent1AtGaussPt(tangent1_at_pts)
  {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
  );

};

/** \brief Calculate tangent vector on edge form HO geometry approximation
 * \ingroup mofem_forces_and_sources
 */
struct OpGetHoTangentOnEdge: public DataOperator {

  MatrixDouble &tAngent;

  OpGetHoTangentOnEdge(MatrixDouble &tangent):
    tAngent(tangent) {}

  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief transform Hcurl base fluxes from reference element to physical edge
 * \ingroup mofem_forces_and_sources
 */
struct OpSetCovariantPiolaTransoformOnEdge: public DataOperator {

  const VectorDouble &tAngent;
  const MatrixDouble &tangentAtGaussPt;

  OpSetCovariantPiolaTransoformOnEdge(
    const VectorDouble &tangent,
    const MatrixDouble &tangent_at_pts
  ):
  tAngent(tangent),
  tangentAtGaussPt(tangent_at_pts)
  {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
  );

};

}

#endif //__DATAOPERATORS_HPP

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources Forces and sources
 ******************************************************************************/
