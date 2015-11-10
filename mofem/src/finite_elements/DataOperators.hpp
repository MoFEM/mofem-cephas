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

  virtual PetscErrorCode opLhs(DataForcesAndSurcesCore &row_data,DataForcesAndSurcesCore &col_data,bool symm = true);


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

  virtual PetscErrorCode opRhs(DataForcesAndSurcesCore &data);

};

/// \brief Transform local reference derivatives of shape function to global derivatives
struct OpSetInvJacH1: public DataOperator {

  MatrixDouble &invJac;
  OpSetInvJacH1(MatrixDouble &_invJac): invJac(_invJac) {}

  MatrixDouble diffNinvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/// \brief Transform local reference derivatives of shape function to global derivatives
struct OpSetInvJacHdiv: public DataOperator {

  MatrixDouble &invJac;
  OpSetInvJacHdiv(MatrixDouble &_invJac): invJac(_invJac) {}

  MatrixDouble diffHdiv_invJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief transform local reference derivatives of shape function to global derivatives if higher order geometry is given
  */
struct OpSetHoInvJacH1: public DataOperator {

  MatrixDouble &invHoJac;
  OpSetHoInvJacH1(MatrixDouble &inv_ho_jac): invHoJac(inv_ho_jac) {}

  MatrixDouble diffNinvJac;
  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};


/** \brief transform local reference derivatives of shape function to global derivatives if higher order geometry is given
  */
struct OpSetHoInvJacHdiv: public DataOperator {

  MatrixDouble &invHoJac;
  OpSetHoInvJacHdiv(MatrixDouble &inv_ho_jac): invHoJac(inv_ho_jac) {}

  MatrixDouble diffHdiv_invJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief apply covariant (Piola) transfer for Hdiv space
*/
struct OpSetPiolaTransform: public DataOperator {

  double &vOlume;
  MatrixDouble &Jac;
  OpSetPiolaTransform(double &_vOlume,MatrixDouble &_Jac):
  vOlume(_vOlume),Jac(_Jac) {}

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;
  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief Apply covariant (Piola) transfer for Hdiv space for HO geometry
*/
struct OpSetHoPiolaTransform: public DataOperator {

  VectorDouble &detHoJac;
  MatrixDouble &hoJac;
  OpSetHoPiolaTransform(VectorDouble &det_jac,MatrixDouble &jac):
  detHoJac(det_jac),hoJac(jac) {}

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

  const unsigned int dim;
  const ApproximationRank rank;

  OpGetDataAndGradient(
    MatrixDouble &data_at_gauss_pt,
    MatrixDouble &data_grad_at_gauss_pt,
    ApproximationRank _rank,
    unsigned int _dim = 3):
      data_at_GaussPt(data_at_gauss_pt),
      dataGrad_at_GaussPt(data_grad_at_gauss_pt),
      dim(_dim),
      rank(_rank) {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data);

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


/** \brief transform Hdiv space fluxes from reference element to physical triangle
 */
struct OpSetPiolaTransoformOnTriangle: public DataOperator {

  const VectorDouble &normal;
  const MatrixDouble &nOrmals_at_GaussPt;

  OpSetPiolaTransoformOnTriangle(
    const VectorDouble &_normal,
    const MatrixDouble &_nOrmals_at_GaussPt):
    normal(_normal),nOrmals_at_GaussPt(_nOrmals_at_GaussPt) {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data);

};

/** \brief Calculate tangent vector on edge form HO geometry approximation
 */
struct OpGetHoTangentOnEdge: public DataOperator {

  MatrixDouble &tAngent;

  OpGetHoTangentOnEdge(MatrixDouble &tangent):
    tAngent(tangent) {}

  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};


}

#endif //__DATAOPERATORS_HPP

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources Forces and sources
 ******************************************************************************/
