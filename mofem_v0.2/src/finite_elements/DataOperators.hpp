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

  /** \brief operator for linear form, usually to calculate values on right hand side
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


  /** \brief operator for linear form, usually to calculate values on left hand side
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

/// \brief transform local reference derivatives of shape function to global derivatives
struct OpSetInvJacH1: public DataOperator {

  ublas::matrix<double> &invJac;
  OpSetInvJacH1(ublas::matrix<double> &_invJac): invJac(_invJac) {}

  ublas::matrix<FieldData> diffNinvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/// \brief transform local reference derivatives of shape function to global derivatives
struct OpSetInvJacHdiv: public DataOperator {

  ublas::matrix<double> &invJac;
  OpSetInvJacHdiv(ublas::matrix<double> &_invJac): invJac(_invJac) {}

  ublas::matrix<FieldData> diffHdiv_invJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief transform local reference derivatives of shape function to global derivatives if higher order geometry is given
  */
struct OpSetHoInvJacH1: public DataOperator {

  ublas::matrix<double> &invHoJac;
  OpSetHoInvJacH1(ublas::matrix<double> &_invHoJac): invHoJac(_invHoJac) {}

  ublas::matrix<FieldData> diffNinvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};


/** \brief transform local reference derivatives of shape function to global derivatives if higher order geometry is given
  */
struct OpSetHoInvJacHdiv: public DataOperator {

  ublas::matrix<double> &invHoJac;
  OpSetHoInvJacHdiv(ublas::matrix<double> &_invHoJac): invHoJac(_invHoJac) {}

  ublas::matrix<FieldData> diffHdiv_invJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief apply covariant (Piola) transfer for Hdiv space
  */
struct OpSetPiolaTransform: public DataOperator {

    double &vOlume;
    ublas::matrix<double> &Jac;
    OpSetPiolaTransform(double &_vOlume,ublas::matrix<double> &_Jac):
      vOlume(_vOlume),Jac(_Jac) {}

    ublas::matrix<FieldData> piolaN;
    ublas::matrix<FieldData> piolaDiffN;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief apply covariant (Piola) transfer for Hdiv space for HO geometry
  */
struct OpSetHoPiolaTransform: public DataOperator {

    ublas::vector<double> &detHoJac;
    ublas::matrix<double> &hoJac;
    OpSetHoPiolaTransform(ublas::vector<double> &_detJac,ublas::matrix<double> &_Jac):
      detHoJac(_detJac),hoJac(_Jac) {}

    ublas::matrix<FieldData> piolaN;
    ublas::matrix<FieldData> piolaDiffN;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};


/** \brief operator to calculate function values and its gradients at Gauss points
  * \ingroup mofem_forces_and_sources
  */
struct OpGetData: public DataOperator {

  ublas::matrix<FieldData> &data_at_GaussPt;
  ublas::matrix<FieldData> &dataGrad_at_GaussPt;

  const unsigned int dim;
  const ApproximationRank rank;

  OpGetData(
    ublas::matrix<FieldData> &_data_at_GaussPt,
    ublas::matrix<FieldData> &_dataGrad_at_GaussPt,
    ApproximationRank _rank,unsigned int _dim = 3):
      data_at_GaussPt(_data_at_GaussPt),
      dataGrad_at_GaussPt(_dataGrad_at_GaussPt),
      dim(_dim),rank(_rank) {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data);

};

/** \brief calculate normals at Gauss points of triangle element
  * \ingroup mofem_forces_and_sources
  */
struct OpGetNormals: public DataOperator {

  ublas::matrix<FieldData> &nOrmals_at_GaussPt;
  ublas::matrix<FieldData> &tAngent1_at_GaussPt;
  ublas::matrix<FieldData> &tAngent2_at_GaussPt;

  OpGetNormals(
    ublas::matrix<FieldData> &_nOrmals_at_GaussPt,
    ublas::matrix<FieldData> &_tAngent1_at_GaussPt,
    ublas::matrix<FieldData> &_tAngent2_at_GaussPt):
    nOrmals_at_GaussPt(_nOrmals_at_GaussPt),
    tAngent1_at_GaussPt(_tAngent1_at_GaussPt),
    tAngent2_at_GaussPt(_tAngent2_at_GaussPt) {}

  ublas::matrix<FieldData> sPin;
  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data);

  PetscErrorCode calculateNormals();

};

/** \brief calculate normals at Gauss points of triangle element
  * \ingroup mofem_forces_and_sources
  */
struct OpGetNormalsOnPrism: public DataOperator {

  ublas::matrix<FieldData> &nOrmals_at_GaussPtF3;
  ublas::matrix<FieldData> &tAngent1_at_GaussPtF3;
  ublas::matrix<FieldData> &tAngent2_at_GaussPtF3;
  ublas::matrix<FieldData> &nOrmals_at_GaussPtF4;
  ublas::matrix<FieldData> &tAngent1_at_GaussPtF4;
  ublas::matrix<FieldData> &tAngent2_at_GaussPtF4;

  OpGetNormalsOnPrism(
    ublas::matrix<FieldData> &_nOrmals_at_GaussPtF3,
    ublas::matrix<FieldData> &_tAngent1_at_GaussPtF3,
    ublas::matrix<FieldData> &_tAngent2_at_GaussPtF3,
    ublas::matrix<FieldData> &_nOrmals_at_GaussPtF4,
    ublas::matrix<FieldData> &_tAngent1_at_GaussPtF4,
    ublas::matrix<FieldData> &_tAngent2_at_GaussPtF4):
    nOrmals_at_GaussPtF3(_nOrmals_at_GaussPtF3),
    tAngent1_at_GaussPtF3(_tAngent1_at_GaussPtF3),
    tAngent2_at_GaussPtF3(_tAngent2_at_GaussPtF3),
    nOrmals_at_GaussPtF4(_nOrmals_at_GaussPtF4),
    tAngent1_at_GaussPtF4(_tAngent1_at_GaussPtF4),
    tAngent2_at_GaussPtF4(_tAngent2_at_GaussPtF4) {}

  ublas::matrix<FieldData> sPin;
  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data);

  PetscErrorCode calculateNormals();

};


/** \brief transform Hdiv space fluxes from reference element to physical triangle
 */
struct OpSetPiolaTransoformOnTriangle: public DataOperator {

  const ublas::vector<double> &normal;
  const ublas::matrix<FieldData> &nOrmals_at_GaussPt;

  OpSetPiolaTransoformOnTriangle(
    const ublas::vector<double> &_normal,
    const ublas::matrix<FieldData> &_nOrmals_at_GaussPt):
    normal(_normal),nOrmals_at_GaussPt(_nOrmals_at_GaussPt) {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data);

};



}

#endif //__DATAOPERATORS_HPP

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources Forces and sources
 ******************************************************************************/
