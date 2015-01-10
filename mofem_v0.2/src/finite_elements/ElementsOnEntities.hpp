/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
* --------------------------------------------------------------
*
* DESCRIPTION: FIXME
*
* This is not exactly procedure for linear elatic dynamics, since jacobian is
* evaluated at every time step and snes procedure is involved. However it is
* implemented like that, to test methodology for general nonlinear problem.
*
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


#ifndef __ELEMENTSONENTITIES_HPP
#define __ELEMENTSONENTITIES_HPP

using namespace boost::numeric;

namespace MoFEM {

/** \brief Tet finite element  
 * \ingroup mofem_forces_and_sources_tet_element 
 *
 * User is implementing own operator at Guass piint level, by own object
 * derived from TetElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to vecUserOpN and
 * vecUserOpSymmNN. 
 *
 */
struct TetElementForcesAndSourcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore dataH1;
  DerivedDataForcesAndSurcesCore derivedDataH1;
  DataForcesAndSurcesCore dataL2;
  DerivedDataForcesAndSurcesCore derivedDataL2;
  DataForcesAndSurcesCore dataHdiv;
  DerivedDataForcesAndSurcesCore derivedDataHdiv;

  OpSetInvJacH1 opSetInvJacH1;
  OpSetPiolaTransform opPiolaTransform;
  OpSetInvJacHdiv opSetInvJacHdiv;

  string meshPositionsFieldName;
  ublas::matrix<FieldData> hoCoordsAtGaussPtsPts;
  ublas::matrix<FieldData> hoGaussPtsJac;
  ublas::matrix<FieldData> hoGaussPtsInvJac;
  ublas::vector<FieldData> hoGaussPtsDetJac;

  OpGetData opHOatGaussPoints; ///< higher order geometry data at Gauss pts
  OpSetHoInvJacH1 opSetHoInvJacH1;
  OpSetHoPiolaTransform opSetHoPiolaTransform;
  OpSetHoInvJacHdiv opSetHoInvJacHdiv;

  TetElementForcesAndSourcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),
    dataH1(MBTET),derivedDataH1(dataH1),
    dataL2(MBTET),derivedDataL2(dataL2),
    dataHdiv(MBTET),derivedDataHdiv(dataHdiv),
    opSetInvJacH1(invJac),
    opPiolaTransform(vOlume,Jac),opSetInvJacHdiv(invJac),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHOatGaussPoints(hoCoordsAtGaussPtsPts,hoGaussPtsJac,3,3),
    opSetHoInvJacH1(hoGaussPtsInvJac),
    opSetHoPiolaTransform(hoGaussPtsDetJac,hoGaussPtsJac),
    opSetHoInvJacHdiv(hoGaussPtsInvJac) {};
    
  virtual ~TetElementForcesAndSourcesCore() {}

  ErrorCode rval;
  PetscErrorCode ierr;
  double vOlume;
  ublas::vector<double> coords;

  ublas::matrix<double> Jac;;
  ublas::matrix<double> invJac;

  ublas::matrix<double> gaussPts;
  ublas::matrix<double> coordsAtGaussPts;

  /** \brief default oparator for TET element
    * \ingroup mofem_forces_and_sources_tet_element
    */
  struct UserDataOperator: public DataOperator {
    string row_field_name;
    string col_field_name;
    bool symm;
    UserDataOperator(
      const string &_field_name):
	row_field_name(_field_name),col_field_name(_field_name),symm(true),ptrFE(NULL) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
	row_field_name(_row_field_name),col_field_name(_col_field_name),symm(true),ptrFE(NULL) {};
    virtual ~UserDataOperator() {}
    inline double getVolume() { return ptrFE->vOlume; }
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }
    inline ublas::matrix<double>& getHoCoordsAtGaussPtsPts() { return ptrFE->hoCoordsAtGaussPtsPts; }
    inline ublas::matrix<double>& getHoGaussPtsInvJac() { return ptrFE->hoGaussPtsInvJac; }
    inline ublas::vector<double>& getHoGaussPtsDetJac() { return ptrFE->hoGaussPtsDetJac; }
    inline const FEMethod* getFEMethod() { return ptrFE; }
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };
    PetscErrorCode setPtrFE(TetElementForcesAndSourcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }

    //differential operators
    PetscErrorCode getDivergenceMatrixOperato_Hdiv(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data,
      int gg,ublas::vector<FieldData> &div);

    private:
    TetElementForcesAndSourcesCore *ptrFE; 

  };

  boost::ptr_vector<UserDataOperator> vecUserOpN; 
  boost::ptr_vector<UserDataOperator> vecUserOpNN;

  /** \brief Use to push back operator for right hand side
   * It can be ussed to calculate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpN; }

  /** \brief Use to push back operator for left hand side
   * It can be ussed to calculate matrices or other quantities on mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Lhs() { return vecUserOpNN; }


  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()();
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  
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

/** \brief transfrom Hdiv space fluxes from reference elemento to physical triangle
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

/** \brief Tri finite element  
 * \ingroup mofem_forces_and_sources_tri_element
 *
 * User is implementing own operator at Guass piint level, by own object
 * derived from TriElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to vecUserOpN and
 * vecUserOpSymmNN. 
 *
 */
struct TriElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  ErrorCode rval;
  PetscErrorCode ierr;
  double aRea;;
  ublas::vector<double> normal;
  ublas::vector<double> coords;
  ublas::matrix<double> gaussPts;
  ublas::matrix<double> coordsAtGaussPts;

  DataForcesAndSurcesCore dataH1;
  DerivedDataForcesAndSurcesCore derivedDataH1;
  DataForcesAndSurcesCore dataHdiv;
  DerivedDataForcesAndSurcesCore derivedDataHdiv;

  string meshPositionsFieldName;

  ublas::matrix<FieldData> nOrmals_at_GaussPt;
  ublas::matrix<FieldData> tAngent1_at_GaussPt;
  ublas::matrix<FieldData> tAngent2_at_GaussPt;
  OpGetNormals opHONormals;
  OpSetPiolaTransoformOnTriangle opSetPiolaTransoformOnTriangle;

  TriElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),
    dataH1(MBTRI),derivedDataH1(dataH1),
    dataHdiv(MBTRI),derivedDataHdiv(dataHdiv),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHONormals(nOrmals_at_GaussPt,tAngent1_at_GaussPt,tAngent2_at_GaussPt),
    opSetPiolaTransoformOnTriangle(normal,nOrmals_at_GaussPt) {};

  /** \brief default oparator for TRI element
    * \ingroup mofem_forces_and_sources_tri_element
    */
  struct UserDataOperator: public DataOperator {
    string row_field_name;
    string col_field_name;
    bool symm;
    UserDataOperator(
      const string &_field_name):
	row_field_name(_field_name),col_field_name(_field_name),symm(true),ptrFE(NULL) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
	row_field_name(_row_field_name),col_field_name(_col_field_name),symm(true),ptrFE(NULL) {};
    virtual ~UserDataOperator() {}
    inline double getArea() { return ptrFE->aRea; }

    /** \bried get triangle normal
     */
    inline ublas::vector<double>& getNormal() { return ptrFE->normal; }

    /** \bried get triangle coords
     */
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }

    /** \bried get triangle Gauss pts.
     */
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }

    /** \bried get coordinates at Gauss pts.
     */
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }

    /** \bried if higher order geometry return normals at Gauss pts.
     */
    inline ublas::matrix<FieldData>& getNormals_at_GaussPt() { return ptrFE->nOrmals_at_GaussPt; }

    /** \bried if higher order geometry return normals at Gauss pts.
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<ublas::matrix<double> > getNormals_at_GaussPt(const int gg) { 
      return ublas::matrix_row<ublas::matrix<double> >(ptrFE->nOrmals_at_GaussPt,gg); 
    }

    /** \bried if higher order geometry return tangent vetor to triangle at Gauss pts.
     */
    inline ublas::matrix<FieldData>& getTangent1_at_GaussPt() { return ptrFE->tAngent1_at_GaussPt; }

    /** \bried if higher order geometry return tangent vetor to triangle at Gauss pts.
     */
    inline ublas::matrix<FieldData>& getTangent2_at_GaussPt() { return ptrFE->tAngent2_at_GaussPt; }

    /** \bried return pointer to triangle finite element object 
     */
    inline const TriElementForcesAndSurcesCore* getTriElementForcesAndSurcesCore() { return ptrFE; }

    /** \bried return pointer to FEMthod object
     */
    inline const FEMethod* getFEMethod() { return ptrFE; }

    /** \bried return pointer to NumeredMoFEMFiniteElement 
     */
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };

    PetscErrorCode setPtrFE(TriElementForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    TriElementForcesAndSurcesCore *ptrFE; 
  };

  boost::ptr_vector<UserDataOperator> vecUserOpN; 
  boost::ptr_vector<UserDataOperator> vecUserOpSymmNN;

  /** \brief Use to push back operator for right hand side
   * It can be ussed to calculate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpN; }

  /** \brief Use to push back operator for left hand side
   * It can be ussed to calculate matrices or other quantities on mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Lhs() { return vecUserOpSymmNN; }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()();
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

/** \brief Edge finite element  
 * \ingroup mofem_forces_and_sources
 *
 * User is implementing own operator at Guass piint level, by own object
 * derived from EdgeElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to vecUserOpN and
 * vecUserOpSymmNN. 
 *
 */
struct EdgeElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derivedData;

  EdgeElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),data(MBEDGE),derivedData(data) {};

  ErrorCode rval;
  PetscErrorCode ierr;
  double lEngth;;
  ublas::vector<double> dIrection;
  ublas::vector<double> coords;
  ublas::matrix<double> gaussPts;
  ublas::matrix<double> coordsAtGaussPts;

  /** \brief default oparator for EDGE element
    */
  struct UserDataOperator: public DataOperator {
    string row_field_name;
    string col_field_name;
    UserDataOperator(
      const string &_field_name):
	row_field_name(_field_name),col_field_name(_field_name),ptrFE(NULL) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
	row_field_name(_row_field_name),col_field_name(_col_field_name),ptrFE(NULL) {};
    virtual ~UserDataOperator() {}
    inline double getLength() { return ptrFE->lEngth; }
    inline ublas::vector<double>& getDirection() { return ptrFE->dIrection; }
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }
    inline const FEMethod* getFEMethod() { return ptrFE; }
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };
    PetscErrorCode setPtrFE(EdgeElementForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    EdgeElementForcesAndSurcesCore *ptrFE; 
  };

  boost::ptr_vector<UserDataOperator> vecUserOpN; 
  boost::ptr_vector<UserDataOperator> vecUserOpSymmNN;

  /** \brief Use to push back operator for right hand side
   * It can be ussed to calculate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpN; }

  /** \brief Use to push back operator for left hand side
   * It can be ussed to calculate matrices or other quantities on mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Lhs() { return vecUserOpSymmNN; }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()();
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

/** \brief Vertex finite element  
 * \ingroup mofem_forces_and_sources
 *
 * User is implementing own operator at Guass piint level, by own object
 * derived from VertexElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to vecUserOpN and
 * vecUserOpSymmNN. 
 *
 */
struct VertexElementForcesAndSourcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derivedData;
  string meshPositionsFieldName;

  VertexElementForcesAndSourcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),data(MBVERTEX),derivedData(data) {};

  ErrorCode rval;
  PetscErrorCode ierr;
  ublas::vector<double> coords;

  /** \brief default oparator for VERTEX element
    */
  struct UserDataOperator: public DataOperator {
    string row_field_name;
    string col_field_name;
    UserDataOperator(
      const string &_field_name):
	row_field_name(_field_name),col_field_name(_field_name),ptrFE(NULL) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
	row_field_name(_row_field_name),col_field_name(_col_field_name),ptrFE(NULL) {};
    virtual ~UserDataOperator() {}
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline const FEMethod* getFEMethod() { return ptrFE; }
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };
    PetscErrorCode setPtrFE(VertexElementForcesAndSourcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    VertexElementForcesAndSourcesCore *ptrFE; 
  };

  boost::ptr_vector<UserDataOperator> vecUserOpN; 
  boost::ptr_vector<UserDataOperator> vecUserOpSymmNN;

  /** \brief Use to push back operator for right hand side
   * It can be ussed to calculate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpN; }

  /** \brief Use to push back operator for left hand side
   * It can be ussed to calculate matrices or other quantities on mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Lhs() { return vecUserOpSymmNN; }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()();
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

}

#endif //__ELEMENTSONENTITIES_HPP

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources Forces and sources
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_tet_element Tetrahedral Element 
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_tri_element Triangular Element 
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/




