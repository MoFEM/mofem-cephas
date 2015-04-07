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
 * User is implementing own operator at Gauss point level, by own object
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
  ublas::matrix<double> hoCoordsAtGaussPts;
  ublas::matrix<double> hoGaussPtsJac;
  ublas::matrix<double> hoGaussPtsInvJac;
  ublas::vector<double> hoGaussPtsDetJac;

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
    opHOatGaussPoints(hoCoordsAtGaussPts,hoGaussPtsJac,3,3),
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

  /** \brief default operator for TET element
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
    
    /** \brief element volume (linear geometry)
      */
    inline double getVolume() { return ptrFE->vOlume; }

    /** \brief nodal coordinates 
      */ 
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }

    /** \brief matrix of Gauss pts
      */
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }

    /** \brief Gauss points and weight, matrix (nb. of points x 4)

      Column 0-3 and 4 represents Gauss pts coordinate and weight, respectively.

      */
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }

    /** \brief coordinate at Gauss points (if hierarchical approximation of element geometry)
      */
    inline ublas::matrix<double>& getHoCoordsAtGaussPts() { return ptrFE->hoCoordsAtGaussPts; }

    inline ublas::matrix<double>& getHoGaussPtsInvJac() { return ptrFE->hoGaussPtsInvJac; }
    inline ublas::vector<double>& getHoGaussPtsDetJac() { return ptrFE->hoGaussPtsDetJac; }

    /** \bried return pointer to Generic Finite Element object
     */ 
    inline const FEMethod* getFEMethod() { return ptrFE; }

    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };

    /** \bried return pointer to Generic Tetrahedral Finite Element object
     */ 
    inline const TetElementForcesAndSourcesCore* getTetFE() { return ptrFE; }

    PetscErrorCode setPtrFE(TetElementForcesAndSourcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }

    //differential operators
    PetscErrorCode getDivergenceMatrixOperato_Hdiv(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data,
      int gg,ublas::vector<double> &div);

    private:
    TetElementForcesAndSourcesCore *ptrFE; 

  };

  boost::ptr_vector<UserDataOperator> vecUserOpN; 
  boost::ptr_vector<UserDataOperator> vecUserOpNN;

  /** \brief Use to push back operator for right hand side
   * It can be used to calculate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpN; }

  /** \brief Use to push back operator for left hand side
   * It can be used to calculate matrices or other quantities on mesh.
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


/** \brief Tri finite element  
 * \ingroup mofem_forces_and_sources_tri_element
 *
 * User is implementing own operator at Gauss point level, by own object
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

  ublas::matrix<double> nOrmals_at_GaussPt;
  ublas::matrix<double> tAngent1_at_GaussPt;
  ublas::matrix<double> tAngent2_at_GaussPt;
  OpGetNormals opHONormals;
  OpSetPiolaTransoformOnTriangle opSetPiolaTransoformOnTriangle;

  TriElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),
    dataH1(MBTRI),derivedDataH1(dataH1),
    dataHdiv(MBTRI),derivedDataHdiv(dataHdiv),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHONormals(nOrmals_at_GaussPt,tAngent1_at_GaussPt,tAngent2_at_GaussPt),
    opSetPiolaTransoformOnTriangle(normal,nOrmals_at_GaussPt) {};

  /** \brief default operator for TRI element
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

    /** \bried get triangle coordinates
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
    inline ublas::matrix<double>& getNormals_at_GaussPt() { return ptrFE->nOrmals_at_GaussPt; }

    /** \bried if higher order geometry return normals at Gauss pts.
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<ublas::matrix<double> > getNormals_at_GaussPt(const int gg) { 
      return ublas::matrix_row<ublas::matrix<double> >(ptrFE->nOrmals_at_GaussPt,gg); 
    }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline ublas::matrix<double>& getTangent1_at_GaussPt() { return ptrFE->tAngent1_at_GaussPt; }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline ublas::matrix<double>& getTangent2_at_GaussPt() { return ptrFE->tAngent2_at_GaussPt; }

    /** \bried return pointer to triangle finite element object 
     */
    inline const TriElementForcesAndSurcesCore* getTriElementForcesAndSurcesCore() { return ptrFE; }

    /** \bried return pointer to FEMthod object
     */
    inline const FEMethod* getFEMethod() { return ptrFE; }

     /** \bried return pointer to Generic Triangle Finite Element object
     */ 
    inline const TriElementForcesAndSurcesCore* getTriFE() { return ptrFE; }

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
   * It can be used to calculate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpN; }

  /** \brief Use to push back operator for left hand side
   * It can be used to calculate matrices or other quantities on mesh.
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
 * User is implementing own operator at Gauss points level, by own object
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

  /** \brief default operator for EDGE element
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
    inline const EdgeElementForcesAndSurcesCore* getEdgeFE() { return ptrFE; }
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
   * It can be used to calculate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpN; }

  /** \brief Use to push back operator for left hand side
   * It can be used to calculate matrices or other quantities on mesh.
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
 * User is implementing own operator at Gauss points level, by own object
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

  /** \brief default operator for VERTEX element
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
   * It can be used to calculate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpN; }

  /** \brief Use to push back operator for left hand side
   * It can be used to calculate matrices or other quantities on mesh.
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


/** \brief FlatPrism finite element  
 * \ingroup mofem_forces_and_sources_prism_element
 *
 * User is implementing own operator at Gauss points level, by own object
 * derived from FlatPrismElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to vecUserOpN and
 * vecUserOpSymmNN. 
 *
 */
struct FlatPrismElementForcesAndSurcesCore: public ForcesAndSurcesCore {

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

  ublas::matrix<double> nOrmals_at_GaussPtF3;
  ublas::matrix<double> tAngent1_at_GaussPtF3;
  ublas::matrix<double> tAngent2_at_GaussPtF3;
  ublas::matrix<double> nOrmals_at_GaussPtF4;
  ublas::matrix<double> tAngent1_at_GaussPtF4;
  ublas::matrix<double> tAngent2_at_GaussPtF4;
  OpGetNormalsOnPrism opHONormals;

  FlatPrismElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),
    dataH1(MBPRISM),derivedDataH1(dataH1),
    dataHdiv(MBPRISM),derivedDataHdiv(dataHdiv),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHONormals(
    nOrmals_at_GaussPtF3,tAngent1_at_GaussPtF3,tAngent2_at_GaussPtF3,
    nOrmals_at_GaussPtF4,tAngent1_at_GaussPtF4,tAngent2_at_GaussPtF4) {};

  /** \brief default operator for TRI element
    * \ingroup mofem_forces_and_sources_prism_element
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

    /** \bried get triangle coordinates
     */
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }

    /** \bried get triangle Gauss pts.
     */
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }

    /** \bried get coordinates at Gauss pts.
     */
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }

    /** \bried if higher order geometry return normals at face F3 at Gauss pts.
     * 
     * Face 3 is top face in canonical triangle numeration, see 
     * Canonical numbering systems for finite-element codes Timothy J. Tautges
     */
    inline ublas::matrix<double>& getNormals_at_GaussPtF3() { return ptrFE->nOrmals_at_GaussPtF3; }

    /** \bried if higher order geometry return normals at face F4 at Gauss pts.
     * 
     * Face 4 is top face in canonical triangle numeration, see 
     * Canonical numbering systems for finite-element codes Timothy J. Tautges
     */
    inline ublas::matrix<double>& getNormals_at_GaussPtF4() { return ptrFE->nOrmals_at_GaussPtF4; }

    /** \bried if higher order geometry return normals at Gauss pts.
      *
      * Face 3 is top face in canonical triangle numeration, see 
      * Canonical numbering systems for finite-element codes Timothy J. Tautges
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<ublas::matrix<double> > getNormals_at_GaussPtF3(const int gg) { 
      return ublas::matrix_row<ublas::matrix<double> >(ptrFE->nOrmals_at_GaussPtF3,gg); 
    }

    /** \bried if higher order geometry return normals at Gauss pts.
      *
      * Face 3 is top face in canonical triangle numeration, see 
      * Canonical numbering systems for finite-element codes Timothy J. Tautges
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<ublas::matrix<double> > getNormals_at_GaussPtF4(const int gg) { 
      return ublas::matrix_row<ublas::matrix<double> >(ptrFE->nOrmals_at_GaussPtF4,gg); 
    }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline ublas::matrix<double>& getTangent1_at_GaussPtF3() { return ptrFE->tAngent1_at_GaussPtF3; }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline ublas::matrix<double>& getTangent2_at_GaussPtF3() { return ptrFE->tAngent2_at_GaussPtF3; }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline ublas::matrix<double>& getTangent1_at_GaussPtF4() { return ptrFE->tAngent1_at_GaussPtF4; }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline ublas::matrix<double>& getTangent2_at_GaussPtF4() { return ptrFE->tAngent2_at_GaussPtF4; }

    /** \bried return pointer to triangle finite element object 
     */
    inline const FlatPrismElementForcesAndSurcesCore* getFlatPrismElementForcesAndSurcesCore() { return ptrFE; }

    /** \bried return pointer to FEMthod object
     */
    inline const FEMethod* getFEMethod() { return ptrFE; }

    /** \bried return pointer to NumeredMoFEMFiniteElement 
     */
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };

    PetscErrorCode setPtrFE(FlatPrismElementForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    FlatPrismElementForcesAndSurcesCore *ptrFE; 
  };

  boost::ptr_vector<UserDataOperator> vecUserOpN; 
  boost::ptr_vector<UserDataOperator> vecUserOpSymmNN;

  /** \brief Use to push back operator for right hand side
   * It can be used to calculate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpN; }

  /** \brief Use to push back operator for left hand side
   * It can be used to calculate matrices or other quantities on mesh.
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

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_prism_element Prism Element 
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/




