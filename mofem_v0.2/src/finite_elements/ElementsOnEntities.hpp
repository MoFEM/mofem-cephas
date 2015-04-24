/** \file ElementsOnEntities.hpp

  \brief Implementation of elements on entities.

  Those element are inherited by user to implement specific implementation of
  particular problem.
  
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

/** \brief structure to get information form mofem into DataForcesAndSurcesCore
  * \ingroup mofem_forces_and_sources
  * 
  */
struct ForcesAndSurcesCore: public FEMethod {

  PetscErrorCode ierr;

  FieldInterface& mField;
  ForcesAndSurcesCore(FieldInterface& _mField): 
    mField(_mField) {};
  virtual ~ForcesAndSurcesCore() {}

  PetscErrorCode getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);
  PetscErrorCode getOrder(const EntityType type,const FieldSpace space,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);
  PetscErrorCode getOrder(const string &field_name,const EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

  PetscErrorCode getEdgesSense(DataForcesAndSurcesCore &data);
  PetscErrorCode getTrisSense(DataForcesAndSurcesCore &data);

  PetscErrorCode getEdgesOrder(DataForcesAndSurcesCore &data,const FieldSpace space);
  PetscErrorCode getTrisOrder(DataForcesAndSurcesCore &data,const FieldSpace space);
  PetscErrorCode getTetsOrder(DataForcesAndSurcesCore &data,const FieldSpace space);
  PetscErrorCode getEdgesOrder(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisOrder(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsOrder(DataForcesAndSurcesCore &data,const string &field_name);	

  // ** Indices **

  /// \brief get node indices
  PetscErrorCode getNodesIndices(const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices);

  /// \brief get indices by type (generic function)
  PetscErrorCode getTypeIndices(const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<int> &indices);

  /// \brief get indices by type (generic function)
  PetscErrorCode getTypeIndices(
    const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,
    EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

  /// \brief get row node indices from FENumeredDofMoFEMEntity_multiIndex
  PetscErrorCode getRowNodesIndices(DataForcesAndSurcesCore &data,const string &field_name); 

  /// \brief get col node indices from FENumeredDofMoFEMEntity_multiIndex
  PetscErrorCode getColNodesIndices(DataForcesAndSurcesCore &data,const string &field_name);

  /// \brief get Edges row indices from FENumeredDofMoFEMEntity_multiIndex
  PetscErrorCode getEdgesRowIndices(DataForcesAndSurcesCore &data,const string &field_name);

  /// \brief get Edges col indices from FENumeredDofMoFEMEntity_multiIndex
  PetscErrorCode getEdgesColIndices(DataForcesAndSurcesCore &data,const string &field_name);

  /// \brief get Tris row indices from FENumeredDofMoFEMEntity_multiIndex
  PetscErrorCode getTrisRowIndices(DataForcesAndSurcesCore &data,const string &field_name);

  /// \brief get Tris col indices from FENumeredDofMoFEMEntity_multiIndex
  PetscErrorCode getTrisColIndices(DataForcesAndSurcesCore &data,const string &field_name);

  /// \brief get Tets row indices from FENumeredDofMoFEMEntity_multiIndex
  PetscErrorCode getTetsRowIndices(DataForcesAndSurcesCore &data,const string &field_name);

  /// \brief get Tets col indices from FENumeredDofMoFEMEntity_multiIndex
  PetscErrorCode getTetsColIndices(DataForcesAndSurcesCore &data,const string &field_name);

  // ** Data **

  PetscErrorCode getNodesFieldData(const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<double> &nodes_data);
  PetscErrorCode getTypeFieldData(const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,ublas::vector<double> &ent_field_data);
  PetscErrorCode getTypeFieldData(const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

  // ** DoFS **

  PetscErrorCode getNodesFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgesFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsFieldData(DataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getNodesFieldDofs(const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<const FEDofMoFEMEntity*> &nodes_dofs);
  PetscErrorCode getTypeFieldDofs(const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<const FEDofMoFEMEntity*> &ent_field_dofs);
  PetscErrorCode getTypeFieldDofs(const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

  PetscErrorCode getNodesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getFaceNodes(DataForcesAndSurcesCore &data);
  PetscErrorCode getSpacesOnEntities(DataForcesAndSurcesCore &data);

  // ** Data form NumeredDofMoFEMEntity_multiIndex **

  /// \brief get indices of nodal indices which are declared for problem but not this particular element
  PetscErrorCode getProblemNodesIndices(const string &field_name,const NumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices) const;

  /// \brief get indices by type (generic function) which are declared for problem but not this particular element
  PetscErrorCode getProblemTypeIndices(
    const string &field_name,const NumeredDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,ublas::vector<int> &indices) const;

  PetscErrorCode getProblemNodesRowIndices(const string &field_name,ublas::vector<int> &nodes_indices) const;
  PetscErrorCode getProblemTypeRowIndices(const string &field_name,EntityType type,int side_number,ublas::vector<int> &indices) const;
  PetscErrorCode getProblemNodesColIndices(const string &field_name,ublas::vector<int> &nodes_indices) const;
  PetscErrorCode getProblemTypeColIndices(const string &field_name,EntityType type,int side_number,ublas::vector<int> &indices) const;

  /** \brief computes approximation functions for tetrahedral and H1 space
    */
  PetscErrorCode shapeTETFunctions_H1(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);

  /** \brief computes approximation functions for tetrahedral and L2 space
    */
  PetscErrorCode shapeTETFunctions_L2(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);


  ublas::matrix<ublas::matrix<double> > N_face_edge;
  ublas::vector<ublas::matrix<double> > N_face_bubble;
  ublas::vector<ublas::matrix<double> > N_volume_edge;
  ublas::vector<ublas::matrix<double> > N_volume_face;
  ublas::matrix<double> N_volume_bubble;

  ublas::matrix<ublas::matrix<double> > diffN_face_edge;
  ublas::vector<ublas::matrix<double> > diffN_face_bubble;
  ublas::vector<ublas::matrix<double> > diffN_volume_edge;
  ublas::vector<ublas::matrix<double> > diffN_volume_face;
  ublas::matrix<double> diffN_volume_bubble;


  /** \brief computes approximation functions for tetrahedral and H1 space
    */
  PetscErrorCode shapeTETFunctions_Hdiv(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);


  /** \brief computes approximation functions for triangle and H1 space
    */
  PetscErrorCode shapeTRIFunctions_H1(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const int G_DIM);


  /** \brief computes approximation functions for triangle and H1 space
    */
  PetscErrorCode shapeTRIFunctions_Hdiv(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const int G_DIM);

  /** \brief computes approximation functions for edge and H1 space
    */
  PetscErrorCode shapeEDGEFunctions_H1(
    DataForcesAndSurcesCore &data,const double *G_X,const int G_DIM);

  /** \brief computes approximation functions for prism and H1 space
    */
  PetscErrorCode shapeFlatPRISMFunctions_H1(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const int G_DIM);


  /** \brief computes approximation functions for prism and H1 space
    */
  PetscErrorCode shapeFlatPRISMFunctions_Hdiv(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const int G_DIM);

  /** \brief it is used to calculate nb. of Gauss integration points
   
   This function in general should be overload, returning integration rank
   depending on operator type. Integration rule should be set to
   integrate matrix and vector exactly.
   
   If function return -1
   \code
   int getRule(int order) { return -1; };
   \endcode
   then, fuction \codes setGaussPts(order) \endcode is called. In setGaussPts
   user can implement own intergartion rule for specific approx. ordrr.  
   
   At this stage of development integration points are weight are calculated following this paper:
   Albert Nijenhuis, Herbert Wilf, Combinatorial Algorithms for Computers and
   Calculators, Second Edition, Academic Press, 1978, ISBN: 0-12-519260-6,
   LC: QA164.N54.

   More details about algorithm
   https://github.com/johannesgerer/jburkardt-m/tree/master/gm_rule
   http://people.sc.fsu.edu/~jburkardt/cpp_src/gm_rule/gm_rule.html
  **/ 
  virtual int getRule(int order) { return order; };

  /** \brief set user specific integration rule
  
    This function allows for user defined integration rule. The key is to
    called matrix gaussPts, which is used by other MoFEM procedures. Matrix has
    number of rows equal to problem dimension plus one, where last index is used to
    store weight values. Number of columns is equal to number of integration points.

    Note: that matrix is called gussPts, however user can keep in it any integration rule. 

    User sets 
    \code 
    ublas::matrix<double> gaussPts;
    \endcode
    where 
    \code
    gaussPts.resize(dim+1,nb_gauss_pts);
    \endcode
    number rows represents local coordinates of integration points
    in reference element, where last index in row is for integration weight.

    */
  virtual PetscErrorCode setGaussPts(int order) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented");
    PetscFunctionReturn(0);
  }



  struct UserDataOperator: public DataOperator {

    string rowFieldName;
    string colFieldName;
    bool sYmm;

    UserDataOperator(
      const string &_field_name):
	rowFieldName(_field_name),colFieldName(_field_name),sYmm(true),ptrFE(NULL) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
	rowFieldName(_row_field_name),colFieldName(_col_field_name),sYmm(true),ptrFE(NULL) {};
    virtual ~UserDataOperator() {}

    /** \bried return pointer to NumeredMoFEMFiniteElement 
     */
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };

    /** \brief Get row indices 
    
    Field could be or not declared for this element but is declared for porblem
  
    \param field_name 
    \param type entity type
    \param side side number, any number if type is MBVERTEX
    \return indices
  
    NOTE: Using those indices to assemble matrix will result in error if new non-zero values need to be created.
  
    */
    PetscErrorCode getPorblemRowIndices(const string filed_name,const EntityType type,const int side,ublas::vector<int>& indices);
  
    /** \brief Get col indices 
    
    Field could be or not declared for this element but is declared for porblem
  
    \param field_name 
    \param type entity type
    \param side side number, any number if type is MBVERTEX
    \return indices
  
    NOTE: Using those indices to assemble matrix will result in error if new non-zero values need to be created.
  
    */
    PetscErrorCode getPorblemColIndices(const string filed_name,const EntityType type,const int side,ublas::vector<int>& indices);
  
    virtual PetscErrorCode setPtrFE(ForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    
    private:
    ForcesAndSurcesCore *ptrFE; 

  };

  boost::ptr_vector<UserDataOperator> rowOpPtrVector; 
  boost::ptr_vector<UserDataOperator> colOpPtrVector; 
  boost::ptr_vector<UserDataOperator> rowColOpPtrVector;

  /** \brief Use to push back operator for row operator

   It can be used to calculate nodal forces or other quantities on the mesh.

   */
  boost::ptr_vector<UserDataOperator>& getRowOpPtrVector() { return rowOpPtrVector; }

  /** \brief Use to push back operator for col operator

   It can be used to calculate nodal forces or other quantities on the mesh.

   */
  boost::ptr_vector<UserDataOperator>& getColOpPtrVector() { return colOpPtrVector; }


  /** \brief use to push back operator for row-col operator

   it can be used to calculate matrices or other quantities on mesh.

   */
  boost::ptr_vector<UserDataOperator>& getRowColOpPtrVector() { return rowColOpPtrVector; }

  DEPRECATED boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return getRowOpPtrVector(); }
  DEPRECATED boost::ptr_vector<UserDataOperator>& get_op_to_do_Lhs() { return getRowColOpPtrVector(); }

  virtual PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  virtual PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  virtual PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

/** \brief Volume finite element  
 \ingroup mofem_forces_and_sources_tet_element 
 
 User is implementing own operator at Gauss point level, by own object
 derived from VolumeElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector. 
 
 */
struct VolumeElementForcesAndSourcesCore: public ForcesAndSurcesCore {

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

  VolumeElementForcesAndSourcesCore(FieldInterface &_mField):
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
    
  virtual ~VolumeElementForcesAndSourcesCore() {}

  ErrorCode rval;
  double vOlume;
  ublas::vector<double> coords;

  ublas::matrix<double> Jac;;
  ublas::matrix<double> invJac;

  ublas::matrix<double> gaussPts;
  ublas::matrix<double> coordsAtGaussPts;

  /** \brief default operator for TET element
    * \ingroup mofem_forces_and_sources_tet_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &_field_name): ForcesAndSurcesCore::UserDataOperator(_field_name) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
      ForcesAndSurcesCore::UserDataOperator(_row_field_name,_col_field_name) {};
    
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


    /** \bried return pointer to Generic Tetrahedral Finite Element object
     */ 
    inline const VolumeElementForcesAndSourcesCore* getTetFE() { return ptrFE; }

    //differential operators
    PetscErrorCode getDivergenceMatrixOperato_Hdiv(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data,
      int gg,ublas::vector<double> &div);

    PetscErrorCode setPtrFE(ForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = dynamic_cast<VolumeElementForcesAndSourcesCore*>(ptr);
      ForcesAndSurcesCore::UserDataOperator::setPtrFE(ptr);
      PetscFunctionReturn(0);
    }

    private:
    VolumeElementForcesAndSourcesCore *ptrFE; 

  };


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

DEPRECATED typedef VolumeElementForcesAndSourcesCore TetElementForcesAndSourcesCore;


/** \brief Face finite element  
 \ingroup mofem_forces_and_sources_tri_element
 
 User is implementing own operator at Gauss point level, by own object
 derived from FaceElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector. 
 
 */
struct FaceElementForcesAndSourcesCore: public ForcesAndSurcesCore {

  ErrorCode rval;
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

  FaceElementForcesAndSourcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),
    dataH1(MBTRI),derivedDataH1(dataH1),
    dataHdiv(MBTRI),derivedDataHdiv(dataHdiv),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHONormals(nOrmals_at_GaussPt,tAngent1_at_GaussPt,tAngent2_at_GaussPt),
    opSetPiolaTransoformOnTriangle(normal,nOrmals_at_GaussPt) {};

  /** \brief default operator for TRI element
    * \ingroup mofem_forces_and_sources_tri_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &_field_name): ForcesAndSurcesCore::UserDataOperator(_field_name) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
      ForcesAndSurcesCore::UserDataOperator(_row_field_name,_col_field_name) {};

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
    inline const FaceElementForcesAndSourcesCore* getFaceElementForcesAndSourcesCore() { return ptrFE; }

    /** \bried return pointer to FEMthod object
     */
    inline const FEMethod* getFEMethod() { return ptrFE; }

     /** \bried return pointer to Generic Triangle Finite Element object
     */ 
    inline const FaceElementForcesAndSourcesCore* getTriFE() { return ptrFE; }

    PetscErrorCode setPtrFE(ForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = dynamic_cast<FaceElementForcesAndSourcesCore*>(ptr);
      ForcesAndSurcesCore::UserDataOperator::setPtrFE(ptr);
      PetscFunctionReturn(0);
    }

    private:
    FaceElementForcesAndSourcesCore *ptrFE; 


  };

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

DEPRECATED typedef FaceElementForcesAndSourcesCore TriElementForcesAndSurcesCore;

/** \brief Edge finite element  
 * \ingroup mofem_forces_and_sources
 *
 * User is implementing own operator at Gauss points level, by own object
 * derived from EdgeElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to rowOpPtrVector and
 * rowColOpPtrVector. 
 *
 */
struct EdgeElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore dataH1;
  DerivedDataForcesAndSurcesCore derivedDataH1;

  EdgeElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),dataH1(MBEDGE),derivedDataH1(dataH1) {};

  ErrorCode rval;
  double lEngth;;
  ublas::vector<double> dIrection;
  ublas::vector<double> coords;
  ublas::matrix<double> gaussPts;
  ublas::matrix<double> coordsAtGaussPts;

  /** \brief default operator for EDGE element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &_field_name): ForcesAndSurcesCore::UserDataOperator(_field_name) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
      ForcesAndSurcesCore::UserDataOperator(_row_field_name,_col_field_name) {};

    inline double getLength() { return ptrFE->lEngth; }
    inline ublas::vector<double>& getDirection() { return ptrFE->dIrection; }
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }
    inline const FEMethod* getFEMethod() { return ptrFE; }
    inline const EdgeElementForcesAndSurcesCore* getEdgeFE() { return ptrFE; }

    PetscErrorCode setPtrFE(ForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = dynamic_cast<EdgeElementForcesAndSurcesCore*>(ptr);
      ForcesAndSurcesCore::UserDataOperator::setPtrFE(ptr);
      PetscFunctionReturn(0);
    }

    private:
    EdgeElementForcesAndSurcesCore *ptrFE; 
  };

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

 User is implementing own operator at Gauss points level, by own object
 derived from VertexElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector. 

 \bug colOp not implemented
 
 */
struct VertexElementForcesAndSourcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derivedData;
  string meshPositionsFieldName;

  VertexElementForcesAndSourcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),data(MBVERTEX),derivedData(data) {};

  ErrorCode rval;
  ublas::vector<double> coords;

  /** \brief default operator for VERTEX element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &_field_name): ForcesAndSurcesCore::UserDataOperator(_field_name) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
      ForcesAndSurcesCore::UserDataOperator(_row_field_name,_col_field_name) {};

    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline const FEMethod* getFEMethod() { return ptrFE; }

    PetscErrorCode setPtrFE(ForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = dynamic_cast<VertexElementForcesAndSourcesCore*>(ptr);
      ForcesAndSurcesCore::UserDataOperator::setPtrFE(ptr);
      PetscFunctionReturn(0);
    }

    private:
    VertexElementForcesAndSourcesCore *ptrFE; 
  };

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
 \ingroup mofem_forces_and_sources_prism_element
 
 User is implementing own operator at Gauss points level, by own object
 derived from FlatPrismElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector.

 \bug colOp not implemented
 
 */
struct FlatPrismElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  ErrorCode rval;
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
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &_field_name): ForcesAndSurcesCore::UserDataOperator(_field_name) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
      ForcesAndSurcesCore::UserDataOperator(_row_field_name,_col_field_name) {};

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

    PetscErrorCode setPtrFE(ForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = dynamic_cast<FlatPrismElementForcesAndSurcesCore*>(ptr);
      ForcesAndSurcesCore::UserDataOperator::setPtrFE(ptr);
      PetscFunctionReturn(0);
    }

    private:
    FlatPrismElementForcesAndSurcesCore *ptrFE; 

  };

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




