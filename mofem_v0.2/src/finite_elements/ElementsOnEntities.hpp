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
    FENumeredDofMoFEMEntity_multiIndex &dofs,
    VectorInt &nodes_indices
  );

  /// \brief get indices by type (generic function)
  PetscErrorCode getTypeIndices(const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,
    VectorInt &indices
  );

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

  /// \brief get NoField indices
  PetscErrorCode getNoFieldIndices(
    const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,VectorInt &nodes_indices
  );

  /// \brief get col NoField indices
  PetscErrorCode getNoFieldRowIndices(DataForcesAndSurcesCore &data,const string &field_name);

  /// \brief get col NoField indices
  PetscErrorCode getNoFieldColIndices(DataForcesAndSurcesCore &data,const string &field_name);

  // ** Data **

  PetscErrorCode getNodesFieldData(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,VectorDouble &nodes_data
  );

  PetscErrorCode getTypeFieldData(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,VectorDouble &ent_field_data
  );

  PetscErrorCode getTypeFieldData(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
  );

  PetscErrorCode getNoFieldFieldData(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,VectorDouble &ent_field_data
  );

  PetscErrorCode getNoFieldFieldData(
    DataForcesAndSurcesCore &data,const string &field_name
  );

  PetscErrorCode getNodesFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgesFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsFieldData(DataForcesAndSurcesCore &data,const string &field_name);

  // ** DoFS **

  PetscErrorCode getNodesFieldDofs(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,VectorDofs &nodes_dofs
  );

  PetscErrorCode getTypeFieldDofs(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,VectorDofs &ent_field_dofs
  );

  PetscErrorCode getTypeFieldDofs(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
  );

  PetscErrorCode getNoFieldFieldDofs(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,VectorDofs &nodes_dofs
  );

  PetscErrorCode getNoFieldFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getNodesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getFaceNodes(DataForcesAndSurcesCore &data);
  PetscErrorCode getSpacesOnEntities(DataForcesAndSurcesCore &data);

  // ** Data form NumeredDofMoFEMEntity_multiIndex **

  /// \brief get indices of nodal indices which are declared for problem but not this particular element
  PetscErrorCode getProblemNodesIndices(const string &field_name,const NumeredDofMoFEMEntity_multiIndex &dofs,VectorInt &nodes_indices) const;

  /// \brief get indices by type (generic function) which are declared for problem but not this particular element
  PetscErrorCode getProblemTypeIndices(
    const string &field_name,const NumeredDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,VectorInt &indices) const;

  PetscErrorCode getProblemNodesRowIndices(const string &field_name,VectorInt &nodes_indices) const;
  PetscErrorCode getProblemTypeRowIndices(const string &field_name,EntityType type,int side_number,VectorInt &indices) const;
  PetscErrorCode getProblemNodesColIndices(const string &field_name,VectorInt &nodes_indices) const;
  PetscErrorCode getProblemTypeColIndices(const string &field_name,EntityType type,int side_number,VectorInt &indices) const;

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

  ublas::matrix<MatrixDouble > N_face_edge;
  ublas::vector<MatrixDouble > N_face_bubble;
  ublas::vector<MatrixDouble > N_volume_edge;
  ublas::vector<MatrixDouble > N_volume_face;
  MatrixDouble N_volume_bubble;

  ublas::matrix<MatrixDouble > diffN_face_edge;
  ublas::vector<MatrixDouble > diffN_face_bubble;
  ublas::vector<MatrixDouble > diffN_volume_edge;
  ublas::vector<MatrixDouble > diffN_volume_face;
  MatrixDouble diffN_volume_bubble;


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
    MatrixDouble gaussPts;
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

  /** \brief Data operator to do calculations at integration points.
    * \ingroup mofem_forces_and_sources

    Is inherited and implemented by user to do calculations. It can be used in many
    different ways but typically is used to integrate matrices (f.e. stiffness matrux) and
    the right hand vector (f.e. force vector).

    Note: It is assumed that operator is executed for symmetric problem. That means that
    is executed for not repeating entities on finite element. For example on triangle we have
    nodes, 3 edges and one face. Because of symmetry calculations are for:
    nodes-nodes,
    nodes-edge0, nodes-edge_1, nodes-edge_2,
    nodes-face,
    edges_1-edges_1, edges_1-edge_1, edge_1-edge_2,
    edge_1-edge_1, edge_1-edge_2,
    edge_2-edge_2,
    edge_1-face, edge_1-face, edges_3-face,
    face - face

    In case of non symmetric problem in addition calculations of lower off diagonal terms. F.e.
    edge_1-edge_0,
    esges_3-edge_0, edge_3-edge_1,

    In that case class variable UserDataOperator::sYmm = false;

    NoteL: By default sYmm is set for symmetric problems

  */
  struct UserDataOperator: public DataOperator {

    string rowFieldName;
    string colFieldName;
    bool sYmm;

    /// set if operator is executed taking in account symmetry
    inline void setSymm() { sYmm = true; }

    /// unset if operator is executed for  non symmetric problem
    inline void unSetSymm() { sYmm = false; }

    enum OpType {
      OPROW = 1<<0,
      OPCOL = 1<<1,
      OPROWCOL = 1<<2
    };
    char opType;

    inline int getOpType() const { return opType; }
    inline void setOpType(const OpType type) { opType = type; }
    inline void addOpType(const OpType type) { opType |= type; }

    UserDataOperator(const string &_field_name,const char type):
      rowFieldName(_field_name),
      colFieldName(_field_name),
      sYmm(true),
      opType(type),
      ptrFE(NULL) {};

    UserDataOperator(
        const string &_row_field_name,const string &_col_field_name,const char type
      ):
      rowFieldName(_row_field_name),
      colFieldName(_col_field_name),
      sYmm(true),
      opType(type),
      ptrFE(NULL) {}
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
    PetscErrorCode getPorblemRowIndices(const string filed_name,const EntityType type,const int side,VectorInt& indices);

    /** \brief Get col indices

    Field could be or not declared for this element but is declared for porblem

    \param field_name
    \param type entity type
    \param side side number, any number if type is MBVERTEX
    \return indices

    NOTE: Using those indices to assemble matrix will result in error if new non-zero values need to be created.

    */
    PetscErrorCode getPorblemColIndices(const string filed_name,const EntityType type,const int side,VectorInt& indices);

    virtual PetscErrorCode setPtrFE(ForcesAndSurcesCore *ptr) {
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }

    /** \bried return pointer to Generic Finite Element object
     */
    inline const FEMethod* getFEMethod() { return ptrFE; }

    private:
    ForcesAndSurcesCore *ptrFE;

  };

  boost::ptr_vector<UserDataOperator> opPtrVector;

  /** \brief Use to push back operator for row operator

   It can be used to calculate nodal forces or other quantities on the mesh.

   */
  boost::ptr_vector<UserDataOperator>& getOpPtrVector() { return opPtrVector; }

  /** \brief Use to push back operator for row operator

   It can be used to calculate nodal forces or other quantities on the mesh.

   This function is DEPRECATED: use getOpPtrVector instead.

   */
  DEPRECATED boost::ptr_vector<UserDataOperator>& getRowOpPtrVector() { return opPtrVector; }

  /** \brief Use to push back operator for col operator

   It can be used to calculate nodal forces or other quantities on the mesh.

   This function is DEPRECATED: use getOpPtrVector instead.

   */
  DEPRECATED boost::ptr_vector<UserDataOperator>& getColOpPtrVector() { return opPtrVector; }


  /** \brief use to push back operator for row-col operator

   It can be used to calculate matrices or other quantities on mesh.

   This function is DEPRECATED: use getOpPtrVector instead.

   */
  DEPRECATED boost::ptr_vector<UserDataOperator>& getRowColOpPtrVector() { return opPtrVector; }

  DEPRECATED boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return getOpPtrVector(); }
  DEPRECATED boost::ptr_vector<UserDataOperator>& get_op_to_do_Lhs() { return getOpPtrVector(); }

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
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;

  OpSetInvJacH1 opSetInvJacH1;
  OpSetPiolaTransform opPiolaTransform;
  OpSetInvJacHdiv opSetInvJacHdiv;

  string meshPositionsFieldName;
  MatrixDouble hoCoordsAtGaussPts;
  MatrixDouble hoGaussPtsJac;
  MatrixDouble hoGaussPtsInvJac;
  VectorDouble hoGaussPtsDetJac;

  OpGetDataAndGradient opHOatGaussPoints; ///< higher order geometry data at Gauss pts
  OpSetHoInvJacH1 opSetHoInvJacH1;
  OpSetHoPiolaTransform opSetHoPiolaTransform;
  OpSetHoInvJacHdiv opSetHoInvJacHdiv;

  VolumeElementForcesAndSourcesCore(FieldInterface &m_field):
    ForcesAndSurcesCore(m_field),
    dataH1(MBTET),derivedDataH1(dataH1),
    dataL2(MBTET),derivedDataL2(dataL2),
    dataHdiv(MBTET),derivedDataHdiv(dataHdiv),
    dataNoField(MBTET),dataNoFieldCol(MBTET),
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
  VectorDouble coords;

  MatrixDouble Jac;;
  MatrixDouble invJac;

  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  /** \brief default operator for TET element
    * \ingroup mofem_forces_and_sources_tet_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(field_name,type) {}


    UserDataOperator(
      const string &row_field_name,const string &col_field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {};

    /** \brief element volume (linear geometry)
      */
    inline double getVolume() { return ptrFE->vOlume; }

    /** \brief nodal coordinates
      */
    inline VectorDouble& getCoords() { return ptrFE->coords; }

    /** \brief matrix of Gauss pts
      */
    inline MatrixDouble& getGaussPts() { return ptrFE->gaussPts; }

    /** \brief Gauss points and weight, matrix (nb. of points x 4)

      Column 0-3 and 4 represents Gauss pts coordinate and weight, respectively.

      */
    inline MatrixDouble& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }

    /** \brief coordinate at Gauss points (if hierarchical approximation of element geometry)
      */
    inline MatrixDouble& getHoCoordsAtGaussPts() { return ptrFE->hoCoordsAtGaussPts; }

    inline MatrixDouble& getHoGaussPtsInvJac() { return ptrFE->hoGaussPtsInvJac; }
    inline VectorDouble& getHoGaussPtsDetJac() { return ptrFE->hoGaussPtsDetJac; }

    /** \bried return pointer to Generic Tetrahedral Finite Element object
     */
    inline const VolumeElementForcesAndSourcesCore* getTetFE() { return ptrFE; }

    //differential operators
    PetscErrorCode getDivergenceMatrixOperato_Hdiv(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data,
      int gg,VectorDouble &div);

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
  VectorDouble normal;
  VectorDouble coords;
  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  DataForcesAndSurcesCore dataH1;
  DerivedDataForcesAndSurcesCore derivedDataH1;
  DataForcesAndSurcesCore dataHdiv;
  DerivedDataForcesAndSurcesCore derivedDataHdiv;
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;

  string meshPositionsFieldName;

  MatrixDouble nOrmals_at_GaussPt;
  MatrixDouble tAngent1_at_GaussPt;
  MatrixDouble tAngent2_at_GaussPt;
  OpGetNormals opHONormals;
  OpSetPiolaTransoformOnTriangle opSetPiolaTransoformOnTriangle;

  FaceElementForcesAndSourcesCore(FieldInterface &m_field):
    ForcesAndSurcesCore(m_field),
    dataH1(MBTRI),derivedDataH1(dataH1),
    dataHdiv(MBTRI),derivedDataHdiv(dataHdiv),
    dataNoField(MBTRI),dataNoFieldCol(MBTRI),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHONormals(nOrmals_at_GaussPt,tAngent1_at_GaussPt,tAngent2_at_GaussPt),
    opSetPiolaTransoformOnTriangle(normal,nOrmals_at_GaussPt) {};

  /** \brief default operator for TRI element
    * \ingroup mofem_forces_and_sources_tri_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(field_name,type) {}

    UserDataOperator(
      const string &row_field_name,const string &col_field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {};

    inline double getArea() { return ptrFE->aRea; }

    /** \bried get triangle normal
     */
    inline VectorDouble& getNormal() { return ptrFE->normal; }

    /** \bried get triangle coordinates
     */
    inline VectorDouble& getCoords() { return ptrFE->coords; }

    /** \bried get triangle Gauss pts.
     */
    inline MatrixDouble& getGaussPts() { return ptrFE->gaussPts; }

    /** \bried get coordinates at Gauss pts.
     */
    inline MatrixDouble& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }

    /** \bried if higher order geometry return normals at Gauss pts.
     */
    inline MatrixDouble& getNormals_at_GaussPt() { return ptrFE->nOrmals_at_GaussPt; }

    /** \bried if higher order geometry return normals at Gauss pts.
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<MatrixDouble > getNormals_at_GaussPt(const int gg) {
      return ublas::matrix_row<MatrixDouble >(ptrFE->nOrmals_at_GaussPt,gg);
    }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent1_at_GaussPt() { return ptrFE->tAngent1_at_GaussPt; }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent2_at_GaussPt() { return ptrFE->tAngent2_at_GaussPt; }

    /** \bried return pointer to triangle finite element object
     */
    inline const FaceElementForcesAndSourcesCore* getFaceElementForcesAndSourcesCore() { return ptrFE; }

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
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;
  string meshPositionsFieldName;

  MatrixDouble tAngent_at_GaussPt;
  OpGetHoTangentOnEdge opGetHoTangentOnEdge;

  EdgeElementForcesAndSurcesCore(FieldInterface &m_field):
    ForcesAndSurcesCore(m_field),
    dataH1(MBEDGE),
    derivedDataH1(dataH1),
    dataNoField(MBEDGE),
    dataNoFieldCol(MBEDGE),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opGetHoTangentOnEdge(tAngent_at_GaussPt)
  {};

  ErrorCode rval;
  double lEngth;;
  VectorDouble dIrection;
  VectorDouble cOords;
  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  /** \brief default operator for EDGE element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(field_name,type) {}

    UserDataOperator(
      const string &row_field_name,const string &col_field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {}

    inline double getLength() { return ptrFE->lEngth; }
    inline VectorDouble& getDirection() { return ptrFE->dIrection; }
    inline VectorDouble& getCoords() { return ptrFE->cOords; }
    inline MatrixDouble& getGaussPts() { return ptrFE->gaussPts; }
    inline MatrixDouble& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }
    inline MatrixDouble& getTangetAtGaussPtrs() { return ptrFE->tAngent_at_GaussPt; }
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

 */
struct VertexElementForcesAndSourcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derivedData;
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;
  string meshPositionsFieldName;

  VertexElementForcesAndSourcesCore(FieldInterface &m_field):
    ForcesAndSurcesCore(m_field),
    data(MBVERTEX),
    derivedData(data),
    dataNoField(MBVERTEX),
    dataNoFieldCol(MBVERTEX)
  {};

  ErrorCode rval;
  VectorDouble coords;

  /** \brief default operator for VERTEX element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(field_name,type) {}

    UserDataOperator(
      const string &row_field_name,const string &col_field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {}

    inline VectorDouble& getCoords() { return ptrFE->coords; }

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

 */
struct FlatPrismElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  ErrorCode rval;
  double aRea;;
  VectorDouble normal;
  VectorDouble coords;
  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  DataForcesAndSurcesCore dataH1;
  DerivedDataForcesAndSurcesCore derivedDataH1;
  DataForcesAndSurcesCore dataHdiv;
  DerivedDataForcesAndSurcesCore derivedDataHdiv;
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;

  string meshPositionsFieldName;

  MatrixDouble nOrmals_at_GaussPtF3;
  MatrixDouble tAngent1_at_GaussPtF3;
  MatrixDouble tAngent2_at_GaussPtF3;
  MatrixDouble nOrmals_at_GaussPtF4;
  MatrixDouble tAngent1_at_GaussPtF4;
  MatrixDouble tAngent2_at_GaussPtF4;
  OpGetNormalsOnPrism opHONormals;

  FlatPrismElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),
    dataH1(MBPRISM),derivedDataH1(dataH1),
    dataHdiv(MBPRISM),derivedDataHdiv(dataHdiv),
    dataNoField(MBPRISM),dataNoFieldCol(MBPRISM),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHONormals(
    nOrmals_at_GaussPtF3,tAngent1_at_GaussPtF3,tAngent2_at_GaussPtF3,
    nOrmals_at_GaussPtF4,tAngent1_at_GaussPtF4,tAngent2_at_GaussPtF4) {};

  /** \brief default operator for TRI element
    * \ingroup mofem_forces_and_sources_prism_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const string &field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(field_name,type) {}

    UserDataOperator(
      const string &row_field_name,const string &col_field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {}

    inline double getArea() { return ptrFE->aRea; }

    /** \bried get triangle normal
     */
    inline VectorDouble& getNormal() { return ptrFE->normal; }

    /** \bried get triangle coordinates
     */
    inline VectorDouble& getCoords() { return ptrFE->coords; }

    /** \bried get triangle Gauss pts.
     */
    inline MatrixDouble& getGaussPts() { return ptrFE->gaussPts; }

    /** \bried get coordinates at Gauss pts.
     */
    inline MatrixDouble& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }

    /** \bried if higher order geometry return normals at face F3 at Gauss pts.
     *
     * Face 3 is top face in canonical triangle numeration, see
     * Canonical numbering systems for finite-element codes Timothy J. Tautges
     */
    inline MatrixDouble& getNormals_at_GaussPtF3() { return ptrFE->nOrmals_at_GaussPtF3; }

    /** \bried if higher order geometry return normals at face F4 at Gauss pts.
     *
     * Face 4 is top face in canonical triangle numeration, see
     * Canonical numbering systems for finite-element codes Timothy J. Tautges
     */
    inline MatrixDouble& getNormals_at_GaussPtF4() { return ptrFE->nOrmals_at_GaussPtF4; }

    /** \bried if higher order geometry return normals at Gauss pts.
      *
      * Face 3 is top face in canonical triangle numeration, see
      * Canonical numbering systems for finite-element codes Timothy J. Tautges
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<MatrixDouble > getNormals_at_GaussPtF3(const int gg) {
      return ublas::matrix_row<MatrixDouble >(ptrFE->nOrmals_at_GaussPtF3,gg);
    }

    /** \bried if higher order geometry return normals at Gauss pts.
      *
      * Face 3 is top face in canonical triangle numeration, see
      * Canonical numbering systems for finite-element codes Timothy J. Tautges
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<MatrixDouble > getNormals_at_GaussPtF4(const int gg) {
      return ublas::matrix_row<MatrixDouble >(ptrFE->nOrmals_at_GaussPtF4,gg);
    }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent1_at_GaussPtF3() { return ptrFE->tAngent1_at_GaussPtF3; }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent2_at_GaussPtF3() { return ptrFE->tAngent2_at_GaussPtF3; }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent1_at_GaussPtF4() { return ptrFE->tAngent1_at_GaussPtF4; }

    /** \bried if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent2_at_GaussPtF4() { return ptrFE->tAngent2_at_GaussPtF4; }

    /** \bried return pointer to triangle finite element object
     */
    inline const FlatPrismElementForcesAndSurcesCore* getFlatPrismElementForcesAndSurcesCore() { return ptrFE; }

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
