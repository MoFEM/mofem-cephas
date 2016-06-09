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
  ForcesAndSurcesCore(FieldInterface& m_field):
    mField(m_field) {};
  virtual ~ForcesAndSurcesCore() {}

  PetscErrorCode getNumberOfNodes(int &num_nodes) const;

  /** \brief Get max order of approximation for data fields

  Method  getMaxDataOrder () return maximal order on entities, for
  all data on the element. So for example if finite element is triangle, and
  triangle base function have order 4 and on edges base function have order 2,
  this function return 4.

  If finite element has for example 2 or more approximated fields, for example
  Pressure (order 3) and displacement field (order 5), this function returns 5.

  */
  int getMaxDataOrder() const;

  /// \brief Get max order of approximation for field in rows
  int getMaxRowOrder() const;

  /// \brief Get max order of approximation for field in columns
  int getMaxColOrder() const;

  /**
   * \brief get sense (orientation) of entity
   * @param  type type of entity
   * @param  data entity data
   * @return      error code
   */
  PetscErrorCode getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) const;

  /// \brief get maximal approximation order of approximation on the entity
  PetscErrorCode getDataOrder(
    const EntityType type,const FieldSpace space,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
  ) const;

  /// \brief get maximal approximation order on entity
  PetscErrorCode getDataOrderSpaceAndBase(
    const std::string &field_name,const EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
  ) const;

  /**
   * get  edge sense )orientation) in respect finite element entity
   * @param  data structure where results are stored
   * @return      error code
   */
  PetscErrorCode getEdgesSense(DataForcesAndSurcesCore &data) const;

  /**
   * get triangle sense (orientation)
   * @param  data structure where results are stored
   * @return      error code
   */
  PetscErrorCode getTrisSense(DataForcesAndSurcesCore &data) const;
  PetscErrorCode getQuadSense(DataForcesAndSurcesCore &data) const;

  PetscErrorCode getEdgesDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const;
  PetscErrorCode getTrisDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const;
  PetscErrorCode getQuadDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const;
  PetscErrorCode getTetDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const;
  PetscErrorCode getPrismDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const;

  PetscErrorCode getEdgesDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const;
  PetscErrorCode getTrisDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const;
  PetscErrorCode getQuadDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const;
  PetscErrorCode getTetDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const;
  PetscErrorCode getPrismDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  // ** Indices **

  /// \brief get node indices
  PetscErrorCode getNodesIndices(const std::string &field_name,
    FENumeredDofEntity_multiIndex &dofs,
    VectorInt &nodes_indices,
    VectorInt &local_nodes_indices
  ) const;

  /// \brief get indices by type (generic function)
  PetscErrorCode getTypeIndices(const std::string &field_name,
    FENumeredDofEntity_multiIndex &dofs,EntityType type,int side_number,
    VectorInt &indices,
    VectorInt &local_indices
  ) const;

  /// \brief get indices by type (generic function)
  PetscErrorCode getTypeIndices(
    const std::string &field_name,FENumeredDofEntity_multiIndex &dofs,
    EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
  ) const;

  /// \brief get row node indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getRowNodesIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get col node indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getColNodesIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Edges row indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getEdgesRowIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Edges col indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getEdgesColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Tris row indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getTrisRowIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Tris col indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getTrisColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Tets row indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getTetsRowIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Tets col indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getTetsColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Quad row indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getQuadRowIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Quad col indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getQuadColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Prism row indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getPrismRowIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get Prism col indices from FENumeredDofEntity_multiIndex
  PetscErrorCode getPrismColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get NoField indices
  PetscErrorCode getNoFieldIndices(
    const std::string &field_name,FENumeredDofEntity_multiIndex &dofs,VectorInt &nodes_indices
  ) const;

  /// \brief get col NoField indices
  PetscErrorCode getNoFieldRowIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief get col NoField indices
  PetscErrorCode getNoFieldColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  // ** Data **

  /**
   * \brief Get field data on nodes
   * @param  field_name Name of field
   * @param  dofs       Dofs (element) multi index
   * @param  nodes_data Returned DOFs values
   * @param  nodes_dofs Vector of pointers to DOFs data structure
   * @param  space      Get space on nodes (Only H! is valid)
   * @param  base       Get base on nodes
   * @return            Error code
   */
  PetscErrorCode getNodesFieldData(
    const std::string &field_name,
    FEDofEntity_multiIndex &dofs,
    VectorDouble &nodes_data,
    VectorDofs &nodes_dofs,
    FieldSpace &space,
    FieldApproximationBase &base
  ) const;

  /**
   * \brief Get field data on entities
   * @param  field_name     Field name
   * @param  dofs           Dofs (element) multi index
   * @param  type           Entity type
   * @param  side_number    Side number (Local number of entity on element in canonical order)
   * @param  ent_field_data Vector of DOFs values on entities
   * @param  ent_field_dofs Vector of pointers to DOFs data structure
   * @return                Error code
   */
  PetscErrorCode getTypeFieldData(
    const std::string &field_name,
    FEDofEntity_multiIndex &dofs,
    EntityType type,
    int side_number,
    VectorDouble &ent_field_data,
    VectorDofs &ent_field_dofs
  ) const;

  PetscErrorCode getTypeFieldData(
    const std::string &field_name,
    FEDofEntity_multiIndex &dofs,
    EntityType type,
    boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
  ) const;

  PetscErrorCode getNoFieldFieldData(
    const std::string &field_name,
    FEDofEntity_multiIndex &dofs,
    VectorDouble &ent_field_data,
    VectorDofs &ent_field_dofs
  ) const;
  PetscErrorCode getNoFieldFieldData(
    DataForcesAndSurcesCore &data,const std::string &field_name
  ) const;

  /**
   * \brief Get data on nodes
   * @param  data       Data structure
   * @param  field_name Field name
   * @return            Error code
   */
  PetscErrorCode getNodesFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  PetscErrorCode getEdgesFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const;
  PetscErrorCode getTrisFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const;
  PetscErrorCode getQuadFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const;
  PetscErrorCode getTetsFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const;
  PetscErrorCode getPrismFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const;

  /// \brief Get nodes on triangles
  PetscErrorCode getFaceTriNodes(DataForcesAndSurcesCore &data) const;

  /// \brief Get field approximation space and base on entities
  PetscErrorCode getSpacesAndBaseOnEntities(DataForcesAndSurcesCore &data) const;

  // ** Data form NumeredDofEntity_multiIndex **

  /// \brief get indices of nodal indices which are declared for problem but not this particular element
  PetscErrorCode getProblemNodesIndices(const std::string &field_name,const NumeredDofEntity_multiIndex &dofs,VectorInt &nodes_indices) const;

  /// \brief get indices by type (generic function) which are declared for problem but not this particular element
  PetscErrorCode getProblemTypeIndices(
    const std::string &field_name,const NumeredDofEntity_multiIndex &dofs,
    EntityType type,int side_number,VectorInt &indices
  ) const;

  PetscErrorCode getProblemNodesRowIndices(const std::string &field_name,VectorInt &nodes_indices) const;
  PetscErrorCode getProblemTypeRowIndices(const std::string &field_name,EntityType type,int side_number,VectorInt &indices) const;
  PetscErrorCode getProblemNodesColIndices(const std::string &field_name,VectorInt &nodes_indices) const;
  PetscErrorCode getProblemTypeColIndices(const std::string &field_name,EntityType type,int side_number,VectorInt &indices) const;

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

  /// \brief It could be be removed in the future use other variant
  virtual int getRule(int order) { return 2*order; };

  virtual int getRule(
    int order_row,int order_col,int order_data
  ) {
    return getRule(order_data);
  };

  /// !\brief It will be removed in the future use other variant
  virtual PetscErrorCode setGaussPts(int order) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented");
    PetscFunctionReturn(0);
  }

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
  virtual PetscErrorCode setGaussPts(int order_row,int order_col,int order_data) {
    PetscFunctionBegin;
    ierr = setGaussPts(order_data); CHKERRQ(ierr);
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

    std::string rowFieldName;
    std::string colFieldName;
    bool doVerticesRow; ///< If false skip vertices
    bool doEdgesRow;    ///< If false skip edges
    bool doQuadsRow;
    bool doTrisRow;
    bool doTetsRow;
    bool doPrismsRow;
    bool doVerticesCol;
    bool doEdgesCol;
    bool doQuadsCol;
    bool doTrisCol;
    bool doTetsCol;
    bool doPrismsCol;
    bool sYmm;          ///< If true assume that matrix is symmetric structure

    /// set if operator is executed taking in account symmetry
    inline void setSymm() { sYmm = true; }

    /// unset if operator is executed for  non symmetric problem
    inline void unSetSymm() { sYmm = false; }

    /**
     * \brief Controls loop over entities on element
     *
     * OPRWO is used if row vector is assembled
     * OPCOL is usually ised if column vector is assembled
     * OPROWCOL is usually used for assemble matrices.
     *
     * For typical problem like Bubnov-Galrekin OPROW and OPCOL are the same. In more
     * general case for example for non-square matrices columns and rows could have
     * different numeration and/or different set of field.
     *
     */
    enum OpType {
      OPROW = 1<<0,
      OPCOL = 1<<1,
      OPROWCOL = 1<<2
    };
    char opType;

    /**
     * \brief Get operator types
     * @return Return operator type
     */
    inline int getOpType() const { return opType; }

    /**
     * \brief Set operator type
     * @param Operator type
     */
    inline void setOpType(const OpType type) { opType = type; }

    /**
     * \brief Add operator type
     */
    inline void addOpType(const OpType type) { opType |= type; }

    UserDataOperator(const std::string &field_name,const char type):
    rowFieldName(field_name),
    colFieldName(field_name),
    doVerticesRow(true),
    doEdgesRow(true),
    doQuadsRow(true),
    doTrisRow(true),
    doTetsRow(true),
    doPrismsRow(true),
    doVerticesCol(true),
    doEdgesCol(true),
    doQuadsCol(true),
    doTrisCol(true),
    doTetsCol(true),
    doPrismsCol(true),
    sYmm(true),
    opType(type),
    ptrFE(NULL) {

    };

    UserDataOperator(
      const std::string &_row_field_name,const std::string &_col_field_name,const char type
    ):
    rowFieldName(_row_field_name),
    colFieldName(_col_field_name),
    doVerticesRow(true),
    doEdgesRow(true),
    doQuadsRow(true),
    doTrisRow(true),
    doTetsRow(true),
    doPrismsRow(true),
    doVerticesCol(true),
    doEdgesCol(true),
    doQuadsCol(true),
    doTrisCol(true),
    doTetsCol(true),
    doPrismsCol(true),
    sYmm(true),
    opType(type),
    ptrFE(NULL) {}
    virtual ~UserDataOperator() {

    }

    /** \brief Return raw pointer to NumeredEntFiniteElement
     */
    inline const NumeredEntFiniteElement* getNumeredEntFiniteElementPtr() const { return ptrFE->numeredEntFiniteElementPtr; };

    /** \brief Get row indices

    Field could be or not declared for this element but is declared for problem

    \param field_name
    \param type entity type
    \param side side number, any number if type is MBVERTEX
    \return indices

    NOTE: Using those indices to assemble matrix will result in error if new non-zero values need to be created.

    */
    PetscErrorCode getPorblemRowIndices(const std::string filed_name,const EntityType type,const int side,VectorInt& indices) const;

    /** \brief Get col indices

    Field could be or not declared for this element but is declared for problem

    \param field_name
    \param type entity type
    \param side side number, any number if type is MBVERTEX
    \return indices

    NOTE: Using those indices to assemble matrix will result in error if new non-zero values need to be created.

    */
    PetscErrorCode getPorblemColIndices(const std::string filed_name,const EntityType type,const int side,VectorInt& indices) const;

    virtual PetscErrorCode setPtrFE(ForcesAndSurcesCore *ptr) {
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }

    /** \brief Return raw pointer to Finite Element Method object
     */
    inline const FEMethod* getFEMethod() { return ptrFE; }

  protected:
    ForcesAndSurcesCore *ptrFE;

  };

  boost::ptr_vector<UserDataOperator> opPtrVector;

  /** \brief Use to push back operator for row operator

   It can be used to calculate nodal forces or other quantities on the mesh.

   */
  boost::ptr_vector<UserDataOperator>& getOpPtrVector() { return opPtrVector; }

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

}

#endif //__ELEMENTSONENTITIES_HPP

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources Forces and sources
 * \ingroup mofem
 ******************************************************************************/
