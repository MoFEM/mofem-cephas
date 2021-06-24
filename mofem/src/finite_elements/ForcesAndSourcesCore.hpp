/** \file ForcesAndSourcesCore.hpp
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

#ifndef __FORCES_AND_SOURCES_CORE__HPP__
#define __FORCES_AND_SOURCES_CORE__HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief structure to get information form mofem into DataForcesAndSourcesCore
 * \ingroup mofem_forces_and_sources
 *
 */
struct ForcesAndSourcesCore : public FEMethod {

  Interface &mField;

  ForcesAndSourcesCore(Interface &m_field);
  typedef boost::function<int(int order_row, int order_col, int order_data)>
      RuleHookFun;

  typedef boost::function<MoFEMErrorCode(ForcesAndSourcesCore *fe_raw_ptr,
                                         int order_row, int order_col,
                                         int order_data)>
      GaussHookFun;

  /**
   * \brief Hook to get rule
   *
   * \todo check preferred format how works with gcc and clang,
   * see
   * <http://www.boost.org/doc/libs/1_64_0/doc/html/function/tutorial.html#idp247873024>
   */
  RuleHookFun getRuleHook;

  /**
   * @brief Set function to calculate integration rule
   *
   */
  GaussHookFun setRuleHook;

  /** \brief Data operator to do calculations at integration points.
    * \ingroup mofem_forces_and_sources

    Is inherited and implemented by user to do calculations. It can be used in
    many different ways but typically is used to integrate matrices (f.e.
    stiffness matrix) and the right hand vector (f.e. force vector).

    Note: It is assumed that operator is executed for symmetric problem. That
    means that is executed for not repeating entities on finite element. For
    example on triangle we have nodes, 3 edges and one face. Because of symmetry
    calculations are for: nodes-nodes, nodes-edge0, nodes-edge_1, nodes-edge_2,
    nodes-face,
    edges_1-edges_1, edges_1-edge_1, edge_1-edge_2,
    edge_1-edge_1, edge_1-edge_2,
    edge_2-edge_2,
    edge_1-face, edge_1-face, edges_3-face,
    face - face

    In case of non symmetric problem in addition calculations of lower off
    diagonal terms. F.e. edge_1-edge_0, esges_3-edge_0, edge_3-edge_1,

    In that case class variable UserDataOperator::sYmm = false;

    NoteL: By default sYmm is set for symmetric problems

  */
  struct UserDataOperator : public DataOperator {

    /**
     * \brief Controls loop over entities on element
     *
     * OPROW is used if row vector is assembled
     * OPCOL is usually used if column vector is assembled
     * OPROWCOL is usually used for assemble matrices.
     *
     * For typical problem like Bubnov-Galerkin OPROW and OPCOL are the same. In
     * more general case for example for non-square matrices columns and rows
     * could have different numeration and/or different set of fields.
     *
     */
    enum OpType {
      OPROW = 1 << 0,
      OPCOL = 1 << 1,
      OPROWCOL = 1 << 2,
      OPLAST = 1 << 3
    };

    char opType;
    std::string rowFieldName;
    std::string colFieldName;
    FieldSpace sPace;

    /**
     * This Constructor is used typically when some modification base shape
     * functions on some approx. space is applied. Operator is run for all data
     * on space.
     *
     * User has no access to field data from this operator.
     */
    UserDataOperator(const FieldSpace space, const char type = OPLAST,
                     const bool symm = true);

    UserDataOperator(const std::string &field_name, const char type,
                     const bool symm = true);

    UserDataOperator(const std::string &row_field_name,
                     const std::string &col_field_name, const char type,
                     const bool symm = true);

    /** \brief Return raw pointer to NumeredEntFiniteElement
     */
    inline boost::shared_ptr<const NumeredEntFiniteElement>
    getNumeredEntFiniteElementPtr() const;

    /**
     * \brief Return finite element entity handle
     * @return Finite element entity handle
     */
    inline EntityHandle getFEEntityHandle() const;

    /**
     * @brief Get dimension of finite element
     *
     * @return int
     */
    inline int getFEDim() const;

    /**
     * @brief Get the side number pointer
     *
     * \note For vertex is expection. Side basses in argument of function doWork
     * is zero. For other entity types side can be used as argument of this
     * function.
     *
     * @param side_number
     * @param type
     * @return boost::weak_ptr<SideNumber>
     */
    inline boost::weak_ptr<SideNumber> getSideNumberPtr(const int side_number,
                                                        const EntityType type);

    /**
     * @brief Get the side entity
     *
     * \note For vertex is expection. Side basses in argument of function
     * doWork is zero. For other entity types side can be used as argument of
     * this function.
     *
     * \code
     * MoFEMErrorCode doWork(int side, EntityType type,
     *                     DataForcesAndSourcesCore::EntData &data) {
     *  MoFEMFunctionBegin;
     *
     *  if (type == MBVERTEX) {
     *    for (int n = 0; n != number_of_nodes; ++n)
     *      EntityHandle ent = getSideEntity(n, type);
     *
     *      // Do somthing
     *
     *  } else {
     *    EntityHandle ent = getSideEntity(side, type);
     *
     *    // Do somthing
     *
     *  }
     *
     *  MoFEMFunctionReturn(0);
     * }
     * \endcode
     *
     * @param side_number
     * @param type
     */
    inline EntityHandle getSideEntity(const int side_number,
                                      const EntityType type);

    /**
     * @brief Get the number of nodes on finite element
     *
     * @return int
     */
    inline int getNumberOfNodesOnElement();

    /** \brief Get row indices

    Field could be or not declared for this element but is declared for
    problem

    \param field_name
    \param type entity type
    \param side side number, any number if type is MBVERTEX
    \return indices

    NOTE: Using those indices to assemble matrix will result in error if new
    non-zero values need to be created.

    */
    MoFEMErrorCode getProblemRowIndices(const std::string filed_name,
                                        const EntityType type, const int side,
                                        VectorInt &indices) const;

    /** \brief Get col indices

    Field could be or not declared for this element but is declared for
    problem

    \param field_name
    \param type entity type
    \param side side number, any number if type is MBVERTEX
    \return indices

    NOTE: Using those indices to assemble matrix will result in error if new
    non-zero values need to be created.

    */
    MoFEMErrorCode getProblemColIndices(const std::string filed_name,
                                        const EntityType type, const int side,
                                        VectorInt &indices) const;

    /** \brief Return raw pointer to Finite Element Method object
     */
    inline const FEMethod *getFEMethod() const;

    /**
     * \brief Get operator types
     * @return Return operator type
     */
    inline int getOpType() const;

    /**
     * \brief Set operator type
     * @param Operator type
     */
    inline void setOpType(const OpType type);

    /**
     * \brief Add operator type
     */
    inline void addOpType(const OpType type);

    /**
     * \brief get number of finite element in the loop
     * @return number of finite element
     */
    inline int getNinTheLoop() const;

    /**
     * \brief get size of elements in the loop
     * @return loop size
     */
    inline int getLoopSize() const;

    /** \brief Get name of the element
     */
    inline const std::string &getFEName() const;

    /** \name Accessing KSP */

    /**@{*/

    inline const PetscData::Switches &getDataCtx() const;

    inline const KspMethod::KSPContext getKSPCtx() const;

    inline const SnesMethod::SNESContext getSNESCtx() const;

    inline const TSMethod::TSContext getTSCtx() const;

    /**@}*/

    /**@{*/

    inline Vec getKSPf() const;

    inline Mat getKSPA() const;

    inline Mat getKSPB() const;

    /**@}*/

    /** \name Accessing SNES */

    /**@{*/

    inline Vec getSNESf() const;

    inline Vec getSNESx() const;

    inline Mat getSNESA() const;

    inline Mat getSNESB() const;

    //! \deprecated Use getSNESF intead
    DEPRECATED inline Vec getSnesF() const { return getSNESf(); }

    //! \deprecated Use getSNESX intead
    DEPRECATED inline Vec getSnesX() const { return getSNESx(); }

    //! \deprecated Use getSNESA intead
    DEPRECATED inline Mat getSnesA() const { return getSNESA(); }

    //! \deprecated Use getSNESB intead
    DEPRECATED inline Mat getSnesB() const { return getSNESB(); }

    /**@}*/

    /** \name Accessing TS */

    /**@{*/

    inline Vec getTSu() const;

    inline Vec getTSu_t() const;

    inline Vec getTSu_tt() const;

    inline Vec getTSf() const;

    inline Mat getTSA() const;

    inline Mat getTSB() const;

    inline int getTSstep() const;

    inline double getTStime() const;

    inline double getTSa() const;

    inline double getTSaa() const;

    /**@}*/

    /**@{*/

    /** \name Base funtions and integration points */

    /** \brief matrix of integration (Gauss) points for Volume Element
     *
     * For triangle: columns 0,1 are x,y coordinates respectively and column
     * 2 is a weight value for example getGaussPts()(1,13) returns y
     * coordinate of 13th Gauss point on particular volume element
     *
     * For tetrahedron: columns 0,1,2 are x,y,z coordinates respectively and
     * column 3 is a weight value for example getGaussPts()(1,13) returns y
     * coordinate of 13th Gauss point on particular volume element
     *
     */
    inline MatrixDouble &getGaussPts();

    /**
     * @brief Get integration weights
     *
     * \code
     * auto t_w = getFTensor0IntegrationWeight();
     * for(int gg = 0; gg!=getGaussPts.size2(); ++gg) {
     *  // integrate something
     *  ++t_w;
     * }
     * \endcode
     *
     * @return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
     */
    inline auto getFTensor0IntegrationWeight();

    /**@}*/

    /** \name Coordinates and access to internal data */

    /**@{*/

    /** \brief Gauss points and weight, matrix (nb. of points x 3)

      Column 0-2 integration points coordinate x, y and z, respectively. At rows
      are integration points.

    */
    inline MatrixDouble &getCoordsAtGaussPts();

    /**@}*/

    /**@{*/

    /** \name Deprecated (do not use) */

    // \deprecated Deprecated function with spelling mistake
    DEPRECATED inline MoFEMErrorCode
    getPorblemRowIndices(const std::string filed_name, const EntityType type,
                         const int side, VectorInt &indices) const;

    /**@}*/

  protected:
    ForcesAndSourcesCore *ptrFE;

    virtual MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);

    inline ForcesAndSourcesCore *getPtrFE() const;

    inline ForcesAndSourcesCore *getSidePtrFE() const;

  private:
    /**
     * @brief User call this function to loop over elements on the side of
     * face. This function calls finite element with is operator to do
     * calculations.
     *
     * @param fe_name       name of the side element
     * @param side_fe       pointer to the side element instance
     * @param dim           dimension the of side element
     * @param ent_for_side  entity handle for which adjacent volume or face will
     * be accessed
     * @return MoFEMErrorCode
     */
    MoFEMErrorCode loopSide(const string &fe_name,
                            ForcesAndSourcesCore *side_fe, const size_t dim,
                            const EntityHandle ent_for_side = 0);

    friend class ForcesAndSourcesCore;
    friend class EdgeElementForcesAndSourcesCoreBase;
    friend class FaceElementForcesAndSourcesCoreBase;
    friend class ContactPrismElementForcesAndSourcesCore;
  };

  /** \brief Use to push back operator for row operator

   It can be used to calculate nodal forces or other quantities on the mesh.

   */
  boost::ptr_vector<UserDataOperator> &getOpPtrVector() { return opPtrVector; }

  /**
   * @brief Get the Entity Polynomial Base object
   *
   * @return boost::shared_ptr<BaseFunction>&&
   */
  auto &getElementPolynomialBase() { return elementPolynomialBasePtr; }

  /**
   * @brief Get the User Polynomial Base object
   *
   * @return boost::shared_ptr<BaseFunction>&
   */
  auto &getUserPolynomialBase() { return userPolynomialBasePtr; }

  /**
   * @brief Matrix of integration points
   *
   * Columns is equal to number of integration points, numver of rows
   * depends on dimension of finite element entity, for example for
   * tetrahedron rows are x,y,z,weight. Last row is integration weight.
   *
   * FIXME: that should be moved to private class data and acessed only by
   * member function
   */
  MatrixDouble gaussPts;

  virtual MoFEMErrorCode preProcess();
  virtual MoFEMErrorCode operator()();
  virtual MoFEMErrorCode postProcess();

public:
  /** \brief Get max order of approximation for data fields

  Method  getMaxDataOrder () return maximal order on entities, for
  all data on the element. So for example if finite element is triangle, and
  triangle base function have order 4 and on edges base function have order
  2, this function return 4.

  If finite element has for example 2 or more approximated fields, for
  example Pressure (order 3) and displacement field (order 5), this function
  returns 5.

  */
  int getMaxDataOrder() const;

  /// \brief Get max order of approximation for field in rows
  int getMaxRowOrder() const;

  /// \brief Get max order of approximation for field in columns
  int getMaxColOrder() const;

  /**
   * @brief Get the entity data
   *
   * @param space
   * @param type
   * @param side
   * @return const DataForcesAndSourcesCore::EntData&
   */
  const DataForcesAndSourcesCore::EntData &getEntData(const FieldSpace space,
                                                      const EntityType type,
                                                      const int side) const {
    return dataOnElement[space]->dataOnEntities[type][side];
  }

  /**
   * @brief Get the entity data
   *
   * @param space
   * @param type
   * @param side
   * @return DataForcesAndSourcesCore::EntData&
   */
  DataForcesAndSourcesCore::EntData &
  getEntData(const FieldSpace space, const EntityType type, const int side) {
    return dataOnElement[space]->dataOnEntities[type][side];
  }

protected:
  /**
   * \brief get sense (orientation) of entity
   * @param  type type of entity
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode getEntitySense(
      const EntityType type,
      boost::ptr_vector<DataForcesAndSourcesCore::EntData> &data) const;

  /**
   * @brief Get the entity data order
   *
   * @param type
   * @param space
   * @param data
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode getEntityDataOrder(
      const EntityType type, const FieldSpace space,
      boost::ptr_vector<DataForcesAndSourcesCore::EntData> &data) const;

  /**
   * @brief Get the entity sense (orientation)
   *
   * @tparam type
   * @param data
   * @return MoFEMErrorCode
   */
  template <EntityType type>
  inline MoFEMErrorCode getEntitySense(DataForcesAndSourcesCore &data) const {
    return getEntitySense(type, data.dataOnEntities[type]);
  }

  /**
   * @brief Get the entity data order for given space
   *
   * @tparam type
   * @param data
   * @param space
   * @return MoFEMErrorCode
   */
  template <EntityType type>
  inline MoFEMErrorCode getEntityDataOrder(DataForcesAndSourcesCore &data,
                                           const FieldSpace space) const {
    return getEntityDataOrder(type, space, data.dataOnEntities[type]);
  }

  /** \name Indices */

  /**@{*/

  /// \brief get node indices
  template <typename EXTRACTOR>
  MoFEMErrorCode
  getNodesIndices(const std::string &field_name,
                  FieldEntity_vector_view &ents_field, VectorInt &nodes_indices,
                  VectorInt &local_nodes_indices, EXTRACTOR &&extractor) const;

  /// \brief get row node indices from FENumeredDofEntity_multiIndex
  MoFEMErrorCode getRowNodesIndices(DataForcesAndSourcesCore &data,
                                    const std::string &field_name) const;

  /// \brief get col node indices from FENumeredDofEntity_multiIndex
  MoFEMErrorCode getColNodesIndices(DataForcesAndSourcesCore &data,
                                    const std::string &field_name) const;

  template <typename EXTRACTOR>
  MoFEMErrorCode getEntityIndices(DataForcesAndSourcesCore &data,
                                  const std::string &field_name,
                                  FieldEntity_vector_view &ents_field,
                                  const EntityType type_lo,
                                  const EntityType type_hi,
                                  EXTRACTOR &&extractor) const;

  MoFEMErrorCode
  getEntityRowIndices(DataForcesAndSourcesCore &data,
                      const std::string &field_name,
                      const EntityType type_lo = MBVERTEX,
                      const EntityType type_hi = MBPOLYHEDRON) const;

  MoFEMErrorCode
  getEntityColIndices(DataForcesAndSourcesCore &data,
                      const std::string &field_name,
                      const EntityType type_lo = MBVERTEX,
                      const EntityType type_hi = MBPOLYHEDRON) const;

  /// \brief get NoField indices
  MoFEMErrorCode
  getNoFieldIndices(const std::string &field_name,
                    boost::shared_ptr<FENumeredDofEntity_multiIndex> dofs,
                    VectorInt &nodes_indices) const;

  /// \brief get col NoField indices
  MoFEMErrorCode getNoFieldRowIndices(DataForcesAndSourcesCore &data,
                                      const std::string &field_name) const;

  /// \brief get col NoField indices
  MoFEMErrorCode getNoFieldColIndices(DataForcesAndSourcesCore &data,
                                      const std::string &field_name) const;

  /**@}*/

  /** \name Data */

  /**@{*/

  /**
   * \brief Get field data on nodes
   */
  MoFEMErrorCode getNoFieldFieldData(const std::string field_name,
                                     VectorDouble &ent_field_data,
                                     VectorDofs &ent_field_dofs,
                                     VectorFieldEntities &ent_field) const;

  MoFEMErrorCode getNoFieldFieldData(DataForcesAndSourcesCore &data,
                                     const std::string field_name) const;

  /**
   * \brief Get data on nodes
   * @param  data       Data structure
   * @param  field_name Field name
   * @return            Error code
   */
  MoFEMErrorCode getNodesFieldData(DataForcesAndSourcesCore &data,
                                   const std::string &field_name) const;

  MoFEMErrorCode
  getEntityFieldData(DataForcesAndSourcesCore &data,
                     const std::string &field_name,
                     const EntityType type_lo = MBVERTEX,
                     const EntityType type_hi = MBPOLYHEDRON) const;

  /**@}*/

  /// \brief Get nodes on triangles
  MoFEMErrorCode getFaceTriNodes(DataForcesAndSourcesCore &data) const;

  /// \brief Get field approximation space and base on entities
  MoFEMErrorCode
  getSpacesAndBaseOnEntities(DataForcesAndSourcesCore &data) const;

  /** \name Data form NumeredDofEntity_multiIndex */

  /**@{*/

  /// \brief get indices of nodal indices which are declared for problem but
  /// not this particular element
  MoFEMErrorCode getProblemNodesIndices(const std::string &field_name,
                                        const NumeredDofEntity_multiIndex &dofs,
                                        VectorInt &nodes_indices) const;

  /// \brief get indices by type (generic function) which are declared for
  /// problem but not this particular element
  MoFEMErrorCode getProblemTypeIndices(const std::string &field_name,
                                       const NumeredDofEntity_multiIndex &dofs,
                                       EntityType type, int side_number,
                                       VectorInt &indices) const;

  MoFEMErrorCode getProblemNodesRowIndices(const std::string &field_name,
                                           VectorInt &nodes_indices) const;
  MoFEMErrorCode getProblemTypeRowIndices(const std::string &field_name,
                                          EntityType type, int side_number,
                                          VectorInt &indices) const;
  MoFEMErrorCode getProblemNodesColIndices(const std::string &field_name,
                                           VectorInt &nodes_indices) const;
  MoFEMErrorCode getProblemTypeColIndices(const std::string &field_name,
                                          EntityType type, int side_number,
                                          VectorInt &indices) const;

  /**@}*/

  /**
   * \brief another variant of getRule
   * @param  order_row  order of base function on row
   * @param  order_col  order of base function on columns
   * @param  order_data order of base function approximating data
   * @return            integration rule
   *
   * This function is overloaded by the user. The integration rule
   * is set such that specific operator implemented by the user is
   * integrated accurately. For example if user implement bilinear operator
   * \f[
   * b(u,v) =
   * \int_\mathcal{T}
   * \frac{\partial u_i}{\partial x_j}\frac{\partial v_i}{\partial x_j}
   * \textrm{d}\mathcal{T}
   * \f]
   * then if \f$u\f$ and \f$v\f$ are polynomial of given \em order, then
   * exact integral would be \code int getRule(int order) { return
   * 2*(order-1); }; \endcode
   *
   * The integration points and weights are set appropriately for given
   * entity type and integration rule from \ref quad.c
   *
   * Method \ref ForcesAndSourcesCore::getRule takes at argument takes
   * maximal polynomial order set on the element on all fields defined on
   * the element. If a user likes to have more control, another variant of
   * this function can be called which distinguishing between field orders
   * on rows, columns and data, the i.e. first argument of a bilinear form,
   * the second argument of bilinear form and field coefficients on the
   * element.
   *
   * \note If user set rule to -1 or any other negative integer, then method
   * \ref ForcesAndSourcesCore::setGaussPts is called. In that method user
   * can implement own (specific) integration method.
   *
   * \bug this function should be const
   */
  virtual int getRule(int order_row, int order_col, int order_data);

  /** \brief set user specific integration rule

    This function allows for user defined integration rule. The key is to
    called matrix gaussPts, which is used by other MoFEM procedures. Matrix
    has number of rows equal to problem dimension plus one, where last index
    is used to store weight values. %Number of columns is equal to number of
    integration points.

    \note This function is called if method \ref
    ForcesAndSourcesCore::getRule is returning integer -1 or any other
    negative integer.

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
  virtual MoFEMErrorCode setGaussPts(int order_row, int order_col,
                                     int order_data);

  /**
   * \brief Calculate base functions
   * @return Error code
   */
  MoFEMErrorCode
  calHierarchicalBaseFunctionsOnElement(const FieldApproximationBase b);

  /**
   * \brief Calculate base functions
   * @return Error code
   */
  MoFEMErrorCode calHierarchicalBaseFunctionsOnElement();

  /**
   * @brief Calculate Bernstein-Bezier base
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode calBernsteinBezierBaseFunctionsOnElement();

  /**
   * @brief Create a entity data on element object
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode createDataOnElement();

  /**
   * @brief Iterate user data operators
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode loopOverOperators();

  /**@{*/

  /** \name Deprecated (do not use) */

  /** \deprecated Use getRule(int row_order, int col_order, int data order)
   */
  virtual int getRule(int order);

  /** \deprecated setGaussPts(int row_order, int col_order, int data order);
   */
  virtual MoFEMErrorCode setGaussPts(int order);

  /**@/}*/

  /**
   * @brief Entity data on element entity rows fields
   *
   */
  const std::array<boost::shared_ptr<DataForcesAndSourcesCore>, LASTSPACE>
      dataOnElement;

  /**
   * @brief Entity data on element entity columns fields
   *
   */
  const std::array<boost::shared_ptr<DataForcesAndSourcesCore>, LASTSPACE>
      derivedDataOnElement;

  DataForcesAndSourcesCore &dataNoField;
  DataForcesAndSourcesCore &dataH1;
  DataForcesAndSourcesCore &dataHcurl;
  DataForcesAndSourcesCore &dataHdiv;
  DataForcesAndSourcesCore &dataL2;

  /**
   * @brief Vector of finite element users data operators
   *
   */
  boost::ptr_vector<UserDataOperator> opPtrVector;

  friend class UserDataOperator;

protected:
  /**
   * @brief Last evaluated type of element entity
   *
   */
  EntityType lastEvaluatedElementEntityType;

private:
  /**
   * @brief Pointer to entity polynomial base
   *
   */
  boost::shared_ptr<BaseFunction> elementPolynomialBasePtr;

  /**
   * @brief Pointer to user polynomail base
   */
  boost::shared_ptr<BaseFunction> userPolynomialBasePtr;

  /**
   * @brief Element to integrate on the sides
   *
   */
  ForcesAndSourcesCore *sidePtrFE;

  /**
   * @brief Set the pointer to face element on the side
   *
   * \note Function is is used by face element, while it iterates over
   * elements on the side
   *
   * @param side_fe_ptr
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode setSideFEPtr(const ForcesAndSourcesCore *side_fe_ptr);

  friend class VolumeElementForcesAndSourcesCoreOnSideBase;
  friend class FaceElementForcesAndSourcesCoreOnSideBase;
  friend class VolumeElementForcesAndSourcesCoreOnContactPrismSideBase;

protected:
  MatrixDouble coordsAtGaussPts; ///< coordinated at gauss points
};

/// \deprecated Used ForcesAndSourcesCore instead
DEPRECATED typedef ForcesAndSourcesCore ForcesAndSurcesCore;

boost::shared_ptr<const NumeredEntFiniteElement>
ForcesAndSourcesCore::UserDataOperator::getNumeredEntFiniteElementPtr() const {
  return ptrFE->numeredEntFiniteElementPtr;
};

EntityHandle ForcesAndSourcesCore::UserDataOperator::getFEEntityHandle() const {
  return getNumeredEntFiniteElementPtr()->getEnt();
}

int ForcesAndSourcesCore::UserDataOperator::getFEDim() const {
  return ptrFE->mField.get_moab().dimension_from_handle(getFEEntityHandle());
};

boost::weak_ptr<SideNumber>
ForcesAndSourcesCore::UserDataOperator::getSideNumberPtr(
    const int side_number, const EntityType type) {
  auto &side_table_by_side_and_type =
      ptrFE->numeredEntFiniteElementPtr->getSideNumberTable().get<1>();
  auto side_it =
      side_table_by_side_and_type.find(boost::make_tuple(type, side_number));
  if (side_it != side_table_by_side_and_type.end())
    return *side_it;
  else
    return boost::weak_ptr<SideNumber>();
}

EntityHandle
ForcesAndSourcesCore::UserDataOperator::getSideEntity(const int side_number,
                                                      const EntityType type) {
  if (auto side_ptr = getSideNumberPtr(side_number, type).lock())
    return side_ptr->ent;
  else
    return 0;
}

int ForcesAndSourcesCore::UserDataOperator::getNumberOfNodesOnElement() {
  int num_nodes;
  CHKERR ptrFE->getNumberOfNodes(num_nodes);
  return num_nodes;
}

const FEMethod *ForcesAndSourcesCore::UserDataOperator::getFEMethod() const {
  return ptrFE;
}

int ForcesAndSourcesCore::UserDataOperator::getOpType() const { return opType; }

void ForcesAndSourcesCore::UserDataOperator::setOpType(const OpType type) {
  opType = type;
}

void ForcesAndSourcesCore::UserDataOperator::addOpType(const OpType type) {
  opType |= type;
}

int ForcesAndSourcesCore::UserDataOperator::getNinTheLoop() const {
  return getFEMethod()->getNinTheLoop();
}

int ForcesAndSourcesCore::UserDataOperator::getLoopSize() const {
  return getFEMethod()->getLoopSize();
}

const std::string &ForcesAndSourcesCore::UserDataOperator::getFEName() const {
  return getFEMethod()->feName;
}

const PetscData::Switches &
ForcesAndSourcesCore::UserDataOperator::getDataCtx() const {
  return getFEMethod()->data_ctx;
}

const KspMethod::KSPContext
ForcesAndSourcesCore::UserDataOperator::getKSPCtx() const {
  return getFEMethod()->ksp_ctx;
}

const SnesMethod::SNESContext
ForcesAndSourcesCore::UserDataOperator::getSNESCtx() const {
  return getFEMethod()->snes_ctx;
}

const TSMethod::TSContext
ForcesAndSourcesCore::UserDataOperator::getTSCtx() const {
  return getFEMethod()->ts_ctx;
}

Vec ForcesAndSourcesCore::UserDataOperator::getKSPf() const {
  return getFEMethod()->ksp_f;
}

Mat ForcesAndSourcesCore::UserDataOperator::getKSPA() const {
  return getFEMethod()->ksp_A;
}

Mat ForcesAndSourcesCore::UserDataOperator::getKSPB() const {
  return getFEMethod()->ksp_B;
}

Vec ForcesAndSourcesCore::UserDataOperator::getSNESf() const {
  return getFEMethod()->snes_f;
}

Vec ForcesAndSourcesCore::UserDataOperator::getSNESx() const {
  return getFEMethod()->snes_x;
}

Mat ForcesAndSourcesCore::UserDataOperator::getSNESA() const {
  return getFEMethod()->snes_A;
}

Mat ForcesAndSourcesCore::UserDataOperator::getSNESB() const {
  return getFEMethod()->snes_B;
}

Vec ForcesAndSourcesCore::UserDataOperator::getTSu() const {
  return getFEMethod()->ts_u;
}

Vec ForcesAndSourcesCore::UserDataOperator::getTSu_t() const {
  return getFEMethod()->ts_u_t;
}

Vec ForcesAndSourcesCore::UserDataOperator::getTSu_tt() const {
  return getFEMethod()->ts_u_tt;
}

Vec ForcesAndSourcesCore::UserDataOperator::getTSf() const {
  return getFEMethod()->ts_F;
}

Mat ForcesAndSourcesCore::UserDataOperator::getTSA() const {
  return getFEMethod()->ts_A;
}

Mat ForcesAndSourcesCore::UserDataOperator::getTSB() const {
  return getFEMethod()->ts_B;
}

int ForcesAndSourcesCore::UserDataOperator::getTSstep() const {
  return getFEMethod()->ts_step;
}

double ForcesAndSourcesCore::UserDataOperator::getTStime() const {
  return getFEMethod()->ts_t;
}

double ForcesAndSourcesCore::UserDataOperator::getTSa() const {
  return getFEMethod()->ts_a;
}

double ForcesAndSourcesCore::UserDataOperator::getTSaa() const {
  return getFEMethod()->ts_aa;
}

MatrixDouble &ForcesAndSourcesCore::UserDataOperator::getGaussPts() {
  return static_cast<ForcesAndSourcesCore *>(ptrFE)->gaussPts;
}

auto ForcesAndSourcesCore::UserDataOperator::getFTensor0IntegrationWeight() {
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
      &(getGaussPts()(getGaussPts().size1() - 1, 0)));
}

MoFEMErrorCode ForcesAndSourcesCore::UserDataOperator::getPorblemRowIndices(
    const std::string filed_name, const EntityType type, const int side,
    VectorInt &indices) const {
  return getProblemRowIndices(filed_name, type, side, indices);
}

ForcesAndSourcesCore *ForcesAndSourcesCore::UserDataOperator::getPtrFE() const {
  return ptrFE;
}

ForcesAndSourcesCore *
ForcesAndSourcesCore::UserDataOperator::getSidePtrFE() const {
  return ptrFE->sidePtrFE;
}

MatrixDouble &ForcesAndSourcesCore::UserDataOperator::getCoordsAtGaussPts() {
  return static_cast<ForcesAndSourcesCore *>(ptrFE)->coordsAtGaussPts;
}

} // namespace MoFEM

#endif //__FORCES_AND_SOURCES_CORE__HPP__

/**
 * \defgroup mofem_forces_and_sources Forces and sources
 * \ingroup mofem
 *
 * \brief Manages complexities related to assembly of vector and matrices at
 * single finite element level.
 *
 **/
