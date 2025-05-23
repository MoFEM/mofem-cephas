/** \file Core.hpp
 * \brief Core interface class for user interface
 *
 * Low level data structures not used directly by user
 *
 * FIXME It is a mess with names of core cpp files need better organization
 *
 */

#ifndef __CORE_HPP__
#define __CORE_HPP__

namespace MoFEM {

/**
 * @brief Wrap MPI communicator such that is destroyed when is out of scope
 *
 */
struct WrapMPIComm {
  WrapMPIComm(MPI_Comm comm, bool petsc);
  ~WrapMPIComm();

  inline auto get_comm() { return duplicatedComm; }

private:
  MPI_Comm comm;
  MPI_Comm duplicatedComm;
  bool isPetscComm;
};

// This is to have obsolete back compatibility
struct MeshsetsManager;

template <int V> struct CoreValue {};

template <int N> struct CoreTmp : public CoreTmp<N - 1> {

  static constexpr const int value = N;
  const int getValue() const { return value; }

  boost::shared_ptr<RefEntityTmp<0>> virtual make_shared_ref_entity(
      const EntityHandle ent);

  using CoreTmp<N - 1>::CoreTmp;

  /**
   * Construct core database
   */
  CoreTmp(moab::Interface &moab,            ///< MoAB interface
          MPI_Comm comm = PETSC_COMM_WORLD, ///< MPI communicator
          const int verbose = VERBOSE       ///< Verbosity level

  );

  MoFEMErrorCode set_moab_interface(moab::Interface &new_moab, int verb);
};

template <int N> constexpr const int CoreTmp<N>::value;

/** \brief Core (interface) class
* \ingroup mofem
* \nosubgrouping

This is the implementation of abstract MoFEM::Interface class. Similarly to the
convention used in MoAB, we use small letters to name function of purely
abstract classes. This is an exception used only here. For more details about
naming functions see \ref coding_practice

This class is not used directly by the user. For internal use only.   It is
database with basic functions to access data. Abstraction of this is MoFEM
Interface structure.

Such deign to hide complexities for users and allow low development
without interfering with users modules programmer work.

\todo Implement static functions for Initialization and Finalization of MoFEM.
Those functions should keep all static variables and initialize/finalize other
libs like PETSc. Moreover initialization functions should set error handlers,
etc.

*/
template <> struct CoreTmp<0> : public Interface {

  static constexpr int value = 0;
  const int getValue() const { return value; }
  RefEntityTmp<0> getRefEntity(const EntityHandle ent) {
    return RefEntityTmp<0>(this->basicEntityDataPtr, ent);
  }

  virtual boost::shared_ptr<RefEntityTmp<0>>
  make_shared_ref_entity(const EntityHandle ent);

  /**
   * Construct core database
   */
  CoreTmp(moab::Interface &moab,            ///< MoAB interface
          MPI_Comm comm = PETSC_COMM_WORLD, ///< MPI communicator
          const int verbose = VERBOSE       ///< Verbosity level

  );

  ~CoreTmp();

  /** \name Global initialisation and finalisation  */

  /**@{*/

  /**
   * @brief Initializes the MoFEM database PETSc, MOAB and MPI.
   *
   * \note This function calls PetscInitialize, for more details see
   * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscInitialize.html>
   *
   * Example:
   * \code
   *
   * int main(int argc, char *argv[]) {
   *
   * // Initailise MoFEM and Petsc
   * MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
   *
   * try {
   *
   *   moab::Core mb_instance; // MoAB database
   *   moab::Interface &moab = mb_instance;
   *   MoFEM::Core core(moab); // MOFEM database
   *   MoFEM::CoreInterface &m_field = core;
   *
   *   CHKERR foo(); // Call function
   *
   * }
   * CATCH_ERRORS;
   *
   * return MoFEM::Core::Finalize();
   *
   * }
   *
   * \endcode
   *
   * @param argc count of number of command line arguments
   * @param args the command line arguments
   * @param file [optional] PETSc database file, also checks ~username/.petscrc
   * * and .petscrc use NULL to not check for code specific file. Use *
   * -skip_petscrc in the code specific file to skip the .petscrc files
   * @param help [optional] Help message to print, use NULL for no message
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode Initialize(int *argc, char ***args, const char file[],
                                   const char help[]);

  /**
   * @brief Checks for options to be called at the conclusion of the program.
   *
   * MPI_Finalize() is called only if the user had not called MPI_Init() before
   * calling Initialize.
   *
   * \note This function calls PetscInitialize, for more details see
   * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscFinalize.html>
   *
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode Finalize();

  /**@}*/

  /**@{*/

  /** \name Static functions */

  static void setRefEntBasicDataPtr(MoFEM::Interface &m_field,
                                    boost::shared_ptr<BasicEntityData> &ptr);

  static boost::shared_ptr<RefEntityTmp<0>>
  makeSharedRefEntity(MoFEM::Interface &m_field, const EntityHandle ent);

  /**@}&*/

  /** \name Assessing interfaces **/

  /**@{*/

  /**
   * \brief Getting interface of core database
   * @param  uuid  unique ID of interface
   * @param  iface returned pointer to interface
   * @return       error code
   */
  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  /**@}*/

  /** \name Get tag handles to data on the mesh */

  /**@{*/

  inline Tag get_th_RefParentHandle() const { return th_RefParentHandle; }
  inline Tag get_th_RefBitLevel() const { return th_RefBitLevel; }
  inline Tag get_th_RefBitEdge() const { return th_RefBitEdge; }

  /**@}*/

  /** \name Auxiliary data and functions */

  /**@{*/

  /**
   * Is used to check consistency. I n future properly this will be removed and
   * replaced by other solution. It is only for internal use.
   */
  enum SemaphoresBuildMofem {
    BUILD_FIELD = 1 << 0,
    BUILD_FE = 1 << 1,
    BUILD_ADJ = 1 << 2,
    BUILD_PROBLEM = 1 << 3,
    PARTITION_PROBLEM = 1 << 4,
    PARTITION_FE = 1 << 5,
    PARTITION_GHOST_DOFS = 1 << 6,
    PARTITION_MESH = 1 << 7
  };

  /**
   * \brief Get flags/semaphores for different stages
   */
  inline int &getBuildMoFEM() const { return *buildMoFEM; }

  /**
   * \brief add prim element
   *
   * FIXME: This is dirt solution, need to be fixed
   *
   * @param  prism prim handle
   * @param  verb  verbosity level
   * @return       error code
   */
  MoFEMErrorCode addPrismToDatabase(const EntityHandle prism,
                                    int verb = DEFAULT_VERBOSITY);

  /**@}*/

protected:
  /**
   * Construct core database
   */
  template <int V>
  CoreTmp(moab::Interface &moab, ///< MoAB interface
          MPI_Comm comm,         ///< MPI communicator
          const int verbose, CoreValue<V>);

  MoFEMErrorCode coreGenericConstructor(moab::Interface &moab, MPI_Comm comm,
                                        const int verbose);

  /** \name Tags to data on mesh and entities */

  /**@{*/

  Tag th_Part; ///< Tag for partition number
  Tag th_RefParentHandle, th_RefBitLevel, th_RefBitLevel_Mask, th_RefBitEdge,
      th_RefFEMeshset;
  Tag th_FieldId, th_FieldName, th_FieldName_DataNamePrefix, th_FieldSpace,
      th_FieldContinuity, th_FieldBase;
  Tag th_FEId, th_FEName;
  Tag th_FEIdCol, th_FEIdRow, th_FEIdData;
  Tag th_ProblemId, th_ProblemName, th_ProblemFEId;
  Tag th_ProblemNbDofsRow, th_ProblemNbDofsCol;
  Tag th_ProblemLocalNbDofRow, th_ProblemGhostNbDofRow;
  Tag th_ProblemLocalNbDofCol, th_ProblemGhostNbDofCol;
  Tag th_ProblemShift;
  Tag th_ElemType;   ///< Needed for VTK files
  Tag th_MoFEMBuild; ///< Internal use storing state, used to detect error and
                     ///< inconsistencies

  /**
   * @return pointer to BasicEntityData structure
   *
   * BasicEntityData is structure which every BasicEntity have. It keeps data
   * about tags to handles on the mesh, in particular tag to BitRefLevel and
   * tag with handle to parent.
   *
   */
  boost::shared_ptr<BasicEntityData> basicEntityDataPtr;

  /**
   * \brief Get pointer to basic entity data.
   *
   * This structure keeps data like tags handlers and other data used to
   * construct mofem entities, dofs and finite elements.
   *
   */
  boost::shared_ptr<BasicEntityData> &get_basic_entity_data_ptr() {
    return basicEntityDataPtr;
  }

  /**@}*/

  /** \name Multi-Indices accessing data on the mesh */

  /**@{*/

  RefEntity_multiIndex refinedEntities;        ///< refined entities
  RefElement_multiIndex refinedFiniteElements; ///< refined elements

  Field_multiIndex fIelds;           ///< fields
  FieldEntity_multiIndex entsFields; ///< entities on fields
  DofEntity_multiIndex dofsField;    ///< dofs on fields

  FiniteElement_multiIndex finiteElements;        ///< finite elements
  EntFiniteElement_multiIndex entsFiniteElements; ///< finite element entities

  FieldEntityEntFiniteElementAdjacencyMap_multiIndex
      entFEAdjacencies; ///< adjacencies of elements to dofs

  Problem_multiIndex pRoblems; ///< problems multi-index

  /**@}*/

  /** \name Get moab database */

  /**@{*/

  std::reference_wrapper<moab::Interface> moab; ///< moab database
  inline moab::Interface &get_moab() { return moab; }
  inline const moab::Interface &get_moab() const { return moab; }

  MoFEMErrorCode set_moab_interface(moab::Interface &new_moab,
                                    int verb = VERBOSE);

  MoFEMErrorCode setMoabInterface(moab::Interface &new_moab,
                                  int verb = VERBOSE);

  /**@}*/

  /** \name Check database consistency */

  /**@{*/

  MoFEMErrorCode
  check_number_of_ents_in_ents_field(const std::string &name) const;
  MoFEMErrorCode check_number_of_ents_in_ents_field() const;
  MoFEMErrorCode
  check_number_of_ents_in_ents_finite_element(const std::string &name) const;
  MoFEMErrorCode check_number_of_ents_in_ents_finite_element() const;

  /**@}*/

  /** \name Clear database */

  /**@{*/

  MoFEMErrorCode clear_database(int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode rebuild_database(int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /** \name Getting access to meshset manager */

  /**@{*/

  MeshsetsManager *get_meshsets_manager_ptr();
  const MeshsetsManager *get_meshsets_manager_ptr() const;
  inline MeshsetsManager &get_meshsets_manager() {
    return *get_meshsets_manager_ptr();
  }
  inline const MeshsetsManager &get_meshsets_manager() const {
    return *get_meshsets_manager_ptr();
  }

  /**@}*/

  /** \name Remove and delete entities */

  /**@{*/

  MoFEMErrorCode remove_parents_by_ents(const Range &ents,
                                        int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode remove_parents_by_bit_ref(const BitRefLevel bit,
                                           const BitRefLevel mask,
                                           int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode remove_parents_by_parents(const Range &ents,
                                           int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode remove_ents(const Range ents, int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode remove_ents_by_bit_ref(const BitRefLevel bit,
                                        const BitRefLevel mask,
                                        int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode delete_ents_by_bit_ref(const BitRefLevel bit,
                                        const BitRefLevel mask,
                                        const bool remove_parent = false,
                                        int verb = DEFAULT_VERBOSITY,
                                        MoFEMTypes mf = MF_ZERO);
  /**@}*/

  /** \name Fields */

  /**@{*/

  /**
  * \brief Add field
  * @param  name           Field name
  * @param  space          Space L2,H1,Hdiv,Hcurl
  * @param  continuity     Field continuity (you can set broken space)
  * @param  base           Approximation base AINSWORTH_LEGENDRE_BASE,
  AINSWORTH_BERNSTEIN_BEZIER_BASE ...
  * @param  nb_coefficients Number of field coefficients
  @ @param  type_dof_side_map  Map of entity type to function returning DofsSideMap
  * @param  tag_type       Tag type, MB_TAG_DENSE or MB_TAG_SPARSE (default)
  * @param  bh             Control behavior, if MF_EXCL throws error if exist
  * @param  verb           Verbosity level
  * @return                Return error code

  TODO: \todo MB_TAG_DENSE will not work properly in general case. It is need to
  separate field tags for each entity separately. That will allow for HO orders
  but homogenous approx. order on each entity. Need some discussion what is
  optimal solution. MB_TAG_SPARSE gives flexibility, but it not memory
  efficient. MB_TAG_DENSE uses memory more efficient and in principle allow for
  better efficiency if properly utilized.


  FIXME: \bug Need to resolve problem of dense tags at this stage of development
  will make only problems

  */
  virtual MoFEMErrorCode add_broken_field(
      const std::string name, const FieldSpace space,
      const FieldApproximationBase base,
      const FieldCoefficientsNumber nb_coefficients,

      std::vector<

          std::pair<EntityType,
                    std::function<MoFEMErrorCode(BaseFunction::DofsSideMap &)>

                    >>
          list_dof_side_map,

      const TagType tag_type = MB_TAG_SPARSE,
      const enum MoFEMTypes bh = MF_EXCL, int verb = DEFAULT_VERBOSITY);
  /**@{*/

  /**
  * \brief Add field
  * @param  name           Field name
  * @param  space          Space L2,H1,Hdiv,Hcurl
  * @param  continuity     Field continuity (you can set broken space)
  * @param  base           Approximation base AINSWORTH_LEGENDRE_BASE,
  AINSWORTH_BERNSTEIN_BEZIER_BASE ...
  * @param  nb_coefficients Number of field coefficients
  * @param  tag_type       Tag type, MB_TAG_DENSE or MB_TAG_SPARSE (default)
  * @param  bh             Control behavior, if MF_EXCL throws error if exist
  * @param  verb           Verbosity level
  * @return                Return error code

  TODO: \todo MB_TAG_DENSE will not work properly in general case. It is need to
  separate field tags for each entity separately. That will allow for HO orders
  but homogenous approx. order on each entity. Need some discussion what is
  optimal solution. MB_TAG_SPARSE gives flexibility, but it not memory
  efficient. MB_TAG_DENSE uses memory more efficient and in principle allow for
  better efficiency if properly utilized.


  FIXME: \bug Need to resolve problem of dense tags at this stage of development
  will make only problems

  */
  virtual MoFEMErrorCode
  add_field(const std::string name, const FieldSpace space,
            const FieldApproximationBase base,
            const FieldCoefficientsNumber nb_coefficients,
            const TagType tag_type = MB_TAG_SPARSE,
            const enum MoFEMTypes bh = MF_EXCL, int verb = DEFAULT_VERBOSITY);

  /**
   * @brief Delete field
   *
   * @param name field name
   * @param verb verbosity level
   * @return error code
   */
  MoFEMErrorCode delete_field(const std::string name,
                              int verb = DEFAULT_VERBOSITY);

  /**
   * @brief Template for add_field
   *
   * @tparam CoreN
   * @param name
   * @param space
   * @param continuity
   * @param base
   * @param nb_coefficients
   * @param tag_type
   * @param bh
   * @param verb
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode addField(const std::string &name, const FieldSpace space,
                          const FieldContinuity continuity,
                          const FieldApproximationBase base,
                          const FieldCoefficientsNumber nb_coefficients,
                          const TagType tag_type, const enum MoFEMTypes bh,
                          int verb);

  MoFEMErrorCode addEntsToFieldByDim(const Range &ents, const int dim,
                                     const std::string &name,
                                     int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode add_ents_to_field_by_dim(const Range &ents, const int dim,
                                          const std::string &name,
                                          int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode add_ents_to_field_by_type(const Range &ents,
                                           const EntityType type,
                                           const std::string &name,
                                           int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode add_ents_to_field_by_dim(const EntityHandle meshset,
                                          const int dim,
                                          const std::string &name,
                                          const bool recursive = true,
                                          int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode add_ents_to_field_by_type(const EntityHandle meshset,
                                           const EntityType type,
                                           const std::string &name,
                                           const bool recursive = true,
                                           int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode create_vertices_and_add_to_field(const std::string name,
                                                  const double coords[],
                                                  int size,
                                                  int verb = DEFAULT_VERBOSITY);

  /// \name Set approximation order

  MoFEMErrorCode setFieldOrder(const Range &ents, const BitFieldId id,
                               const ApproximationOrder order, int ver);

  MoFEMErrorCode setFieldOrderImpl(boost::shared_ptr<Field> field_ptr,
                                   const Range &ents,
                                   const ApproximationOrder order, int verb);

  MoFEMErrorCode set_field_order(const Range &ents, const BitFieldId id,
                                 const ApproximationOrder order,
                                 int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode set_field_order(const EntityHandle meshset,
                                 const EntityType type, const BitFieldId id,
                                 const ApproximationOrder order,
                                 int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode set_field_order(const Range &ents, const std::string &name,
                                 const ApproximationOrder order,
                                 int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode set_field_order(const EntityHandle meshset,
                                 const EntityType type, const std::string &name,
                                 const ApproximationOrder order,
                                 int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode set_field_order_by_entity_type_and_bit_ref(
      const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
      const BitFieldId id, const ApproximationOrder order,
      int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode set_field_order_by_entity_type_and_bit_ref(
      const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
      const std::string &name, const ApproximationOrder order,
      int verb = DEFAULT_VERBOSITY);

  /// \name Build fields

  MoFEMErrorCode
  buildFieldForNoFieldImpl(boost::shared_ptr<Field> field_ptr,
                           std::map<EntityType, int> &dof_counter, int verb);

  MoFEMErrorCode buildFieldForNoField(const BitFieldId id,
                                      std::map<EntityType, int> &dof_counter,
                                      int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode
  buildFieldForL2H1HcurlHdiv(const BitFieldId id,
                             std::map<EntityType, int> &dof_counter,
                             std::map<EntityType, int> &inactive_dof_counter,
                             int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode buildField(const boost::shared_ptr<Field> &field,
                            int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode build_fields(int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode build_field(const std::string field_name,
                             int verb = DEFAULT_VERBOSITY);

  /// \name Clear DOFs
  MoFEMErrorCode clear_inactive_dofs(int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_dofs_fields_by_bit_ref(const BitRefLevel bit,
                                              const BitRefLevel mask,
                                              int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_dofs_fields(const Range ents,
                                   int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_dofs_fields(const std::string name, const Range ents,
                                   int verb = DEFAULT_VERBOSITY);

  /// \name Clear ENTs
  MoFEMErrorCode clear_ents_fields_by_bit_ref(const BitRefLevel bit,
                                              const BitRefLevel mask,
                                              int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_ents_fields(const Range ents,
                                   int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_ents_fields(const std::string name, const Range ents,
                                   int verb = DEFAULT_VERBOSITY);

  /// \name Remove field entities

  MoFEMErrorCode
  remove_ents_from_field_by_bit_ref(const BitRefLevel bit,
                                    const BitRefLevel mask,
                                    int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode remove_ents_from_field(const std::string name,
                                        const EntityHandle meshset,
                                        const EntityType type,
                                        int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode remove_ents_from_field(const std::string name,
                                        const Range ents,
                                        int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode remove_ents_from_field(const Range ents,
                                        int verb = DEFAULT_VERBOSITY);

  /// \name Other auxiliary functions for fields

  MoFEMErrorCode list_dofs_by_field_name(const std::string &name) const;
  MoFEMErrorCode list_fields() const;
  BitFieldId get_field_id(const std::string &name) const;
  FieldBitNumber get_field_bit_number(const std::string name) const;
  std::string get_field_name(const BitFieldId id) const;
  EntityHandle get_field_meshset(const BitFieldId id) const;
  EntityHandle get_field_meshset(const std::string name) const;
  MoFEMErrorCode get_field_entities_by_dimension(const std::string name,
                                                 int dim, Range &ents) const;
  MoFEMErrorCode get_field_entities_by_type(const std::string name,
                                            EntityType type, Range &ents) const;
  MoFEMErrorCode get_field_entities_by_handle(const std::string name,
                                              Range &ents) const;
  bool check_field(const std::string &name) const;

  const Field *get_field_structure(const std::string &name,
                                   enum MoFEMTypes bh = MF_EXIST) const;

  /**@}*/

  /** \name Finite elements */

  /**@{*/

  const FiniteElement *
  get_finite_element_structure(const std::string &name,
                               enum MoFEMTypes bh = MF_EXCL) const;

  bool check_finite_element(const std::string &name) const;

  MoFEMErrorCode add_finite_element(const std::string &fe_name,
                                    enum MoFEMTypes bh = MF_EXCL,
                                    int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode
  modify_finite_element_adjacency_table(const std::string &fe_name,
                                        const EntityType type,
                                        ElementAdjacencyFunct function);
  MoFEMErrorCode
  modify_finite_element_add_field_data(const std::string &fe_name,
                                       const std::string name_filed);
  MoFEMErrorCode
  modify_finite_element_add_field_row(const std::string &fe_name,
                                      const std::string name_row);
  MoFEMErrorCode
  modify_finite_element_add_field_col(const std::string &fe_name,
                                      const std::string name_col);
  MoFEMErrorCode
  modify_finite_element_off_field_data(const std::string &fe_name,
                                       const std::string name_filed);
  MoFEMErrorCode
  modify_finite_element_off_field_row(const std::string &fe_name,
                                      const std::string name_row);
  MoFEMErrorCode
  modify_finite_element_off_field_col(const std::string &fe_name,
                                      const std::string name_col);
  MoFEMErrorCode add_ents_to_finite_element_by_type(
      const EntityHandle meshset, const EntityType type,
      const std::string name, const bool recursive = true);
  MoFEMErrorCode add_ents_to_finite_element_by_dim(const EntityHandle meshset,
                                                   const int dim,
                                                   const std::string name,
                                                   const bool recursive = true);
  MoFEMErrorCode add_ents_to_finite_element_by_type(const Range ents,
                                                    const EntityType type,
                                                    const std::string name);
  MoFEMErrorCode add_ents_to_finite_element_by_dim(const Range ents,
                                                   const int dim,
                                                   const std::string name);
  MoFEMErrorCode add_ents_to_finite_element_by_bit_ref(
      const BitRefLevel bit, const BitRefLevel mask, const std::string name,
      EntityType type, int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode
  add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,
                                        const std::string &name,
                                        const bool recursive = false);
  DEPRECATED MoFEMErrorCode add_ents_to_finite_element_EntType_by_bit_ref(
      const BitRefLevel &bit, const std::string &name, EntityType type,
      int verb = DEFAULT_VERBOSITY);
  DEPRECATED MoFEMErrorCode add_ents_to_finite_element_EntType_by_bit_ref(
      const BitRefLevel &bit, const BitRefLevel &mask, const std::string &name,
      EntityType type, int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode
  remove_ents_from_finite_element_by_bit_ref(const BitRefLevel bit,
                                             const BitRefLevel mask,
                                             int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode remove_ents_from_finite_element(const std::string name,
                                                 const EntityHandle meshset,
                                                 const EntityType type,
                                                 int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode remove_ents_from_finite_element(const std::string name,
                                                 const Range ents,
                                                 int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode remove_ents_from_finite_element(const Range ents,
                                                 int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode delete_finite_element(const std::string name,
                                       int verb = DEFAULT_VERBOSITY);

  // \name Other auxiliary functions for finite element

  /**
   * \brief Get field Id
   * @param  name field name
   * @return      field id
   */
  BitFEId getBitFEId(const std::string &fe_name) const;

  /**
   * \brief Get field name
   * @param  id field id
   * @return    field name
   */
  std::string getBitFEIdName(const BitFEId id) const;

  EntityHandle get_finite_element_meshset(const BitFEId id) const;
  EntityHandle get_finite_element_meshset(const std::string name) const;
  MoFEMErrorCode
  get_finite_element_entities_by_dimension(const std::string name, int dim,
                                           Range &ents) const;
  MoFEMErrorCode get_finite_element_entities_by_type(const std::string name,
                                                     EntityType type,
                                                     Range &ents) const;
  MoFEMErrorCode get_finite_element_entities_by_handle(const std::string name,
                                                       Range &ents) const;
  MoFEMErrorCode list_finite_elements() const;

  /**@}*/

  /** \name Problems */

  /**@{*/

  MoFEMErrorCode add_problem(const std::string &name,
                             enum MoFEMTypes bh = MF_EXCL,
                             int verb = DEFAULT_VERBOSITY);
  bool check_problem(const std::string name);
  MoFEMErrorCode delete_problem(const std::string name);
  MoFEMErrorCode
  modify_problem_add_finite_element(const std::string name_problem,
                                    const std::string &fe_name);
  MoFEMErrorCode
  modify_problem_unset_finite_element(const std::string name_problem,
                                      const std::string &fe_name);
  MoFEMErrorCode
  modify_problem_ref_level_add_bit(const std::string &name_problem,
                                   const BitRefLevel &bit);
  MoFEMErrorCode
  modify_problem_ref_level_set_bit(const std::string &name_problem,
                                   const BitRefLevel &bit);
  MoFEMErrorCode
  modify_problem_mask_ref_level_add_bit(const std::string &name_problem,
                                        const BitRefLevel &bit);
  MoFEMErrorCode
  modify_problem_mask_ref_level_set_bit(const std::string &name_problem,
                                        const BitRefLevel &bit);
  BitProblemId getBitProblemId(const std::string &name) const;
  MoFEMErrorCode list_problem() const;
  MoFEMErrorCode clear_problem(const std::string name,
                               int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_problems(int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode build_finite_elements(int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode build_finite_elements(const BitRefLevel &bit,
                                       int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode build_finite_elements(const string fe_name,
                                       const Range *const ents_ptr = nullptr,
                                       int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode buildFiniteElements(const boost::shared_ptr<FiniteElement> &fe,
                                     const Range *ents_ptr = NULL,
                                     int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_finite_elements_by_bit_ref(const BitRefLevel bit,
                                                  const BitRefLevel mask,
                                                  int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_finite_elements(const Range &ents,
                                       int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_finite_elements(const std::string &fe_name,
                                       const Range &ents,
                                       int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode
  get_problem_finite_elements_entities(const std::string name,
                                       const std::string &fe_name,
                                       const EntityHandle meshset);

  // \name Problem building (deprecated)

  DEPRECATED MoFEMErrorCode
  build_problem_on_distributed_mesh(int verb = DEFAULT_VERBOSITY);
  DEPRECATED MoFEMErrorCode build_problems(int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /** \name Adjacencies */

  /**@{*/

  MoFEMErrorCode build_adjacencies(const Range &ents,
                                   int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode build_adjacencies(const BitRefLevel &bit,
                                   int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode build_adjacencies(const BitRefLevel &bit,
                                   const BitRefLevel &mask,
                                   int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_adjacencies_entities(const BitRefLevel bit,
                                            const BitRefLevel mask,
                                            int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_adjacencies_entities(const Range ents,
                                            int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode clear_adjacencies_entities(const std::string name,
                                            const Range ents,
                                            int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode
  clear_adjacencies_finite_elements(const BitRefLevel bit,
                                    const BitRefLevel mask,
                                    int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode
  clear_adjacencies_finite_elements(const Range ents,
                                    int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode
  clear_adjacencies_finite_elements(const std::string name, const Range ents,
                                    int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /** \name Methods for preforming operations on elements */

  /**@{*/

  MoFEMErrorCode problem_basic_method_preProcess(const Problem *problem_ptr,
                                                 BasicMethod &method,
                                                 int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode
  problem_basic_method_preProcess(const std::string &problem_name,
                                  BasicMethod &method,
                                  int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode problem_basic_method_postProcess(const Problem *problem_ptr,
                                                  BasicMethod &method,
                                                  int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode
  problem_basic_method_postProcess(const std::string &problem_name,
                                   BasicMethod &method,
                                   int verb = DEFAULT_VERBOSITY);

  /**
   * @copydoc MoFEM::CoreInterface::cache_problem_entities
   */
  MoFEMErrorCode cache_problem_entities(const std::string prb_name,
                                        CacheTupleWeakPtr cache_ptr);

  MoFEMErrorCode loop_finite_elements(
      const Problem *problem_ptr, const std::string &fe_name, FEMethod &method,
      int lower_rank, int upper_rank,
      boost::shared_ptr<NumeredEntFiniteElement_multiIndex> fe_ptr = nullptr,
      MoFEMTypes bh = MF_EXIST,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr(),
      int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode loop_finite_elements(
      const std::string problem_name, const std::string &fe_name,
      FEMethod &method, int lower_rank, int upper_rank,
      boost::shared_ptr<NumeredEntFiniteElement_multiIndex> fe_ptr = nullptr,
      MoFEMTypes bh = MF_EXIST,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr(),
      int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode loop_finite_elements(
      const std::string problem_name, const std::string &fe_name,
      FEMethod &method,
      boost::shared_ptr<NumeredEntFiniteElement_multiIndex> fe_ptr = nullptr,
      MoFEMTypes bh = MF_EXIST,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr(),
      int verb = DEFAULT_VERBOSITY);

  MoFEMErrorCode loop_dofs(const Problem *problem_ptr,
                           const std::string &field_name, RowColData rc,
                           DofMethod &method, int lower_rank, int upper_rank,
                           int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode loop_dofs(const std::string &problem_name,
                           const std::string &field_name, RowColData rc,
                           DofMethod &method, int lower_rank, int upper_rank,
                           int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode loop_dofs(const std::string &problem_name,
                           const std::string &field_name, RowColData rc,
                           DofMethod &method, int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode loop_dofs(const std::string &field_name, DofMethod &method,
                           int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode loop_entities(const Problem *problem_ptr,
                               const std::string field_name, RowColData rc,
                               EntityMethod &method, int lower_rank,
                               int upper_rank, int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode loop_entities(const std::string problem_name,
                               const std::string field_name, RowColData rc,
                               EntityMethod &method, int lower_rank,
                               int upper_rank, int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode loop_entities(const std::string problem_name,
                               const std::string field_name, RowColData rc,
                               EntityMethod &method,
                               int verb = DEFAULT_VERBOSITY);
  MoFEMErrorCode loop_entities(const std::string field_name,
                               EntityMethod &method,
                               Range const *const ents = nullptr,
                               int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /** \name Accessing multi-indices */

  /**@{*/

  MoFEMErrorCode get_fields(const Field_multiIndex **fields_ptr) const;
  MoFEMErrorCode
  get_ref_ents(const RefEntity_multiIndex **refined_entities_ptr) const;
  MoFEMErrorCode get_ref_finite_elements(
      const RefElement_multiIndex **refined_finite_elements_ptr) const;
  MoFEMErrorCode
  get_finite_elements(const FiniteElement_multiIndex **fe_ptr) const;
  MoFEMErrorCode get_ents_finite_elements(
      const EntFiniteElement_multiIndex **fe_ent_ptr) const;
  MoFEMErrorCode
  get_field_ents(const FieldEntity_multiIndex **field_ents) const;
  MoFEMErrorCode get_dofs(const DofEntity_multiIndex **dofs_ptr) const;
  MoFEMErrorCode get_problem(const std::string &problem_name,
                             const Problem **problem_ptr) const;
  MoFEMErrorCode get_problems(const Problem_multiIndex **problems_ptr) const;
  MoFEMErrorCode get_ents_elements_adjacency(
      const FieldEntityEntFiniteElementAdjacencyMap_multiIndex *
          *dofs_elements_adjacency) const;

  const Field_multiIndex *get_fields() const;
  const RefEntity_multiIndex *get_ref_ents() const;
  const RefElement_multiIndex *get_ref_finite_elements() const;
  const FiniteElement_multiIndex *get_finite_elements() const;
  const EntFiniteElement_multiIndex *get_ents_finite_elements() const;
  const FieldEntity_multiIndex *get_field_ents() const;
  const DofEntity_multiIndex *get_dofs() const;
  const Problem *get_problem(const std::string problem_name) const;
  const Problem_multiIndex *get_problems() const;
  const FieldEntityEntFiniteElementAdjacencyMap_multiIndex *
  get_ents_elements_adjacency() const;

  FieldEntityByUId::iterator
  get_ent_field_by_name_begin(const std::string &field_name) const;
  FieldEntityByUId::iterator
  get_ent_field_by_name_end(const std::string &field_name) const;
  DofEntityByUId::iterator
  get_dofs_by_name_begin(const std::string &field_name) const;
  DofEntityByUId::iterator
  get_dofs_by_name_end(const std::string &field_name) const;
  DofEntityByUId::iterator
  get_dofs_by_name_and_ent_begin(const std::string &field_name,
                                 const EntityHandle ent) const;
  DofEntityByUId::iterator
  get_dofs_by_name_and_ent_end(const std::string &field_name,
                               const EntityHandle ent) const;
  DofEntityByUId::iterator
  get_dofs_by_name_and_type_begin(const std::string &field_name,
                                  const EntityType type) const;
  DofEntityByUId::iterator
  get_dofs_by_name_and_type_end(const std::string &field_name,
                                const EntityType ent) const;

  EntFiniteElement_multiIndex::index<Unique_mi_tag>::type::iterator
  get_fe_by_name_begin(const std::string &fe_name) const;
  EntFiniteElement_multiIndex::index<Unique_mi_tag>::type::iterator
  get_fe_by_name_end(const std::string &fe_name) const;

  /**@}*/

  /** \name Log events */

  /**@{*/

  // Events are are using for logging and hailed by PETSc

  PetscLogEvent MOFEM_EVENT_preProcess; ///< Event for preProcess finite element
  PetscLogEvent
      MOFEM_EVENT_operator; ///< Event for evaluating operator of finite element
  PetscLogEvent
      MOFEM_EVENT_postProcess; ///< Event for postProcess finite element
  PetscLogEvent MOFEM_EVENT_createMat;

  /**@}*/

  /** \name Communicator */

  /**@{*/

  mutable MPI_Comm mofemComm;  ///< MoFEM communicator
  mutable ParallelComm *pComm; ///< MOAB communicator structure

  int sIze; ///< MoFEM communicator size
  int rAnk; ///< MOFEM communicator rank

  /**
   * @return return communicator
   */
  inline MPI_Comm &get_comm() const { return mofemComm; }

  /**
   * @return return communicator size
   */
  inline int get_comm_size() const { return sIze; }

  /**
   * @return return communicator rank/processor
   */
  inline int get_comm_rank() const { return rAnk; }

  /**@}*/

protected:
  boost::shared_ptr<WrapMPIComm>
      wrapMPIMOABComm; ///< manage creation and destruction of MOAB communicator

  int verbose; ///< Verbosity level

  /**
   * \brief Hash map of pointers to interfaces
   */
  mutable boost::ptr_map<boost::typeindex::type_index, UnknownInterface> iFaces;

  mutable int *buildMoFEM; ///< keeps flags/semaphores for different stages

  std::string optionsPrefix; ///< Prefix for options on command line

  PetscBool initaliseAndBuildField; ///< If true build field on database
                                    ///< initialisation

  PetscBool initaliseAndBuildFiniteElements; // If true build finite elements on
                                             // database initialisation

  static bool isGloballyInitialised; ///< Core base globally initialized
  static int mpiInitialised;         ///< mpi was initialised by other agent
  static PetscBool isInitialized;    ///< petsc was initialised by other agent

  /**
   * @brief add problem
   *
   * @param id  problem id
   * @param name problem name
   * @param verb verbosity level
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode addProblem(const BitProblemId id, const std::string &name,
                            int verb = DEFAULT_VERBOSITY);

  /**
   * \brief Get tag handles
   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode getTags(int verb = DEFAULT_VERBOSITY);

  /**
   * \brief Cleaning database
   */
  MoFEMErrorCode clearMap();

  /**
   * @brief Register insterfaces
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode registerSubInterfaces();

  /**
   * \brief Return unique problem Id
   *
   * Each time this function is called, it gives new unit problem Id for bit.
   *
   */
  BitProblemId getProblemShift();

  /**
   * \brief Initialize database getting information on mesh
   */
  MoFEMErrorCode initialiseDatabaseFromMesh(int verb = DEFAULT_VERBOSITY);

  /**
   * @brief Get core options from command line
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode getOptions(int verb = DEFAULT_VERBOSITY);

  /**
   * @brief Register sub-interfaces in core interface
   *
   * @tparam IFACE
   * @return MoFEMErrorCode
   */
  template <class IFACE> MoFEMErrorCode regSubInterface();

  /**
   * @brief Register petsc events
   *
   * @tparam IFACE
   * @return MoFEMErrorCode
   */
  template <class IFACE> MoFEMErrorCode regEvents();
};

template <> struct CoreTmp<-1> : public CoreTmp<0> {

  static constexpr const int value = -1;
  const int getValue() const { return value; }

  virtual boost::shared_ptr<RefEntityTmp<0>>
  make_shared_ref_entity(const EntityHandle ent);

  /**
   * Construct core database
   */
  CoreTmp(moab::Interface &moab,            ///< MoAB interface
          MPI_Comm comm = PETSC_COMM_WORLD, ///< MPI communicator
          const int verbose = VERBOSE       ///< Verbosity level

  );

  virtual MoFEMErrorCode set_moab_interface(moab::Interface &new_moab,
                                            int verb = VERBOSE);
};

using Core = CoreTmp<0>;

} // namespace MoFEM

#include <CoreTemplates.hpp>

#endif // __CORE_HPP__
