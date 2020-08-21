/** \file Interface.hpp
 * \brief MoFEM interface
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __INTERFACE_HPP__
#define __INTERFACE_HPP__

/** \brief name space of MoFEM library functions and classes
 */
namespace MoFEM {

struct MeshsetsManager;

static const MOFEMuuid IDD_MOFEMCoreInterface =
    MOFEMuuid(BitIntefaceId(CORE_INTERFACE));
static const MOFEMuuid IDD_MOFEMDeprecatedCoreInterface =
    MOFEMuuid(BitIntefaceId(DEPRECATED_CORE_INTERFACE));

/**
 * \brief Interface
 * \nosubgrouping
 * \ingroup mofem
 *
 * This interface is used by user to: <br>
 * - create approximation fields,  <br>
 * - declare elements, <br>
 * - declare problems, <br>
 *
 *  \todo Clean interface, organize groups outsource some functions to
 * independent interface
 */
struct CoreInterface : public UnknownInterface {

  virtual ~CoreInterface() = default;

  /**
   * @brief Get the core
   *
   * @return int
   */
  virtual const int getValue() const = 0;

  /**
   * @brief Get RefEntity
   *
   * @param ent
   * @return boost::shared_ptr<RefEntityTmp<0>>
   */
  virtual boost::shared_ptr<RefEntityTmp<0>>
  make_shared_ref_entity(const EntityHandle ent) = 0;

  /** \name Interfaces */

  /**@{*/

  /**
   * get moab instance
   */
  virtual moab::Interface &get_moab() = 0;

  /**
   * get moab instance interface
   */
  virtual const moab::Interface &get_moab() const = 0;

  /**
   * @brief Set the moab interface object
   *
   * @param new_moab
   * @param verb
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode set_moab_interface(moab::Interface &new_moab,
                                            int verb = VERBOSE) = 0;

  /** \brief get MeshsetsManager pointer
   */
  virtual MeshsetsManager *get_meshsets_manager_ptr() = 0;

  /** \brief get MeshsetsManager pointer
   */
  virtual const MeshsetsManager *get_meshsets_manager_ptr() const = 0;

  /** \brief get MeshsetsManager pointer
   */
  virtual MeshsetsManager &get_meshsets_manager() = 0;

  /** \brief get MeshsetsManager pointer
   */
  virtual const MeshsetsManager &get_meshsets_manager() const = 0;

  /**@}*/

  /** \name Basic entity */

  /**@{*/

  /**
   * \brief Get pointer to basic entity data.
   *
   * This structure keeps data like tags handlers and other data used to
   * construct mofem entities, dofs and finite elements.
   *
   */
  virtual boost::shared_ptr<BasicEntityData> &get_basic_entity_data_ptr() = 0;

  /**@}*/

  /** \name Communicator */

  /**@{*/

  /**
   * get MPI communicator
   *
   */
  virtual MPI_Comm &get_comm() const = 0;

  /**
   * get communicator size
   */
  virtual int get_comm_size() const = 0;

  /**
   * get comm rank
   */
  virtual int get_comm_rank() const = 0;

  /** \name Check consistency */

  /**@{*/

  /**
   * \brief check data consistency in entitiesPtr
   *
   */
  virtual MoFEMErrorCode
  check_number_of_ents_in_ents_field(const std::string &name) const = 0;

  /**
   * \brief check data consistency in entitiesPtr
   *
   */
  virtual MoFEMErrorCode check_number_of_ents_in_ents_field() const = 0;

  /**
   * \brief check data consistency in entsFiniteElements
   *
   */
  virtual MoFEMErrorCode check_number_of_ents_in_ents_finite_element(
      const std::string &name) const = 0;

  /**
   * \brief check data consistency in entsFiniteElements
   *
   */
  virtual MoFEMErrorCode
  check_number_of_ents_in_ents_finite_element() const = 0;

  /** \name Database */

  /**@{*/

  /**
   * \brief Clear database
   * @param  verb Verbosity level
   * @return      Error code
   */
  virtual MoFEMErrorCode clear_database(int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief Clear database and initialize it once again
   * @param  verb Verbosity level
   * @return      Error code
   */
  virtual MoFEMErrorCode rebuild_database(int verb = DEFAULT_VERBOSITY) = 0;

  /**@}*/

  /** \name Delete and remove */

  /**@{*/

  /**
   * @brief Remove parents from entities
   *
   */
  virtual MoFEMErrorCode
  remove_parents_by_ents(const Range &ents, int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * @brief Remove parent from entities on bit level
   *
   * Evert entity created by refinement, split or any other change on the mesh
   * can have parent. This function remove parent from entity,
   *
   * \note Functions makeing topological changes on entities should repsect
   * parents child relation. This erase that relation. If you going to split
   * faces and create interface is recommended to call this function before
   * split opeartion.
   *
   * @param bit level
   * @param mask of bit level
   * @param verb verbosity level
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode
  remove_parents_by_bit_ref(const BitRefLevel bit, const BitRefLevel mask,
                            int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * @brief Remove paremts from entities having parents in passed range
   *
   */
  virtual MoFEMErrorCode
  remove_parents_by_parents(const Range &ents,
                            int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief delete entities form mofem and moab database
   *
   */
  virtual MoFEMErrorCode
  delete_ents_by_bit_ref(const BitRefLevel bit, const BitRefLevel mask,
                         const bool remove_parent = false,
                         int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief remove entities form mofem database
   */
  virtual MoFEMErrorCode remove_ents(const Range ents,
                                     int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief remove entities form mofem database
   */
  virtual MoFEMErrorCode
  remove_ents_by_bit_ref(const BitRefLevel bit, const BitRefLevel mask,
                         int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief delete finite element from mofem database
   */
  virtual MoFEMErrorCode
  delete_finite_element(const std::string name,
                        int verb = DEFAULT_VERBOSITY) = 0;

  /** \name Fields */

  /**@{*/

  /**
   * \brief Add field
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * @param  name              name of the filed
   * @param  space             space (L2,H1,Hdiv,Hcurl)
   * @param  base              approximation base, see FieldApproximationBase
   * @param  nb_of_coefficients number of field coefficients
   * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE
   * (DENSE is faster and uses less memory, SPARSE is more flexible if you
   * define field on subdomains)
   * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO
   * no error if field exist
   * @param  verb              verbosity level
   * @return                   error code
   */
  virtual MoFEMErrorCode
  add_field(const std::string &name, const FieldSpace space,
            const FieldApproximationBase base,
            const FieldCoefficientsNumber nb_of_coefficients,
            const TagType tag_type = MB_TAG_SPARSE,
            const enum MoFEMTypes bh = MF_EXCL,
            int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief Add entities to field meshset
   * \ingroup mofem_field
   *
   * \note not collective
   *
   * The lower dimension entities are added depending on the space type
   *
   * @param  ents range of entities
   * @param  dim  dimension of entities
   * @param  name name of field
   * @param  verb verbosity level
   * @return      error code
   */
  virtual MoFEMErrorCode
  add_ents_to_field_by_dim(const Range &ents, const int dim,
                           const std::string &name,
                           int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief Add entities to field meshset
   * \ingroup mofem_field
   *
   * \note not collective
   *
   * The lower dimension entities are added depending on the space type
   *
   * @param  ents range of entities
   * @param  type  type of entities
   * @param  name name of field
   * @param  verb verbosity level
   * @return      error code
   */
  virtual MoFEMErrorCode
  add_ents_to_field_by_type(const Range &ents, const EntityType type,
                            const std::string &name,
                            int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief Add entities to field meshset
   * \ingroup mofem_field
   *
   * \not collective
   *
   * The lower dimension entities are added depending on the space type
   *
   * @param  meshset
   * @param  dim  dimension
   * @param  name name of field
   * @param  recursive take entities recursively from embedded entities
   * @param  verb verbosity level
   * @return      error code
   */
  virtual MoFEMErrorCode
  add_ents_to_field_by_dim(const EntityHandle meshset, const int dim,
                           const std::string &name, const bool recursive = true,
                           int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief Add entities to field meshset
   * \ingroup mofem_field
   *
   * \note not collective
   *
   * The lower dimension entities are added depending on the space type
   *
   * @param  meshset
   * @param  type of entities
   * @param  name name of field
   * @param  recursive take entities recursively from embedded entities
   * @param  verb verbosity level
   * @return      error code
   */
  virtual MoFEMErrorCode
  add_ents_to_field_by_type(const EntityHandle meshset, const EntityType type,
                            const std::string &name,
                            const bool recursive = true,
                            int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * @brief Create a vertices and add to field object
   *
   * Create vertices and add them to field. Those vertices would be carring
   * DOFs of the filed.
   *
   * \note This function is typically used when NOFIELD is created, for example
   * load factor in arc-length control.
   *
   * @param name name of the field
   * @param bit bit ref level of the created vertices
   * @param coords of the vertices
   * @param size number of vertices
   * @param verb verbosity level
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode
  create_vertices_and_add_to_field(const std::string name,
                                   const double coords[], int size,
                                   int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief remove entities from field
   * \ingroup mofem_field
   *
   * \note not collective
   *
   */
  virtual MoFEMErrorCode
  remove_ents_from_field_by_bit_ref(const BitRefLevel bit,
                                    const BitRefLevel mask,
                                    int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief remove entities from field
   * \ingroup mofem_field
   *
   * \note not collective
   *
   */
  virtual MoFEMErrorCode
  remove_ents_from_field(const std::string name, const EntityHandle meshset,
                         const EntityType type,
                         int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief remove entities from field
   * \ingroup mofem_field
   *
   * \note not collective
   *
   */
  virtual MoFEMErrorCode
  remove_ents_from_field(const std::string name, const Range ents,
                         int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief remove entities from all fields
   * \ingroup mofem_field
   *
   * \note not collective
   *
   */
  virtual MoFEMErrorCode
  remove_ents_from_field(const Range ents, int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief Set order approximation of the entities in the field
   * \ingroup mofem_field
   *
   * \note not collective
   *
   * \param meshset containing set of the entities (use 0 for all the entities
   * in the meshset) \param type selected type of the entities f.e. MBTET,
   * MBTRI, MBEDGE, MBVERTEX, see moab documentation \param order
   * approximation order
   */
  virtual MoFEMErrorCode set_field_order(const EntityHandle meshset,
                                         const EntityType type,
                                         const std::string &name,
                                         const ApproximationOrder order,
                                         int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief Set order approximation of the entities in the field
   * \ingroup mofem_field
   *
   * \note not collective
   *
   * \param entities
   * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE,
   * MBVERTEX, see moab documentation \param order approximation order
   */
  virtual MoFEMErrorCode set_field_order(const Range &ents,
                                         const std::string &name,
                                         const ApproximationOrder order,
                                         int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief Set order approximation of the entities in the field
   * \ingroup mofem_field
   *
   * \note not collective
   *
   * \param bit refinement level
   * \param mask bit mask
   * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE,
   * MBVERTEX, see moab documentation \param order approximation order
   */
  virtual MoFEMErrorCode set_field_order_by_entity_type_and_bit_ref(
      const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
      const std::string &name, const ApproximationOrder order,
      int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief list entities in the field
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode list_fields() const = 0;

  /** build fields
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode build_fields(int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * @brief build field by name
   *
   * @param field_name
   * @param verb
   * m@return MoFEMErrorCode
   */
  virtual MoFEMErrorCode build_field(const std::string field_name,
                                     int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief get field bit number
  *
  * \param name of Field
  * Example:\code
  auto field_number = mField.get_field_bit_number("DISPLACEMENT");
  * \endcode
  */
  virtual FieldBitNumber get_field_bit_number(const std::string name) const = 0;

  /** \brief get field meshset
  *
  * \param name of Field
  * Example:\code
  EntityHandle disp_files_meshset = mField.get_field_meshset("DISPLACEMENT");
  * \endcode
  */
  virtual EntityHandle get_field_meshset(const std::string name) const = 0;

  /**
  * \brief get entities in the field by dimension
  * @param  name field name
  * @param  dim  dim
  * @param  ents ents
  * @return      error code

  * \ingroup mofem_field
  */
  virtual MoFEMErrorCode get_field_entities_by_dimension(const std::string name,
                                                         int dim,
                                                         Range &ents) const = 0;

  /**
  * \brief get entities in the field by type
  * @param  name field name
  * @param  type entity type
  * @param  ents ents
  * @return      error code

  * \ingroup mofem_field
  */
  virtual MoFEMErrorCode get_field_entities_by_type(const std::string name,
                                                    EntityType type,
                                                    Range &ents) const = 0;

  /**
  * \brief get entities in the field by handle
  * @param  name field name
  * @param  ents ents
  * @return      error code

  * \ingroup mofem_field
  */
  virtual MoFEMErrorCode get_field_entities_by_handle(const std::string name,
                                                      Range &ents) const = 0;

  /** \brief check if field is in database
   * \ingroup mofem_field
   *
   * \param name field name
   * \return true if field exist
   *
   */
  virtual bool check_field(const std::string &name) const = 0;

  /** \brief get field structure
   * \ingroup mofem_field
   *
   * \param name field name
   * \return const Field*
   *
   */
  virtual const Field *get_field_structure(const std::string &name) = 0;

  /**@}*/

  /** \name Finite elements */

  /**@{*/

  /**
   * \brief Check if finite element is in database
   * @param  name Name of finite element
   * @return      true if element is declared
   */
  virtual bool check_finite_element(const std::string &name) const = 0;

  /**
  * \brief add finite element
  * \ingroup mofem_fe
  * \param name finite element name
  *
  * \note add_file is a collective, should be executed on all processors.
  * Otherwise could lead to deadlock.
  *
  * Example \code
  CHKERR mField.add_finite_element("ELASTIC");
  CHKERR mField.add_finite_element("PLASTIC");
  \endcode
  */
  virtual MoFEMErrorCode add_finite_element(const std::string &fe_name,
                                            enum MoFEMTypes bh = MF_EXCL,
                                            int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief modify finite element table, only for advanced user
   * \ingroup mofem_fe
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * Using that functions means that you like to do something not usual.
   *
   */
  virtual MoFEMErrorCode
  modify_finite_element_adjacency_table(const std::string &fe_name,
                                        const EntityType type,
                                        ElementAdjacencyFunct function) = 0;

  /** \brief set finite element field data
   * \ingroup mofem_fe
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * \param name finite element name
   * \param name field name
   *
   * This function will set memory in the form of a vector
   */
  virtual MoFEMErrorCode
  modify_finite_element_add_field_data(const std::string &fe_name,
                                       const std::string &name_filed) = 0;

  /** \brief unset finite element field data
   * \ingroup mofem_fe
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * \param name finite element name
   * \param name field name
   *
   * This function will set memory in the form of a vector
   */
  virtual MoFEMErrorCode
  modify_finite_element_off_field_data(const std::string &fe_name,
                                       const std::string &name_filed) = 0;

  /** \brief set field row which finite element use
   * \ingroup mofem_fe
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * \param name finite element name
   * \param name field name
   */
  virtual MoFEMErrorCode
  modify_finite_element_add_field_row(const std::string &fe_name,
                                      const std::string &name_row) = 0;

  /** \brief unset field row which finite element use
   * \ingroup mofem_fe
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * \param name finite element name
   * \param name field name
   */
  virtual MoFEMErrorCode
  modify_finite_element_off_field_row(const std::string &fe_name,
                                      const std::string &name_row) = 0;

  /** \brief set field col which finite element use
   * \ingroup mofem_fe
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * \param name finite element name
   * \param name field name
   */
  virtual MoFEMErrorCode
  modify_finite_element_add_field_col(const std::string &fe_name,
                                      const std::string &name_row) = 0;

  /** \brief unset field col which finite element use
   * \ingroup mofem_fe
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * \param name finite element name
   * \param name field name
   */
  virtual MoFEMErrorCode
  modify_finite_element_off_field_col(const std::string &fe_name,
                                      const std::string &name_row) = 0;

  /**
   * \brief add entities to finite element
   * \ingroup mofem_fe
   *
   * \note not collective
   *
   * @param  entities   meshset or range form were entities taken
   * @param  type      type of entity
   * @param  name      name of field
   * @param  recursive take entities from meshsets in meshset
   * @return           error code
   */
  virtual MoFEMErrorCode add_ents_to_finite_element_by_type(
      const EntityHandle entities, const EntityType type,
      const std::string &name, const bool recursive = true) = 0;

  /**
   * \brief add entities to finite element
   * \ingroup mofem_fe
   *
   * \note not collective
   *
   * @param  entities  meshset or range form were entities taken
   * @param  dim       dimension
   * @param  name      name of field
   * @param  recursive take entities from meshsets in meshset
   * @return           error code
   */
  virtual MoFEMErrorCode
  add_ents_to_finite_element_by_dim(const EntityHandle entities, const int dim,
                                    const std::string &name,
                                    const bool recursive = true) = 0;

  /**
   * \brief add entities to finite elements
   * \ingroup mofem_fe
   *
   * \note not collective
   *
   * @param  ents range of entities
   * @param  type type of entity (MBVERTEX, MBEDGE, MBTRI, ...)
   * @param  name name of finite element
   * @return      error code
   */
  virtual MoFEMErrorCode
  add_ents_to_finite_element_by_type(const Range &ents, const EntityType type,
                                     const std::string &name) = 0;

  /**
   * \brief add entities to finite elements
   * \ingroup mofem_fe
   *
   * \note not collective
   *
   * @param  ents range of entities
   * @param  dim dimension of entities
   * @param  name name of finite element
   * @return      error code
   */
  virtual MoFEMErrorCode
  add_ents_to_finite_element_by_dim(const Range &ents, const int dim,
                                    const std::string &name) = 0;

  /** \brief add TET entities from given refinement level to finite element
   * database given by name
   * \ingroup mofem_fe
   *
   * \note not collective
   *
   * \param BitRefLevel bit
   * \param BitRefLevel mask
   * \param finite element name
   * \param finite element type
   * \param verbose level
   */
  virtual MoFEMErrorCode add_ents_to_finite_element_by_bit_ref(
      const BitRefLevel &bit, const BitRefLevel &mask, const std::string &name,
      EntityType type, int verb = DEFAULT_VERBOSITY) = 0;

  /** get finite element meshset
   * \ingroup mofem_fe
   *
   */
  virtual EntityHandle
  get_finite_element_meshset(const std::string &name) const = 0;

  /**
  * \brief get entities in the finite element by dimension
  * @param  name finite element name
  * @param  dim  dim
  * @param  ents ents
  * @return      error code

  * \ingroup mofem_field
  */
  virtual MoFEMErrorCode
  get_finite_element_entities_by_dimension(const std::string name, int dim,
                                           Range &ents) const = 0;

  /**
  * \brief get entities in the finite element by type
  * @param  name finite element name
  * @param  type entity type
  * @param  ents ents
  * @return      error code

  * \ingroup mofem_field
  */
  virtual MoFEMErrorCode
  get_finite_element_entities_by_type(const std::string name, EntityType type,
                                      Range &ents) const = 0;

  /**
  * \brief get entities in the finite element by handle
  * @param  name finite element name
  * @param  ents ents
  * @return      error code

  * \ingroup mofem_field
  */
  virtual MoFEMErrorCode
  get_finite_element_entities_by_handle(const std::string name,
                                        Range &ents) const = 0;

  /** \brief remove elements from given refinement level to finite element
   * database \ingroup mofem_fe
   *
   * \param BitRefLevel bit
   * \param BitRefLevel mask
   * \param verbose level
   */
  virtual MoFEMErrorCode
  remove_ents_from_finite_element_by_bit_ref(const BitRefLevel bit,
                                             const BitRefLevel mask,
                                             int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief remove entities from given refinement level to finite element
   * database
   *
   */
  virtual MoFEMErrorCode remove_ents_from_finite_element(
      const std::string name, const EntityHandle meshset, const EntityType type,
      int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief remove entities from finite element database
   * \ingroup mofem_fe
   *
   */
  virtual MoFEMErrorCode
  remove_ents_from_finite_element(const std::string name, const Range ents,
                                  int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief remove entities from finite elements in database
   * \ingroup mofem_fe
   *
   */
  virtual MoFEMErrorCode
  remove_ents_from_finite_element(const Range ents,
                                  int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief add MESHSET element to finite element database given by name
   * \ingroup mofem_fe
   *
   * \note not collective
   *
   * \param meshset contains all entities that could be used for finite
   * element \param name Finite Element name
   */
  virtual MoFEMErrorCode
  add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,
                                        const std::string &name,
                                        const bool recursive = false) = 0;

  /** \brief list finite elements in database
   * \ingroup mofem_fe
   */
  virtual MoFEMErrorCode list_finite_elements() const = 0;

  /**@}*/

  /** \name Problems */

  /**@{*/

  /** \brief Add problem
   * \ingroup mofem_problems
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   */
  virtual MoFEMErrorCode add_problem(const std::string &name,
                                     enum MoFEMTypes bh = MF_EXCL,
                                     int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief check if problem exist
   * @param  name problem name
   * @return      true if problem is in database
   */
  virtual bool check_problem(const std::string name) = 0;

  /** \brief Delete problem
   * \ingroup mofem_problems
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   */
  virtual MoFEMErrorCode delete_problem(const std::string name) = 0;

  /** \brief add finite element to problem, this add entities assigned to
   * finite element to a particular problem \ingroup mofem_problems
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * \param name Problem name
   * \param name Finite Element name
   */
  virtual MoFEMErrorCode
  modify_problem_add_finite_element(const std::string &name_problem,
                                    const std::string &fe_name) = 0;

  /** \brief unset finite element from problem, this remove entities assigned
   * to finite element to a particular problem \ingroup mofem_problems
   *
   *  Note: If problem is build, it need to be cleaned to make this effective
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * \param name Problem name
   * \param name Finite Element name
   */
  virtual MoFEMErrorCode
  modify_problem_unset_finite_element(const std::string &name_problem,
                                      const std::string &fe_name) = 0;

  /** \brief add ref level to problem
   * \ingroup mofem_problems
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * if same finite element is solved using different level of refinements,
   * than the level of refinement has to be specificied to problem in query
   *
   * \param name Problem name
   * \param BitRefLevel bitLevel
   * Example: \code
    CHKERR
  mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF1","ELASTIC");
    CHKERR
  mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF2","ELASTIC");
    CHKERR
  mField.modify_problem_ref_level_add_bit("BEAM_BENDING_ON_MESH_REF1",bit_level1);
    CHKERR
  mField.modify_problem_ref_level_add_bit("BEAM_BENDING_ON_MESH_REF2",bit_level2);
  *\endcode
  * Two Problems exist and solved independently, both are elastic, but solved
  using different mesh refinement <br>
  */
  virtual MoFEMErrorCode
  modify_problem_ref_level_add_bit(const std::string &name_problem,
                                   const BitRefLevel &bit) = 0;

  /** \brief set dof mask ref level for problem
   * \ingroup mofem_problems
   *
   */
  virtual MoFEMErrorCode
  modify_problem_mask_ref_level_add_bit(const std::string &name_problem,
                                        const BitRefLevel &bit) = 0;

  /** \brief set ref level for problem
   * \ingroup mofem_problems
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   * if same finite element is solved using different level of refinements,
   * than the level of refinement has to be specificied to problem in query
   *
   * \param name Problem name
   * \param BitRefLevel bitLevel
   * Example: \code
   CHKERR
   mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF1","ELASTIC");
   CHKERR
   mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF2","ELASTIC");

   CHKERR
   mField.modify_problem_ref_level_set_bit("BEAM_BENDING_ON_MESH_REF1",bit_level1);
   CHKERR
   mField.modify_problem_ref_level_set_bit("BEAM_BENDING_ON_MESH_REF2",bit_level2);
   * \endcode
   * Two Problems exist and solved independently, both are elastic, but solved
   * using different mesh refinement <br>
   *
   */
  virtual MoFEMErrorCode
  modify_problem_ref_level_set_bit(const std::string &name_problem,
                                   const BitRefLevel &bit) = 0;

  /** \brief set dof mask ref level for problem
   * \ingroup mofem_problems
   *
   * \note add_file is a collective, should be executed on all processors.
   * Otherwise could lead to deadlock.
   *
   */
  virtual MoFEMErrorCode
  modify_problem_mask_ref_level_set_bit(const std::string &name_problem,
                                        const BitRefLevel &bit) = 0;

  /** \brief list problems
   * \ingroup mofem_problems
   */
  virtual MoFEMErrorCode list_problem() const = 0;

  /**@}*/

  /** \name Clear dofs and entities */

  /**@{*/

  /** list dofs
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode
  list_dofs_by_field_name(const std::string &name) const = 0;

  /** Clear inactive dofs
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode clear_inactive_dofs(int verb = DEFAULT_VERBOSITY) = 0;

  /** Clear dofs by bit level
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode
  clear_dofs_fields_by_bit_ref(const BitRefLevel bit, const BitRefLevel mask,
                               int verb = DEFAULT_VERBOSITY) = 0;

  /** Clear dofs by ents
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode clear_dofs_fields(const Range ents,
                                           int verb = DEFAULT_VERBOSITY) = 0;

  /** Clear dofs by field name and ents
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode clear_dofs_fields(const std::string name,
                                           const Range ents,
                                           int verb = DEFAULT_VERBOSITY) = 0;

  /** Clear entities by field name
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode clear_ents_fields(const Range ents,
                                           int verb = DEFAULT_VERBOSITY) = 0;

  /** Clear entities by field name
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode clear_ents_fields(const std::string name,
                                           const Range ents,
                                           int verb = DEFAULT_VERBOSITY) = 0;

  /** Clear ents by bit level
   * \ingroup mofem_field
   */
  virtual MoFEMErrorCode
  clear_ents_fields_by_bit_ref(const BitRefLevel bit, const BitRefLevel mask,
                               int verb = DEFAULT_VERBOSITY) = 0;

  /**@}*/

  /** \name Build fields, finite elements and problems */

  /**@{*/

  /**
   * \brief Build finite elements
   * \ingroup mofem_fe
   *
   * Build finite element data structures. Have to be run before problem and
   * adjacencies are constructed.
   *
   * @param  verb Verbosity level
   * @return      Error code
   */
  virtual MoFEMErrorCode
  build_finite_elements(int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief Build finite elements
   * \ingroup mofem_fe
   *
   * Build finite element data structures. Have to be run before problem and
   * adjacencies are constructed.
   *
   * @param  fe_name  Name of finite element
   * @param  ents_ptr Pointer to range of finite elements
   * @param  verb     Verbosity level
   * @return      Error code
   */
  virtual MoFEMErrorCode
  build_finite_elements(const string fe_name,
                        const Range *const ents_ptr = nullptr,
                        int verb = DEFAULT_VERBOSITY) = 0;

  /**@}*/

  /** \name Clear finite elements */

  /**@{*/

  /** clear finite elements
   */
  virtual MoFEMErrorCode
  clear_finite_elements_by_bit_ref(const BitRefLevel bit,
                                   const BitRefLevel mask,
                                   int verb = DEFAULT_VERBOSITY) = 0;

  /** clear finite elements
   */
  virtual MoFEMErrorCode
  clear_finite_elements(const Range ents, int verb = DEFAULT_VERBOSITY) = 0;

  /** clear finite elements
   */
  virtual MoFEMErrorCode
  clear_finite_elements(const std::string name, const Range ents,
                        int verb = DEFAULT_VERBOSITY) = 0;

  /**@}*/

  /** \name Build adjacencies */

  /**@{*/

  /** \brief build adjacencies
   *
   * \param list of entities
   *
   * This function will get information of adjacent finite elements and fields
   * of all entities. If this is not executed, partitioning the problem is not
   * possible. Adjacency map is based on degrees of freedom adjacent to
   * elements. This linked to geometric element connectivity.
   *
   * If new degrees of freedom or new finite elements are added to the
   * database, adjacency map has to be rebuild.
   *
   */
  virtual MoFEMErrorCode build_adjacencies(const Range &ents,
                                           int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief build adjacencies
   *
   * \param bit adjacencies for refine level
   *
   * This function will get information of adjacent finite elements and fields
   * of all entities. If this is not executed, partitioning the problem is not
   * possible. Adjacency map is based on degrees of freedom adjacent to
   * elements. This linked to geometric element connectivity.
   *
   * If new degrees of freedom or new finite elements are added to the
   * database, adjacency map has to be rebuild.
   *
   */
  virtual MoFEMErrorCode build_adjacencies(const BitRefLevel &bit,
                                           int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief build adjacencies
   *
   * \param bit adjacencies for refine level
   * \param mask mask for bit level
   *
   * This function will get information of adjacent finite elements and fields
   * of all entities. If this is not executed, partitioning the problem is not
   * possible. Adjacency map is based on degrees of freedom adjacent to
   * elements. This linked to geometric element connectivity.
   *
   * If new degrees of freedom or new finite elements are added to the
   * database, adjacency map has to be rebuild.
   *
   */
  virtual MoFEMErrorCode build_adjacencies(const BitRefLevel &bit,
                                           const BitRefLevel &mask,
                                           int verb = DEFAULT_VERBOSITY) = 0;

  /**@}*/

  /** \name Clear adjacencies */

  /**@{*/

  /** \brief clear adjacency map for finite elements on given bit level
   *
   * \param bit
   * \param mask
   */
  virtual MoFEMErrorCode
  clear_adjacencies_finite_elements(const BitRefLevel bit,
                                    const BitRefLevel mask,
                                    int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief clear adjacency map for entities on given bit level
   *
   * \param bit
   * \param mask
   */
  virtual MoFEMErrorCode
  clear_adjacencies_entities(const BitRefLevel bit, const BitRefLevel mask,
                             int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief clear adjacencies for field entities by entities
   */
  virtual MoFEMErrorCode
  clear_adjacencies_entities(const Range ents,
                             int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief clear adjacencies for field entities by entities and field namd
   */
  virtual MoFEMErrorCode
  clear_adjacencies_entities(const std::string name, const Range ents,
                             int verb = DEFAULT_VERBOSITY) = 0;

  /**@}*/

  /** \name Problems */

  /**@{*/

  /** \brief clear problem
   * \ingroup mofem_problems
   */
  virtual MoFEMErrorCode clear_problem(const std::string name,
                                       int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief clear problems
   * \ingroup mofem_problems
   */
  virtual MoFEMErrorCode clear_problems(int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * \brief add finite elements to the meshset
   * \ingroup mofem_problems
   *
   * \param name is problem name
   * \param fe_name
   * \param meshset
   */
  virtual MoFEMErrorCode
  get_problem_finite_elements_entities(const std::string &name,
                                       const std::string &fe_name,
                                       const EntityHandle meshset) = 0;

  /**@}*/

  /** \name Making loops on elements and entities */

  /**@{*/

  /** \brief Set data for BasicMethod
   *
   * This function set data about problem, adjacencies and other multi-indices
   * in database. This function can be used a special case when user need to
   * do some pre- and post-processing before matrix or vector is initiated, or
   * to assemble matrix for group of FEMethods. Is used by classes
   * SnesCtx and TsCtx. Look for more details there.
   *
   * FIXME: Here we need example
   *
   * \param pointer to problem data structure
   * \param method user method derived from BasicMethod
   *
   **/
  virtual MoFEMErrorCode
  problem_basic_method_preProcess(const Problem *problem_ptr,
                                  BasicMethod &method,
                                  int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief Set data for BasicMethod
   *
   * This function set data about problem, adjacencies and other multi-indices
   * in database. This function can be used a special case when user need to
   * do some pre- and post-processing before matrix or vector is initiated, or
   * to assemble matrix for group of FEMethods. Is used by classes
   * SnesCtx and TsCtx. Look for more details there.
   *
   * FIXME: Here we need example
   *
   * \param problem_name name of the problem
   * \param method user method derived from BasicMethod
   *
   **/
  virtual MoFEMErrorCode
  problem_basic_method_preProcess(const std::string &problem_name,
                                  BasicMethod &method,
                                  int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief Set data for BasicMethod
   * \ingroup mofem_loops
   *
   * This function set data about problem, adjacencies and other multi-indices
   * in database. This function can be used a special case when user need to
   * do some pre- and post-processing before matrix or vector is initiated, or
   * to assemble matrix for group of FEMethods. Is used by classes
   * SnesCtx and TsCtx. Look for more details there.
   *
   * FIXME: Here we need example
   *
   * \param pointer to problem data structure
   * \param method user method derived from BasicMethod
   *
   **/
  virtual MoFEMErrorCode
  problem_basic_method_postProcess(const Problem *problem_ptr,
                                   BasicMethod &method,
                                   int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief Set data for BasicMethod
   * \ingroup mofem_loops
   *
   * This function set data about problem, adjacencies and other multi-indices
   * in database. This function can be used a special case when user need to
   * do some pre- and post-processing before matrix or vector is initiated, or
   * to assemble matrix for group of FEMethods. Is used by classes
   * SnesCtx and TsCtx. Look for more details there.
   *
   * FIXME: Here we need example
   *
   * \param problem_name name of the problem
   * \param method user method derived from BasicMethod
   *
   **/
  virtual MoFEMErrorCode
  problem_basic_method_postProcess(const std::string &problem_name,
                                   BasicMethod &method,
                                   int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * @brief Cache variables
   *
   * @param prb_name
   * @param cache_ptr 
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode
  cache_problem_entities(const std::string prb_name,
                         CacheTupleWeakPtr cache_ptr) = 0;

  /** \brief Make a loop over finite elements.
  *
  * This function is like swiss knife, is can be used to post-processing or
  * matrix and vectors assembly. It makes loop over given finite element for
  * given
  * problem. The particular methods executed on each element are given by
  * class derived form Interface::FEMethod. At beginning of each loop user
  * defined function (method)  preProcess() is called, for each element
  * operator() is
  * executed, at the end loop finalizes with user defined function (method)
  * postProcess().
  *
  * Methods are executed only for local elements at given processor.
  *
  * For more details pleas look to examples.
  *
  * \note If fe_ptr is given it is expected that multi-index is supbes of
  * problem multi-index. If this is not the case behavior of the code is
  * undetermined.
  *
  * @param  problem_name problem consisting set of elements
  * @param  fe_name      name of element in problem
  * @param  method       class derived form Interface::FEMethod
  * @param  fe_ptr       pointer to finite elements multi-index
  * @param  bh           if bH = MF_EXIST, throws error if fe_name does not
  exist
  * @param  cache_tuple_ptr  cache
  * @param  verb         verbosity level
  * @return              error code

  * \ingroup mofem_loops
  **/
  virtual MoFEMErrorCode loop_finite_elements(
      const std::string &problem_name, const std::string &fe_name,
      FEMethod &method,
      boost::shared_ptr<NumeredEntFiniteElement_multiIndex> fe_ptr = nullptr,
      MoFEMTypes bh = MF_EXIST,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr(),
      int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief Make a loop over finite elements on partitions from upper to
  lower rank.
  *
  * This function is like swiss knife, is can be used to post-processing or
  matrix
  * and vectors assembly. It makes loop over given finite element for given
  * problem. The particular methods executed on each element are given by
  * class derived form Interface::FEMethod. At beginning of each loop user
  defined
  * function (method)  preProcess() is called, for each element operator() is
  * executed, at the end loop finalizes with user defined function (method)
  * postProcess().
  *
  * \note If fe_ptr is given it is expected that multi-index is supbes of
  * problem multi-index. If this is not the case behavior of the code is
  * undetermined.
  *
  * For more details please look to examples.
  *
  * Interface::FEMethod

  * @param  problem_ptr pointer to problem consisting set of elements
  * @param  fe_name      name of element in problem
  * @param  method       class derived form Interface::FEMethod
  * @param  lower_rank   lower rank of process owned by finite element
  * @param  upper_rank   lower rank of process owned by finite element
  * @param  fe_ptr       pointer to finite elements multi-index
  * @param  bh           if bH = MF_EXIST, throws error if fe_name does not
  exist
  * @param  cache_data   cache data vector
  * @param  cache_row    cache row vector
  * @param  cache_col    cache row vector
  * @param  verb         verbosity level
  * @return              error code

  * \ingroup mofem_loops
  **/
  virtual MoFEMErrorCode loop_finite_elements(
      const Problem *problem_ptr, const std::string &fe_name, FEMethod &method,
      int lower_rank, int upper_rank,
      boost::shared_ptr<NumeredEntFiniteElement_multiIndex> fe_ptr = nullptr,
      MoFEMTypes bh = MF_EXIST,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr(),
      int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief Make a loop over finite elements on partitions from upper to
  lower rank.
  *
  * This function is like swiss knife, is can be used to post-processing or
  matrix
  * and vectors assembly. It makes loop over given finite element for given
  * problem. The particular methods executed on each element are given by
  * class derived form Interface::FEMethod. At beginning of each loop user
  defined
  * function (method)  preProcess() is called, for each element operator() is
  * executed, at the end loop finalizes with user defined function (method)
  * postProcess().
  *
  * \note If fe_ptr is given it is expected that multi-index is supbes of
  * problem multi-index. If this is not the case behavior of the code is
  * undetermined.
  *
  * For more details please look to examples.
  *
  * @param  problem_name pointer to problem consisting set of elements
  * @param  fe_name      name of element in problem
  * @param  method       class derived form Interface::FEMethod
  * @param  lower_rank   lower rank of process owned by finite element
  * @param  upper_rank   lower rank of process owned by finite element
  * @param  fe_ptr       pointer to finite elements multi-index
  * @param  bh           if bH = MF_EXIST, throws error if fe_name does not
  exist
  * @param  cache_data   cache data vector
  * @param  cache_row    cache row vector
  * @param  cache_col    cache row vector
  * @param  verb         verbosity level
  * @return              error code

  * \ingroup mofem_loops
  **/
  virtual MoFEMErrorCode loop_finite_elements(
      const std::string &problem_name, const std::string &fe_name,
      FEMethod &method, int lower_rank, int upper_rank,
      boost::shared_ptr<NumeredEntFiniteElement_multiIndex> fe_ptr = nullptr,
      MoFEMTypes bh = MF_EXIST,
      CacheTupleWeakPtr cache_ptr = CacheTupleSharedPtr(),
      int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief Make a loop over dofs

  * \ingroup mofem_loops
  */
  virtual MoFEMErrorCode loop_dofs(const Problem *problem_ptr,
                                   const std::string &field_name, RowColData rc,
                                   DofMethod &method, int lower_rank,
                                   int upper_rank,
                                   int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief Make a loop over dofs

  * \ingroup mofem_loops
  */
  virtual MoFEMErrorCode loop_dofs(const std::string &problem_name,
                                   const std::string &field_name, RowColData rc,
                                   DofMethod &method, int lower_rank,
                                   int upper_rank,
                                   int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief Make a loop over dofs

  * \ingroup mofem_loops
  */
  virtual MoFEMErrorCode loop_dofs(const std::string &problem_name,
                                   const std::string &field_name, RowColData rc,
                                   DofMethod &method,
                                   int verb = DEFAULT_VERBOSITY) = 0;

  /** \brief Make a loop over dofs

  * \ingroup mofem_field
  */
  virtual MoFEMErrorCode loop_dofs(const std::string &field_name,
                                   DofMethod &method,
                                   int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * @brief Loop over field entities
   * @ingroup mofem_field
   *
   * @param field_name  field entities
   * @param method user method
   * @param ents if given loop only on subset of entities in the field
   * @param verb
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode loop_entities(const std::string field_name,
                                       EntityMethod &method,
                                       Range const *const ents = nullptr,
                                       int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * @brief Loop over field entities in the problem
   * @ingroup mofem_field
   *
   * @param problem_ptr
   * @param field_name
   * @param rc
   * @param method
   * @param lower_rank
   * @param upper_rank
   * @param verb
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode loop_entities(const Problem *problem_ptr,
                                       const std::string field_name,
                                       RowColData rc, EntityMethod &method,
                                       int lower_rank, int upper_rank,
                                       int verb = DEFAULT_VERBOSITY) = 0;

  /**
   * @brief Loop over field entities in the problem
   * @ingroup mofem_field
   *
   * @param problem_name
   * @param field_name
   * @param rc
   * @param method
   * @param lower_rank
   * @param upper_rank
   * @param verb
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode loop_entities(const std::string problem_name,
                                       const std::string field_name,
                                       RowColData rc, EntityMethod &method,
                                       int lower_rank, int upper_rank,
                                       int verb = DEFAULT_VERBOSITY) = 0;
  /**
   * @brief Loop over field entities in the problem
   * @ingroup mofem_field
   *
   * @param problem_name
   * @param field_name
   * @param rc
   * @param method
   * @param verb
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode loop_entities(const std::string problem_name,
                                       const std::string field_name,
                                       RowColData rc, EntityMethod &method,
                                       int verb = DEFAULT_VERBOSITY) = 0;

  /**@}*/

  /** \name Get pointers to multi-index database */

  /**@{*/

  /**
   * @brief Get the fields object
   * @ingroup mofem_access
   *
   * @return const Field_multiIndex*
   */
  virtual const Field_multiIndex *get_fields() const = 0;

  /**
   * @brief Get the ref ents object
   * @ingroup mofem_access
   *
   * @return const RefEntity_multiIndex*
   */
  virtual const RefEntity_multiIndex *get_ref_ents() const = 0;

  /**
   * @brief Get the ref finite elements object
   * @ingroup mofem_access
   *
   * @return const RefElement_multiIndex*
   */
  virtual const RefElement_multiIndex *get_ref_finite_elements() const = 0;

  /**
   * @brief Get the finite elements object
   * @ingroup mofem_access
   *
   * @return const FiniteElement_multiIndex*
   */
  virtual const FiniteElement_multiIndex *get_finite_elements() const = 0;

  /**
   * @brief Get the ents finite elements object
   * @ingroup mofem_access
   *
   * @return const EntFiniteElement_multiIndex*
   */
  virtual const EntFiniteElement_multiIndex *
  get_ents_finite_elements() const = 0;

  /**
   * @brief Get the field ents object
   * @ingroup mofem_access
   *
   * @return const FieldEntity_multiIndex*
   */
  virtual const FieldEntity_multiIndex *get_field_ents() const = 0;

  /**
   * @brief Get the dofs object
   * @ingroup mofem_access
   *
   * @return const DofEntity_multiIndex*
   */
  virtual const DofEntity_multiIndex *get_dofs() const = 0;

  /**
   * @brief Get the problem object
   * @ingroup mofem_access
   *
   * @param problem_name
   * @return const Problem*
   */
  virtual const Problem *get_problem(const std::string &problem_name) const = 0;

  /**
   * @brief Get the problems object
   * @ingroup mofem_access
   *
   * @return const Problem_multiIndex*
   */
  virtual const Problem_multiIndex *get_problems() const = 0;

  /**
   * @brief Get the dofs elements adjacency object
   *
   * @param dofs_elements_adjacency
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode get_ents_elements_adjacency(
      const FieldEntityEntFiniteElementAdjacencyMap_multiIndex *
          *dofs_elements_adjacency) const = 0;

  /**
   * @brief Get the dofs elements adjacency object
   *
   * @return const FieldEntityEntFiniteElementAdjacencyMap_multiIndex*
   */
  virtual const FieldEntityEntFiniteElementAdjacencyMap_multiIndex *
  get_ents_elements_adjacency() const = 0;

  /** \brief Get fields multi-index from database
   * \ingroup mofem_access
   */
  virtual MoFEMErrorCode
  get_fields(const Field_multiIndex **fields_ptr) const = 0;

  /** \brief Get ref entities multi-index from database
   * \ingroup mofem_access
   */
  virtual MoFEMErrorCode
  get_ref_ents(const RefEntity_multiIndex **refined_ents_ptr) const = 0;

  /** \brief Get ref finite elements multi-index form database
   * \ingroup mofem_access
   */
  virtual MoFEMErrorCode get_ref_finite_elements(
      const RefElement_multiIndex **refined_finite_elements_ptr) const = 0;

  /** \brief Get finite elements multi-index
   * \ingroup mofem_access
   */
  virtual MoFEMErrorCode
  get_finite_elements(const FiniteElement_multiIndex **fe_ptr) const = 0;

  /** \brief Get entities finite elements multi-index
   * \ingroup mofem_access
   */
  virtual MoFEMErrorCode get_ents_finite_elements(
      const EntFiniteElement_multiIndex **fe_ent_ptr) const = 0;

  /** \brief Get problem database (data structure)
   * \ingroup mofem_access
   */
  virtual MoFEMErrorCode get_problem(const std::string &problem_name,
                                     const Problem **problem_ptr) const = 0;

  /**
   * \brief Get pointer to problems multi-index
   * \ingroup mofem_access
   */
  virtual MoFEMErrorCode
  get_problems(const Problem_multiIndex **problems_ptr) const = 0;

  /** \brief Get field multi index
   *
   * \ingroup mofem_access
   *
   */
  virtual MoFEMErrorCode
  get_field_ents(const FieldEntity_multiIndex **field_ents) const = 0;

  /** \brief Get dofs multi index
   *
   * \ingroup mofem_access
   *
   */
  virtual MoFEMErrorCode
  get_dofs(const DofEntity_multiIndex **dofs_ptr) const = 0;

  /**
   * \brief get begin iterator of filed ents of given name (instead you can
   * use _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
   *
   * \ingroup mofem_field
   *
   * for(_IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
   * 	...
   * }
   *
   * \param field_name
   */
  virtual FieldEntityByUId::iterator
  get_ent_field_by_name_begin(const std::string &field_name) const = 0;

  /**
   * \brief get begin iterator of filed dofs of given name (instead you can
   * use _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
   * \ingroup mofem_field
   *
   * for(_IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
   * 	...
   * }
   *
   * \param field_name
   */
  virtual FieldEntityByUId::iterator
  get_ent_field_by_name_end(const std::string &field_name) const = 0;

  /** \brief loop over all dofs from a moFEM field and particular field
   * \ingroup mofem_field
   */
#define _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD, NAME, IT)                  \
  auto IT = (MFIELD).get_ent_field_by_name_begin(NAME);                        \
  IT != (MFIELD).get_ent_field_by_name_end(NAME);                              \
  IT++

  /**
   * \brief get begin iterator of filed dofs of given name (instead you can
   * use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
   * \ingroup mofem_field
   *
   * for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
   * 	...
   * }
   *
   * \param field_name
   */
  virtual DofEntityByUId::iterator
  get_dofs_by_name_begin(const std::string &field_name) const = 0;

  /**
   * \brief get begin iterator of filed dofs of given name (instead you can
   * use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
   * \ingroup mofem_field
   *
   * for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
   * 	...
   * }
   *
   * \param field_name
   */
  virtual DofEntityByUId::iterator
  get_dofs_by_name_end(const std::string &field_name) const = 0;

  /** loop over all dofs from a moFEM field and particular field
   * \ingroup mofem_field
   */
#define _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD, NAME, IT)                 \
  DofEntityByUId::iterator IT = (MFIELD).get_dofs_by_name_begin(NAME);         \
  IT != (MFIELD).get_dofs_by_name_end(NAME);                                   \
  IT++

  /**
   * \brief get begin iterator of filed dofs of given name and ent(instead you
   * can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,ENT,IT)
   *
   * \ingroup mofem_field
   *
   * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(MFIELD,NAME,ENT,IT)) {
   * 	...
   * }
   *
   * \param field_name
   */
  virtual DofEntityByUId::iterator
  get_dofs_by_name_and_ent_begin(const std::string &field_name,
                                 const EntityHandle ent) const = 0;

  /**
   * \brief get begin iterator of filed dofs of given name and ent (instead
   * you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,ENT,IT)
   * \ingroup mofem_field
   *
   * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(MFIELD,NAME,ENT,IT)) {
   * 	...
   * }
   *
   * \param field_name
   */
  virtual DofEntityByUId::iterator
  get_dofs_by_name_and_ent_end(const std::string &field_name,
                               const EntityHandle ent) const = 0;

  /** \brief loop over all dofs from a moFEM field and particular field
   * \ingroup mofem_access
   */
#define _IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(MFIELD, NAME, ENT, IT)    \
  DofEntityByUId::iterator IT =                                                \
      (MFIELD).get_dofs_by_name_and_ent_begin(NAME, ENT);                      \
  IT != (MFIELD).get_dofs_by_name_and_ent_end(NAME, ENT);                      \
  IT++

  /**
   * \brief get begin iterator of filed dofs of given name and ent type
   * (instead you can use
   * _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,TYPE,IT)
   *
   * \ingroup mofem_field
   *
   * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(MFIELD,NAME,TYPE,IT)) {
   * 	...
   * }
   *
   * \param field_name
   */
  virtual DofEntityByUId::iterator
  get_dofs_by_name_and_type_begin(const std::string &field_name,
                                  const EntityType type) const = 0;

  /**
   * \brief get begin iterator of filed dofs of given name end ent
   * type(instead you can use
   * _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,TYPE,IT)
   *
   * \ingroup mofem_field
   *
   * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(MFIELD,NAME,TYPE,IT)) {
   * 	...
   * }
   *
   * \param field_name
   */
  virtual DofEntityByUId::iterator
  get_dofs_by_name_and_type_end(const std::string &field_name,
                                const EntityType type) const = 0;

  /** \brief loop over all dofs from a moFEM field and particular field
   * \ingroup mofem_field
   */
#define _IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(MFIELD, NAME, TYPE, IT)  \
  DofEntityByUId::iterator IT =                                                \
      (MFIELD).get_dofs_by_name_and_type_begin(NAME, TYPE);                    \
  IT != (MFIELD).get_dofs_by_name_and_type_end(NAME, TYPE);                    \
  IT++

  /**
   * \brief get begin iterator of finite elements of given name (instead you
   * can use _IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
   *
   * \ingroup mofem_access
   *
   * for(_IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
   * 	...
   * }
   *
   * \param fe_name
   */
  virtual EntFiniteElementByName::iterator
  get_fe_by_name_begin(const std::string &fe_name) const = 0;

  /**
   * \brief get end iterator of finite elements of given name (instead you can
   * use _IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
   *
   * \ingroup mofem_access
   *
   * for(_IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
   * 	...
   * }
   *
   * \param fe_name
   */
  virtual EntFiniteElementByName::iterator
  get_fe_by_name_end(const std::string &fe_name) const = 0;

  /** \brief loop over all finite elements from a moFEM field and FE
   * \ingroup mofem_access
   */
#define _IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD, NAME, IT)                        \
  EntFiniteElementByName::iterator IT = (MFIELD).get_fe_by_name_begin(NAME);   \
  IT != (MFIELD).get_fe_by_name_end(NAME);                                     \
  IT++
};
} // namespace MoFEM

#include <DeprecatedCoreInterface.hpp>

namespace MoFEM {

using Interface = DeprecatedCoreInterface;

} // namespace MoFEM

#endif // __INTERFACE_HPP__

/**
 * \defgroup mofem_field Fields
 * \brief Data structure for adding and managing fields
 *
 * \ingroup mofem
 ******************************************************************************/

/**
 * \defgroup mofem_ref_ents Get entities and adjacencies
 * \brief Get adjacencies/entities for given BitRefLevel (mesh refinement)
 *
 * \ingroup mofem
 ******************************************************************************/

/**
 * \defgroup mofem_fe Finite elements
 * \brief Adding and managing finite elements
 *
 * \ingroup mofem
 ******************************************************************************/

/**
 * \defgroup mofem_problems Problems
 * \brief Adding and managing problems
 *
 * \ingroup mofem
 ******************************************************************************/

/**
 * \defgroup mofem_loops Loops
 * \brief Manages complexities for integrating over finite elements and dofs.
 *
 * \ingroup mofem
 ******************************************************************************/

/**
 * \defgroup mofem_access Pointers to multi-indices
 * \brief Get direct access to multi-indices in database
 *
 * \ingroup mofem
 ******************************************************************************/
