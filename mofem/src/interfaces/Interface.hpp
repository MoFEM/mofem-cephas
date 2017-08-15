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

static const MOFEMuuid IDD_MOFEMInterface = MOFEMuuid( BitIntefaceId(CORE_INTERFACE) );

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
 *  \todo Clean interface, organize groups outsource some functions to independent interface
 */
struct Interface: public UnknownInterface {

  /** \name Interface */

 /**@{*/

  virtual PetscErrorCode query_interface_type(const std::type_info& type,void*& ptr) const = 0;

  template <class IFace>
  PetscErrorCode query_interface(IFace*& ptr) const {
    PetscFunctionBegin;
    void* tmp_ptr;
    ierr = query_interface_type(typeid(IFace),tmp_ptr); CHKERRQ(ierr);
    ptr = reinterpret_cast<IFace*>(tmp_ptr);
    PetscFunctionReturn(0);
  }

  template <class IFace>
  IFace* query_interface() const {
    void* tmp_ptr;
    ierr = query_interface_type(typeid(IFace),tmp_ptr); CHKERRABORT(PETSC_COMM_SELF,ierr);
    return reinterpret_cast<IFace*>(tmp_ptr);
  }

  /**
   * get moab instance
   * \copydoc MoFEM::Core::get_moab
   */
  virtual moab::Interface& get_moab() = 0;

  /**
   * get moab instance interface
   */
  virtual const moab::Interface& get_moab() const = 0;

  /** \brief get MeshsetsManager pointer
    */
  virtual MeshsetsManager* get_meshsets_manager_ptr() = 0;

  /** \brief get MeshsetsManager pointer
  */
  virtual const MeshsetsManager* get_meshsets_manager_ptr() const = 0;

  /** \brief get MeshsetsManager pointer
  */
  virtual MeshsetsManager& get_meshsets_manager() = 0;

  /** \brief get MeshsetsManager pointer
    */
  virtual const MeshsetsManager& get_meshsets_manager() const = 0;

 /**@}*/

  /** \name Basic entity */

 /**@{*/

  /**
   * \brief Get pointer to basic entity data.
   *
   * This structure keeps data like tags handlers and other data used to construct
   * mofem entities, dofs and finite elements.
   *
   */
  virtual boost::shared_ptr<BasicEntityData> get_basic_entity_data_ptr() = 0;

 /**@}*/

  /** \name Communicator */

 /**@{*/

  /**
   * get MPI communicator
   *
   */
  virtual MPI_Comm& get_comm() const = 0;

  // /**
  //  * get communicator size
  //  * \deprecated use get_comm_size
  //  */
  // DEPRECATED virtual int getCommSize() const = 0;

  // /**
  //  * get comm rank
  //  * \deprecated use get_comm_rank
  //  */
  // DEPRECATED virtual int getCommRank() const = 0;

  /**
   * get communicator size
   */
  virtual int get_comm_size() const = 0;

  /**
   * get comm rank
   */
  virtual int get_comm_rank() const = 0;

 /**@}*/

  /** \name Check consistency */

 /**@{*/

  /**
    * \brief check data consistency in entitiesPtr
    *
    */
  virtual PetscErrorCode check_number_of_ents_in_ents_field(const std::string& name) const = 0;

  /**
    * \brief check data consistency in entitiesPtr
    *
    */
  virtual PetscErrorCode check_number_of_ents_in_ents_field() const = 0;

  /**
    * \brief check data consistency in entsFiniteElements
    *
    */
  virtual PetscErrorCode check_number_of_ents_in_ents_finite_element(const std::string& name) const = 0;

  /**
    * \brief check data consistency in entsFiniteElements
    *
    */
  virtual PetscErrorCode check_number_of_ents_in_ents_finite_element() const = 0;

 /**@}*/

  /** \name Mange meshest (ALL FUNCTION DEPRECATED - DO NOT USE THIS) */

 /**@{*/

  /**
    * \brief check for CUBIT Id and CUBIT type

    \deprecated use MeshsetsManager
    \todo All cubit interface functions should be outsourced to dedicated interface

    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    */
  DEPRECATED virtual bool check_msId_meshset(const int msId,const CubitBCType cubit_bc_type) = 0;


  /**
    * \brief add cubit meshset

    \deprecated use MeshsetsManager

    * \param see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET
    * \param name of set

    */
  DEPRECATED virtual PetscErrorCode add_cubit_msId(const CubitBCType cubit_bc_tyep,const int msId,const std::string name = "") = 0;

  /**
   * \brief set attributes to cubit meshset

   \deprecated use MeshsetsManager

   * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET, NODESET, SIDESET
   * @param  ms_id         id of meshset
   * @param  attributes    attributes
   * @return               error code
   */
  DEPRECATED virtual PetscErrorCode set_cubit_msId_attribites(
    const CubitBCType cubit_bc_type,const int ms_id,const std::vector<double> &attributes,const std::string name = ""
  ) = 0;

  /**
   * \brief set (material) data structure to cubit meshset

   \deprecated use MeshsetsManager
   \todo All cubit interface functions should be outsourced to dedicated interface

   * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET, NODESET, SIDESET
   * @param  ms_id         id of meshset
   * @param  attributes    attributes
   * @return               error code
   */
  DEPRECATED virtual PetscErrorCode set_cubit_msId_attribites_data_structure(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericAttributeData &data,const std::string name = ""
  ) = 0;

  /**
   * \brief set boundary data structure to meshset

   \deprecated use MeshsetsManager
   \todo All cubit interface functions should be outsourced to dedicated interface

   * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET, NODESET, SIDESET
   * @param  ms_id         id of meshset
   * @param  data          data structure
   * @return               error code
   */
  DEPRECATED virtual PetscErrorCode set_cubit_msId_bc_data_structure(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericCubitBcData &data
  ) = 0;


  /**
    * \brief delete cubit meshset
    * \ingroup mopfem_bc

   \deprecated use MeshsetsManager
   \todo All cubit interface functions should be outsourced to dedicated interface

    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    *
    */
  DEPRECATED virtual PetscErrorCode delete_cubit_msId(const CubitBCType cubit_bc_type,const int msId) = 0;

  /**
    * \brief get cubit meshset

   \deprecated use MeshsetsManager
   \todo All cubit interface functions should be outsourced to dedicated interface

    */
  DEPRECATED virtual PetscErrorCode get_cubit_msId(const int msId,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr) = 0;

  /**
    * \brief get entities from CUBIT/meshset of a particular entity dimension

    * Nodeset can contain nodes, edges, triangles and tets. This applies to other  meshsets too.
    * The nodeset's meshset contain the nodes in the MIDDLE of the surface or volume which is done by default in Cubit,
    * Hence if all nodes on a particular nodeset are required,
    * one should get all triangles or tetrahedron for which the nodeset was create in Cubit,
    * and get all the connectivities of tris/tets.

   \deprecated use MeshsetsManager
   \todo All cubit interface functions should be outsourced to dedicated interface

   \code
   MeshsetsManager *meshset_manager_ptr;
   ierr = m_field.query_interface(meshset_manager_ptr); CHKERRQ(ierr);
   ierr = meshset_manager_ptr->getEntitiesByDimension(ms_id,cubit_bc_type,dimension,entities,true); CHKERRQ(ierr);
   \endcode

    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    * \param dimensions (0 - Nodes, 1 - Edges, 2 - Faces, 3 - Volume(tetrahedral))
    * \param Range containing the retreived entities
    * \param recursive If true, meshsets containing meshsets are queried recursively. Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  DEPRECATED virtual PetscErrorCode get_cubit_msId_entities_by_dimension(const int msId,const unsigned int cubit_bc_type, const int dimension,Range &entities,const bool recursive = false) = 0;

  /**
    * \brief get entities related to CUBIT/meshset,

    * NODESET will get Vertices only, even if the NODESET contains edges, tris and tets
    * SIDESET will get Tris, BLOCKSET will get Tets, DISPLACEMENTSET and FORCESET are stored in NODESET, PRESSURESET is stored in Sideset.

   \deprecated use MeshsetsManager
   \todo All cubit interface functions should be outsourced to dedicated interface

    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    * \param Range containing the retreived entities related to the
    * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  DEPRECATED virtual PetscErrorCode get_cubit_msId_entities_by_dimension(const int msId,const unsigned int cubit_bc_type, Range &entities,const bool recursive = false) = 0;

  /**
    * \brief get meshset from CUBIT Id and CUBIT type

   \deprecated use MeshsetsManager
   \todo All cubit interface functions should be outsourced to dedicated interface

    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    * \param meshset where to store the retrieved entities
    */
  DEPRECATED virtual PetscErrorCode get_cubit_msId_meshset(const int msId,const unsigned int cubit_bc_type,EntityHandle &meshset) = 0;

  /**
    * \brief get all CUBIT meshsets by CUBIT type

   \deprecated use MeshsetsManager
   \todo All cubit interface functions should be outsourced to dedicated interface

    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more).
    * \param meshsets is range of meshsets
    */
  DEPRECATED virtual PetscErrorCode get_cubit_meshsets(const unsigned int cubit_bc_type,Range &meshsets) = 0;

  /** \deprecated use MeshsetsManager interface instead
  */
  DEPRECATED virtual PetscErrorCode print_cubit_displacement_set() const = 0;

  /** \deprecated use MeshsetsManager interface instead
  */
  DEPRECATED virtual PetscErrorCode print_cubit_pressure_set() const = 0;

  /** \deprecated use MeshsetsManager interface instead
  */
  DEPRECATED virtual PetscErrorCode print_cubit_force_set() const = 0;

  /** \deprecated use MeshsetsManager interface instead
  */
  DEPRECATED virtual PetscErrorCode print_cubit_materials_set() const = 0;

 /**@}*/

  /** \name Database */

 /**@{*/

  /**
   * \brief Clear database
   * @param  verb Verbosity level
   * @return      Error code
   */
  virtual PetscErrorCode clear_database(int verb  = -1) = 0;

  /**
   * \brief Clear database and initialize it once again
   * @param  verb Verbosity level
   * @return      Error code
   */
  virtual PetscErrorCode rebuild_database(int verb = -1) = 0;

 /**@}*/

  /** \name Synchronize */

 /**@{*/

  /** synchronize entity range on processors (collective)

    collective - need tu be run on all processors in communicator

    */
  virtual PetscErrorCode synchronise_entities(Range &ent,int verb = -1) = 0;

  /** synchronize entity range on processors (collective)
    * \ingroup mofem_field

    collective - need tu be run on all processors in communicator

    \param id field
    \param verbose level

    */
  virtual PetscErrorCode synchronise_field_entities(const BitFieldId id,int verb = -1) = 0;

  /** synchronize entity range on processors (collective)
    * \ingroup mofem_field

    collective - need tu be run on all processors in communicator

    \param name field
    \param verbose level

    */
  virtual PetscErrorCode synchronise_field_entities(const std::string& name,int verb = -1) = 0;

 /**@}*/

  /** \name Seed entities */

 /**@{*/

  /**
  * Create finite elements based from entities in meshses. Throw error if entity is not in database
  * \todo Should be outsourced to separate interface, i.e. BitLevelManager
  *
  * \param EntityHandle meshset
  *
  */
  virtual PetscErrorCode seed_finite_elements(const EntityHandle meshset,int verb = -1) = 0;

  /**
  * Create finite elements based from entities in meshsets. Throw error if entity is not in database
  * \todo Should be outsourced to separate interface, i.e. BitLevelManager
  *
  * \param Range entities
  *
  */
  virtual PetscErrorCode seed_finite_elements(const Range &entities,int verb = -1) = 0;

  /**
  * \brief seed 2D entities (Triangles entities only) in the meshset and their adjacencies (only TRIs adjacencies) in a particular BitRefLevel
  * \todo Should be outsourced to separate interface, i.e. BitLevelManager
  *
  * \param EntityHandle MeshSet
  * \param BitRefLevel bitLevel
  *
  */
  virtual PetscErrorCode seed_ref_level_2D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1) = 0;

  /**
  * \brief seed 2D entities in the meshset and their adjacencies (only TETs adjacencies) in a particular BitRefLevel
  * \todo Should be outsourced to separate interface, i.e. BitLevelManager
  *
  * \param EntityHandle MeshSet
  * \param BitRefLevel bitLevel
  *
  * \brief Example:\code
  EntityHandle meshset1; //contains ent1,ent2,ent3
  BitRefLevel myLevel0;
  myLevel0.set(0);
  seed_ref_level_3D(meshset1,myLevel0);
  //refine meshset1 into meshset2 and get new ents which are ent4, ent5
  EntityHandle meshset2; //contains ent1,ent2,ent3,ent4,ent5
  BitRefLevel myLevel1;
  myLevel1.set(1);
  seed_ref_level_3D(meshset2,myLevel1);  \endcode

  * So entities 1,2,3 would be assigned to bit level 0 and 1 <br>
  * ent1[1,1,0,0,0,0,0], ent2[1,1,0,0,0,0,0], ent3[1,1,0,0,0,0,0], <br>
  * and entities 4 and 5 are assigned to bit level 1 only <br>
  * ent4[0,1,0,0,0,0,0], ent5[0,1,0,0,0,0,0] <br>
  *
  */
  virtual PetscErrorCode seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1) = 0;

  /**
   * \brief seed entities in the range and their adjacencies in a particular BitRefLevel
   * \todo Should be outsourced to separate interface, i.e. BitLevelManager
   */
  virtual PetscErrorCode seed_ref_level(const Range &ents,const BitRefLevel &bit,const bool only_tets = true,int verb = -1) = 0;

  /** brief seed ref level by MESHSET that contains entities other than volumes
   *
   * \param EntityHandle MeshSet
   * \param BitRefLevel bitLevel
   */
  virtual PetscErrorCode seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1) = 0;

 /**@}*/

  /** \name Getting entities by BitRefLevel */

 /**@{*/

  /**\brief add all ents from ref level given by bit to meshset
    * \todo Should be outsourced to separate interface, i.e. BitLevelManager
    * \ingroup mofem_ref_ents
    *
    * \param BitRefLevel bitLevel
    * \param BitRefLevel mask
    * \param EntityType type of entities
    * \retval EntityHandle meshset
    *
    */
  virtual PetscErrorCode get_entities_by_type_and_ref_level(
    const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb = -1
  ) = 0;

  /**\brief add all ents from ref level given by bit to meshset
   * \todo Should be outsourced to separate interface, i.e. BitLevelManager
   * \ingroup mofem_ref_ents
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityType type of entities
   * \retval ents
   *
   */
  virtual PetscErrorCode get_entities_by_type_and_ref_level(
    const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb = -1
  ) = 0;


  /**\brief add all ents from ref level given by bit to meshset
   * \todo Should be outsourced to separate interface, i.e. BitLevelManager
   * \ingroup mofem_ref_ents
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityHandle meshset
   *
   */
  virtual PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) = 0;

  /**\brief add all ents from ref level given by bit to meshset
   * \todo Should be outsourced to separate interface, i.e. BitLevelManager
   * \ingroup mofem_ref_ents
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \retval ents
   */
  virtual PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents) = 0;

 /**@}*/

  /** \name Get adjacencies */

 /**@{*/

  /** \brief Get the adjacencies associated with a entity to entities of a specified dimension.
    * \ingroup mofem_ref_ents
    * \todo Should be outsourced to separate interface, i.e. BitLevelManager
    *
    * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
    */
  virtual PetscErrorCode get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const = 0;

  /** \brief Get the adjacencies associated with a entity to entities of a specified dimension.
    * \ingroup mofem_ref_ents
    *
    * bit ref level of adjacent entities is any of bit ref level of adjacent entities
    */
  virtual PetscErrorCode get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const = 0;

  /** \brief Get the adjacencies associated with a entity to entities of a specified dimension.
    * \ingroup mofem_ref_ents
    * \todo Should be outsourced to separate interface, i.e. BitLevelManage
    *
    * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
    */
  virtual PetscErrorCode get_adjacencies(
    const Problem *problem_ptr,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type = moab::Interface::INTERSECT,
    const int verb = 0
  ) const = 0;

  /** \brief Get the adjacencies associated with a entity to entities of a specified dimension.
    * \ingroup mofem_ref_ents
    * \todo Should be outsourced to separate interface, i.e. BitLevelManage
    *
    * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
    */
  virtual PetscErrorCode get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type = moab::Interface::INTERSECT,
    const int verb = 0
  ) const = 0;

 /**@}*/

  /** \name Updating entities */

 /**@{*/

  /** \brief Get child entities form meshset containing parent entities
    * \todo Should be outsourced to separate interface, i.e. BitLevelManage
    *
    * Search for refined entities of given type whose parent are entities in the
    * parent meshset. It can be used for example to transfer information about
    * boundary conditions to refined mesh or split mesh by interface
    * elements. It is used by function refine_MESHSET, to update MESHSET finite elements.
    *
    * \param parent meshset
    * \param child_bit refinement level
    * \param type of refined entity
    * \param child_type meshset where child entities are stored (if the child meshset is set to be the parent meshset, the parent would be updated with the refined entities)
    * \param recursive if true parent meshset is searched recursively
    *
   **/
  virtual PetscErrorCode update_meshset_by_entities_children(
    const EntityHandle parent,
    const BitRefLevel &child_bit,
    const EntityHandle child,
    EntityType child_type,
    const bool recursive = false,
    int verb = -1
  ) = 0;

  /** \brief update fields meshesets by child entities
    * \ingroup mofem_field
    * \todo Should be outsourced to separate interface, i.e. BitLevelManage
    *
    */
  virtual PetscErrorCode update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb = -1) = 0;

  /** \brief update field mesheset by child entities
    * \ingroup mofem_field
    * \todo Should be outsourced to separate interface, i.e. BitLevelManage
    */
  virtual PetscErrorCode update_field_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,int verb = -1) = 0;

  /** \brief update finite element mesheset by child entities
   * \todo Should be outsourced to separate interface, i.e. BitLevelManage
   */
  virtual PetscErrorCode update_finite_element_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb = -1) = 0;

 /**@}*/

  /** \name Delete and remove */

 /**@{*/

  /** \brief delete entities form mofem and moab database
    *
    */
  virtual PetscErrorCode delete_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const bool remove_parent = false,int verb = -1) = 0;

  /** \brief remove entities form mofem database
    */
  virtual PetscErrorCode remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** \brief remove finite element from mofem database
    */
  virtual PetscErrorCode delete_finite_elements_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** \brief delete finite element from mofem database
    */
  virtual PetscErrorCode delete_finite_element(const std::string name,int verb = -1) = 0;


 /**@}*/

  /** \name Shift BitRefLevl */

 /**@{*/

  /** \brief left shift bit ref level
    * this results of deletion of entities on far left side
    * \todo Should be outsourced to separate interface, i.e. BitLevelManage
    */
  virtual PetscErrorCode shift_left_bit_ref(const int shif,int verb = -1) = 0;

  /** \brief right shift bit ref level
    * \todo Should be outsourced to separate interface, i.e. BitLevelManage
    */
  virtual PetscErrorCode shift_right_bit_ref(const int shift,int verb = -1) = 0;

 /**@}*/

  /** \name Fields */

 /**@{*/

  /**
   * \brief Add field
   * @param  name              name of the filed
   * @param  space             space (L2,H1,Hdiv,Hcurl)
   * @param  base              approximation base, see FieldApproximationBase
   * @param  nb_of_cooficients number of field coefficients
   * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE (DENSE is faster and uses less memory, SPARSE is more flexible if you define field on subdomains)
   * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO no error if field exist
   * @param  verb              verbosity leve
   * @return                   error code
   */
  virtual PetscErrorCode add_field(
    const std::string& name,
    const FieldSpace space,
    const FieldApproximationBase base,
    const FieldCoefficientsNumber nb_of_cooficients,
    const TagType tag_type = MB_TAG_SPARSE,
    const enum MoFEMTypes bh = MF_EXCL,
    int verb = -1
  ) = 0;

  // /**
  // * \deprecated
  // * \brief Add field with default AINSWORTH_LEGENDRE_BASE approximation base
  // *
  // * \note This function is deprecated, use version where base argument is spevified
  // * explicityly.
  // *
  // * @param  name              name of the filed
  // * @param  space             space (L2,H1,Hdiv,Hcurl)
  // * @param  nb_of_cooficients number of field coefficients
  // * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO no error if field exist
  // * @param  verb              verbosity level
  // * @return                   error code
  // */
  // DEPRECATED inline PetscErrorCode add_field(
  //   const std::string& name,
  //   const FieldSpace space,
  //   const FieldCoefficientsNumber nb_of_cooficients,
  //   const enum MoFEMTypes bh = MF_EXCL,
  //   int verb = -1
  // ) {
  //   return add_field(
  //     name,
  //     space,
  //     AINSWORTH_LEGENDRE_BASE,
  //     nb_of_cooficients,
  //     MB_TAG_SPARSE,bh,verb
  //   );
  // }

  /**
   * \brief Add entities to field meshset
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   *
   * @param  ents rnage of entities
   * @param  dim  dimension of entities
   * @param  name name of field
   * @param  verb verbosity level
   * @return      error code
   */
  virtual PetscErrorCode add_ents_to_field_by_dim(
    const Range &ents,const int dim,const std::string& name,int verb = -1
  ) = 0;

  /**
   * \brief Add entities to field meshset
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   *
   * @param  ents rnage of entities
   * @param  type  type of entities
   * @param  name name of field
   * @param  verb verbosity level
   * @return      error code
   */
  virtual PetscErrorCode add_ents_to_field_by_type(
    const Range &ents,const EntityType type,const std::string& name,int verb = -1
  ) = 0;

  /**
   * \brief Add entities to field meshset
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   *
   * @param  meshset
   * @param  dim  diemsipm
   * @param  name name of field
   * @param  recursive take entities recursively from embaded entities
   * @param  verb verbosity level
   * @return      error code
   */
  virtual PetscErrorCode add_ents_to_field_by_dim(
    const EntityHandle meshset,const int dim,const std::string& name,const bool recursive = true,int verb = -1
  ) = 0;

  /**
   * \brief Add entities to field meshset
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   *
   * @param  meshset
   * @param  type of entoties
   * @param  name name of field
   * @param  recursive take entities recursively from embaded entities
   * @param  verb verbosity level
   * @return      error code
   */
  virtual PetscErrorCode add_ents_to_field_by_type(
    const EntityHandle meshset,const EntityType type,const std::string& name,const bool recursive = true,int verb = -1
  ) = 0;

  /**
    * \brief set field entities on vertices
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param nodes contains set vertices
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_VERTICEs(const Range &nodes,const std::string& name,int verb = -1) = 0;

  /**
    * \brief set field entities on vertices
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set vertices
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const std::string& name,int verb = -1) = 0;

  /**
    * \brief set field entities form adjacencies of edges
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param tange contains set edges
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_EDGEs(const Range &edges,const std::string& name,int verb = -1) = 0;

  /**
    * \brief set field entities form adjacencies of edges
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set edges
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_EDGEs(const EntityHandle meshset,const std::string& name,int verb = -1) = 0;

  /**
    * \brief set field entities form adjacencies of triangles
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set triangles
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TRIs(const EntityHandle meshset,const std::string& name,int verb = -1) = 0;

  /**
    * \brief set field entities form adjacencies of triangles
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param range triangles
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TRIs(const Range &tris,const std::string& name,int verb = -1) = 0;

  /**
    * \brief set field entities from adjacencies of tetrahedron
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set tetrahedron
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const std::string& name,int verb = -1) = 0;

  /**
    * \brief set field entities from adjacencies of tetrahedron
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param range contains set tetrahedron
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TETs(const Range &tets,const std::string& name,int verb = -1) = 0;

  // /**
  //   * \brief set field entities from adjacencies of quads
  //   * \ingroup mofem_field
  //   *
  //   * The lower dimension entities are added depending on the space type
  //   * \param quads range of quads
  //   * \param id field id
  //   */
  // virtual PetscErrorCode add_ents_to_field_by_QUADs(const Range &quads,const BitFieldId id,int verb = -1) = 0;

  /**
    * \brief set field entities from adjacencies of quads
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param quads range contains set quads
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_QUADs(const Range &quads,const std::string& name,int verb = -1) = 0;

  /**
    * \brief set field entities from adjacencies of quads
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set quads
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_QUADs(EntityHandle meshset,const std::string& name,int verb = -1) = 0;

  // /**
  //   * \brief set field entities from adjacencies of prisms
  //   * \ingroup mofem_field
  //   *
  //   * The lower dimension entities are added depending on the space type
  //   * \param prisms range of prisms
  //   * \param id field id
  //   */
  // virtual PetscErrorCode add_ents_to_field_by_PRISMs(const Range &prisms,const BitFieldId id,int verb = -1) = 0;

  /**
    * \brief set field entities from adjacencies of prisms
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param prisms range contains set prisms
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_PRISMs(const Range &prisms,const std::string& name,int verb = -1) = 0;

  /**
    * \brief set field entities from adjacencies of prisms
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set prisms
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_PRISMs(EntityHandle meshset,const std::string& name,int verb = -1) = 0;

  /**
    * \brief remove entities from field
    * \ingroup mofem_field
    *
    */
  virtual PetscErrorCode remove_ents_from_field_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /**
    * \brief remove entities from field
    * \ingroup mofem_field
    *
    */
  virtual PetscErrorCode remove_ents_from_field(const std::string& name,const EntityHandle meshset,const EntityType type,int verb = -1) = 0;

  /**
    * \brief remove entities from field
    * \ingroup mofem_field
    *
    */
  virtual PetscErrorCode remove_ents_from_field(const std::string& name,const Range &ents,int verb = -1) = 0;

  /**
    * \brief Set order approximation of the entities in the field
    * \ingroup mofem_field
    *
    * \param meshset containing set of the entities (use 0 for all the entities in the meshset)
    * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE, MBVERTEX, see moab documentation
    * \param order approximation order
    */
  virtual PetscErrorCode set_field_order(
    const EntityHandle meshset,
    const EntityType type,
    const std::string& name,
    const ApproximationOrder order,
    int verb = -1
  ) = 0;

  /**
    * \brief Set order approximation of the entities in the field
    * \ingroup mofem_field
    *
    * \param entities
    * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE, MBVERTEX, see moab documentation
    * \param order approximation order
    */
  virtual PetscErrorCode set_field_order(
    const Range &ents,
    const std::string& name,
    const ApproximationOrder order,
    int verb = -1
  ) = 0;

  /**
    * \brief Set order approximation of the entities in the field
    * \ingroup mofem_field
    *
    * \param bit refinement level
    * \param mask bit mask
    * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE, MBVERTEX, see moab documentation
    * \param order approximation order
    */
  virtual PetscErrorCode set_field_order_by_entity_type_and_bit_ref(
    const BitRefLevel &bit,
    const BitRefLevel &mask,
    const EntityType type,
    const std::string& name,
    const ApproximationOrder order,
    int verb = -1
  ) = 0;


  /** \brief list entities in the field
    * \ingroup mofem_field
    */
  virtual PetscErrorCode list_fields() const = 0;

  /** \brief get field meshset
   *
   * \param name of Field
   * Example:\code
   EntityHandle disp_files_meshset = mField.get_field_meshset("DISPLACEMENT");
   * \endcode
   */
  virtual EntityHandle get_field_meshset(const std::string& name) const = 0;

  /**
   * \brief get entities in the field by dimension
   * @param  name field name
   * @param  dim  dim
   * @param  ents ents
   * @return      error code

   * \ingroup mofem_field
   */
  virtual PetscErrorCode get_field_entities_by_dimension(
    const std::string name,int dim,Range &ents
  ) const = 0;

  /**
   * \brief get entities in the field by type
   * @param  name field name
   * @param  type entity type
   * @param  ents ents
   * @return      error code

   * \ingroup mofem_field
   */
  virtual PetscErrorCode get_field_entities_by_type(
    const std::string name,EntityType type,Range &ents
  ) const = 0;

  /**
   * \brief get entities in the field by handle
   * @param  name field name
   * @param  ents ents
   * @return      error code

   * \ingroup mofem_field
   */
  virtual PetscErrorCode get_field_entities_by_handle(
    const std::string name,Range &ents
  ) const = 0;


  /** \brief check if field is in database
   * \ingroup mofem_field
   *
   * \param name field name
   * \return true if field exist
   *
   */
  virtual bool check_field(const std::string& name) const = 0;

  /** \brief get field structure
   * \ingroup mofem_field
   *
   * \param name field name
   * \return const Field*
   *
   */
  virtual const Field* get_field_structure(const std::string& name) = 0;

 /**@}*/

  /** \name Finite elements */

 /**@{*/

  /**
   * \brief Check if finite element is in database
   * @param  name Name of finite element
   * @return      true if element is declared
   */
  virtual bool check_finite_element(const std::string& name) const = 0;

  /**
    * \brief add finite element
    * \ingroup mofem_fe
    * \param name finite element name
    *
    * Example \code
      ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
      ierr = mField.add_finite_element("PLASTIC"); CHKERRQ(ierr);
   \endcode
    */
  virtual PetscErrorCode add_finite_element(const std::string &fe_name,enum MoFEMTypes bh = MF_EXCL) = 0;

  /**
    * \brief modify finite element table, only for advanced user
    * \ingroup mofem_fe
    *
    * Using that functions means that you like to do something not usual.
    *
    */
  virtual PetscErrorCode modify_finite_element_adjacency_table(const std::string &fe_name,const EntityType type,ElementAdjacencyFunct function) = 0;


  /** \brief set finite element field data
   * \ingroup mofem_fe
   *
   * \param name finite element name
   * \param name field name
   *
   * This function will set memory in the form of a vector
   */
  virtual PetscErrorCode modify_finite_element_add_field_data(const std::string &fe_name,const std::string &name_filed) = 0;

  /** \brief unset finite element field data
   * \ingroup mofem_fe
   *
   * \param name finite element name
   * \param name field name
   *
   * This function will set memory in the form of a vector
   */
  virtual PetscErrorCode modify_finite_element_off_field_data(const std::string &fe_name,const std::string &name_filed) = 0;

    /** \brief set field row which finite element use
     * \ingroup mofem_fe
     *
     * \param name finite element name
     * \param name field name
     */
  virtual PetscErrorCode modify_finite_element_add_field_row(const std::string &fe_name,const std::string &name_row) = 0;

    /** \brief unset field row which finite element use
     * \ingroup mofem_fe
     *
     * \param name finite element name
     * \param name field name
     */
  virtual PetscErrorCode modify_finite_element_off_field_row(const std::string &fe_name,const std::string &name_row) = 0;


    /** \brief set field col which finite element use
     *
     * \param name finite element name
     * \param name field name
     */
  virtual PetscErrorCode modify_finite_element_add_field_col(const std::string &fe_name,const std::string &name_row) = 0;

    /** \brief unset field col which finite element use
     * \ingroup mofem_fe
     *
     * \param name finite element name
     * \param name field name
     */
  virtual PetscErrorCode modify_finite_element_off_field_col(const std::string &fe_name,const std::string &name_row) = 0;

  /**
   * \brief add entities to finite element
   * \ingroup mofem_fe
   * @param  entities   meshset or range form were entities taken
   * @param  type      type of entity
   * @param  name      name of field
   * @param  recursive take entities from meshsets in meshset
   * @return           error code
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_type(
    const EntityHandle entities,const EntityType type,const std::string &name,const bool recursive = true
  ) = 0;

  /**
   * \brief add entities to finite element
   * \ingroup mofem_fe
   * @param  entities  meshset or range form were entities taken
   * @param  dim       dimension
   * @param  name      name of field
   * @param  recursive take entities from meshsets in meshset
   * @return           error code
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_dim(
    const EntityHandle entities,const int dim,const std::string &name,const bool recursive = true
  ) = 0;

  /**
   * \brief add entities to finite elements
   * \ingroup mofem_fe
   * @param  ents range of entities
   * @param  type type of entity (MBVERTEX, MBEDGE, MBTRI, ...)
   * @param  name name of finite element
   * @return      error code
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_type(
    const Range& ents,const EntityType type,const std::string &name
  ) = 0;

  /**
   * \brief add entities to finite elements
   * \ingroup mofem_fe
   * @param  ents range of entities
   * @param  dim dimension of entities
   * @param  name name of finite element
   * @return      error code
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_dim(
    const Range& ents,const int dim,const std::string &name
  ) = 0;

  /** \brief add TET entities from given refinement level to finite element database given by name
   * \ingroup mofem_fe
   *
   * \param BitRefLevel bit
   * \param BitRefLevel mask
   * \param finite element name
   * \param finite element type
   * \param verrbose level
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_bit_ref(
    const BitRefLevel &bit,const BitRefLevel &mask,const std::string &name,EntityType type,int verb = -1
  ) = 0;

  /** \brief add EDGES entities from range to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_type
   *
   * \param range contains tetrahedron
   * \param name Finite Element name
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_by_EDGEs(const Range& edge,const std::string &name) = 0;

  /**
   * \brief add EDGES finite elements
   *
   * \deprecated use add_ents_to_finite_element_by_type
   *
   * @param  meshset
   * @param  name      name of finite element
   * @param  recursive take entities from meshsets in meshset
   * @return           error code
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_by_EDGEs(
    const EntityHandle meshset,const std::string &name,const bool recursive = false
  ) = 0;

  /** \brief add VERTICES entities from range to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_type
   *
   * \param range contains tetrahedron
   * \param name Finite Element name
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_by_VERTICEs(const Range& vert,const std::string &name) = 0;

  /** \brief add TRI entities from range to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_type
   *
   * \param range contains tetrahedron
   * \param name Finite Element name
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_by_TRIs(const Range& tris,const std::string &name) = 0;

  /** \brief add TRI entities from meshset to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_type
   *
   * \param range contains tetrahedron
   * \param name Finite Element name
   * \param recursive if true parent meshset is searched recursively
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_by_TRIs(
    const EntityHandle meshset,const std::string &name,const bool recursive = false
  ) = 0;

  /** \brief add TET entities from range to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_type
   *
   * \param range contains tetrahedron
   * \param name Finite Element name
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_by_TETs(const Range& tets,const std::string &name) = 0;

  /** \brief add TET entities from meshset to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_type
   *
   * \param meshset contains tetrahedron
   * \param name Finite Element name
   * \param recursive if true parent meshset is searched recursively
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_by_TETs(
    const EntityHandle meshset,const std::string &name,const bool recursive = false
  ) = 0;

  /** \brief add PRISM entities from meshset to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_type
   *
   * \param range contains tetrahedron
   * \param name Finite Element name
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_by_PRISMs(const Range& prims,const std::string &name) = 0;

  /** \brief add TET entities from meshset to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_type
   *
   * \param meshset contains tetrahedron
   * \param name Finite Element name
   * \param recursive if true parent meshset is searched recursively
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_by_PRISMs(
    const EntityHandle meshset,const std::string &name,const bool recursive = false
  ) = 0;

  /** \brief add TET elements from given refinement level to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_bit_ref with mask explicitly given
   *
   * \param BitRefLevel bit
   * \param finite element name
   * \param finite elenent type
   * \param verrbose level
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(
    const BitRefLevel &bit,const std::string &name,EntityType type,int verb = -1
  ) = 0;

  /** \brief add TET entities from given refinement level to finite element database given by name
   * \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_bit_ref with mask explicitly given
   *
   * \param BitRefLevel bit
   * \param BitRefLevel mask
   * \param finite element name
   * \param finite element type
   * \param verrbose level
   */
  DEPRECATED virtual PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(
    const BitRefLevel &bit,const BitRefLevel &mask,const std::string &name,EntityType type,int verb = -1
  ) = 0;


  /** get finite element meshset
   * \ingroup mofem_fe
   *
   */
  virtual EntityHandle get_finite_element_meshset(const std::string& name) const = 0;

  /**
   * \brief get entities in the finite element by dimension
   * @param  name finite element name
   * @param  dim  dim
   * @param  ents ents
   * @return      error code

   * \ingroup mofem_field
   */
  virtual PetscErrorCode get_finite_element_entities_by_dimension(
    const std::string name,int dim,Range &ents
  ) const = 0;

  /**
   * \brief get entities in the finite element by type
   * @param  name finite element name
   * @param  type entity type
   * @param  ents ents
   * @return      error code

   * \ingroup mofem_field
   */
  virtual PetscErrorCode get_finite_element_entities_by_type(
    const std::string name,EntityType type,Range &ents
  ) const = 0;

  /**
   * \brief get entities in the finite element by handle
   * @param  name finite element name
   * @param  ents ents
   * @return      error code

   * \ingroup mofem_field
   */
  virtual PetscErrorCode get_finite_element_entities_by_handle(
    const std::string name,Range &ents
  ) const = 0;


  /** \brief remove elements from given refinement level to finite element database
   * \ingroup mofem_fe
   *
   * \param BitRefLevel bit
   * \param BitRefLevel mask
   * \param verbose level
   */
  virtual PetscErrorCode remove_ents_from_finite_element_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** \brief remove entities from given refinement level to finite element database
   *
   */
  virtual PetscErrorCode remove_ents_from_finite_element(const std::string &name,const EntityHandle meshset,const EntityType type,int verb = -1) = 0;

  /** \brief remove entities from finite element database
   * \ingroup mofem_fe
   *
   */
  virtual PetscErrorCode remove_ents_from_finite_element(const std::string &name,const Range &ents,int verb = -1) = 0;

  /** \brief add MESHSET element to finite element database given by name
   * \ingroup mofem_fe
   *
   * \param meshset contains all entities that could be used for finite element
   * \param name Finite Element name
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const std::string& name,const bool recursive = false) = 0;

  /** \brief list finite elements in database
   * \ingroup mofem_fe
   */
  virtual PetscErrorCode list_finite_elements() const = 0;

  /// list adjacencies
  virtual PetscErrorCode list_adjacencies() const = 0;

 /**@}*/

  /** \name Problems */

 /**@{*/

  /** \brief Add problem
   * \ingroup mofem_problems
   */
  virtual PetscErrorCode add_problem(const std::string& name,enum MoFEMTypes bh = MF_EXCL,int verb = -1) = 0;

  /**
   * \brief check if problem exist
   * @param  name problem name
   * @return      true if problem is in database
   */
  virtual bool check_problem(const std::string name) = 0;

  /** \brief Delete problem
  * \ingroup mofem_problems
  */
  virtual PetscErrorCode delete_problem(const std::string name) = 0;

  /** \brief add finite element to problem, this add entities assigned to finite element to a particular problem
   * \ingroup mofem_problems
   *
   * \param name Problem name
   * \param name Finite Element name
   */
  virtual PetscErrorCode modify_problem_add_finite_element(const std::string &name_problem,const std::string &fe_name) = 0;

  /** \brief unset finite element from problem, this remove entities assigned to finite element to a particular problem
   * \ingroup mofem_problems
   *
   *  Note: If problem is build, it need to be cleaned to make this effective
   *
   * \param name Problem name
   * \param name Finite Element name
   */
  virtual PetscErrorCode modify_problem_unset_finite_element(const std::string &name_problem,const std::string &fe_name) = 0;


  /** \brief add ref level to problem
   * \ingroup mofem_problems
   *
   * if same finite element is solved using different level of refinements, than the level of refinement has to be specificied to problem in query
   *
   * \param name Problem name
   * \param BitRefLevel bitLevel
   * Example: \code
   ierr = mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF1","ELASTIC"); CHKERRQ(ierr);
   ierr = mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF2","ELASTIC"); CHKERRQ(ierr);

   ierr = mField.modify_problem_ref_level_add_bit("BEAM_BENDING_ON_MESH_REF1",bit_level1); CHKERRQ(ierr);
   ierr = mField.modify_problem_ref_level_add_bit("BEAM_BENDING_ON_MESH_REF2",bit_level2); CHKERRQ(ierr);
   *\endcode
   * Two Problems exist and solved independently, both are elastic, but solved using different mesh refinement <br>
  */
  virtual PetscErrorCode modify_problem_ref_level_add_bit(const std::string &name_problem,const BitRefLevel &bit) = 0;

  /** \brief set ref level for problem
   * \ingroup mofem_problems
   *
   * if same finite element is solved using different level of refinements, than the level of refinement has to be specificied to problem in query
   *
   * \param name Problem name
   * \param BitRefLevel bitLevel
   * Example: \code
   ierr = mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF1","ELASTIC"); CHKERRQ(ierr);
   ierr = mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF2","ELASTIC"); CHKERRQ(ierr);

   ierr = mField.modify_problem_ref_level_set_bit("BEAM_BENDING_ON_MESH_REF1",bit_level1); CHKERRQ(ierr);
   ierr = mField.modify_problem_ref_level_set_bit("BEAM_BENDING_ON_MESH_REF2",bit_level2); CHKERRQ(ierr);
   *\endcode
   * Two Problems exist and solved independently, both are elastic, but solved using different mesh refinement <br>

   \bug Problem bit level should be defined by bit and mask for better flexibility

   */
  virtual PetscErrorCode modify_problem_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit) = 0;

  /** \brief set dof mask ref level for problem
   * \ingroup mofem_problems
   *
   */
  virtual PetscErrorCode modify_problem_mask_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit) = 0;

  /** \brief list problems
   * \ingroup mofem_problems
   */
  virtual PetscErrorCode list_problem() const = 0;

  /** build fields
    * \ingroup mofem_field
   */
  virtual PetscErrorCode build_fields(int verb = -1) = 0;

  /** list dofs
    * \ingroup mofem_field
   */
  virtual PetscErrorCode list_dofs_by_field_name(const std::string &name) const = 0;

 /**@}*/

  /** \name Clear dofs and entities */

 /**@{*/

  /** Clear inactive dofs
    * \ingroup mofem_field
   */
  virtual PetscErrorCode clear_inactive_dofs(int verb = -1) = 0;

  /** Clear dofs by bit level
    * \ingroup mofem_field
   */
  virtual PetscErrorCode clear_dofs_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** Clear ents by bit level
    * \ingroup mofem_field
   */
  virtual PetscErrorCode clear_ents_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** Clear dofs by field name
    * \ingroup mofem_field
   */
  virtual PetscErrorCode clear_dofs_fields(const std::string &name,const Range ents,int verb = -1) = 0;

  /** Clear entities by field name
    * \ingroup mofem_field
   */
  virtual PetscErrorCode clear_ents_fields(const std::string &name,const Range enst,int verb = -1) = 0;

 /**@}*/

  /** \name Build fields, finite elements and problems */

 /**@{*/

  /**
   * \brief Build finite elements
   * \ingroup mofem_fe
   *
   * Build finite element data structures. Have to be run before problem and adjacencies
   * are constructed.
   *
   * @param  verb Verbosity level
   * @return      Error code
   */
  virtual PetscErrorCode build_finite_elements(int verb = -1) = 0;

  /**
   * \brief Build finite elements
   * \ingroup mofem_fe
   *
   * Build finite element data structures. Have to be run before problem and adjacencies
   * are constructed.
   *
   * @param  fe_name  Name of finite element
   * @param  ents_ptr Pointer to range of finite elements
   * @param  verb     Verbosity level
   * @return      Error code
   */
  virtual PetscErrorCode build_finite_elements(const string fe_name,const Range *ents_ptr = NULL,int verb = -1) = 0;

 /**@}*/

  /** \name Clear finite elements */

 /**@{*/

  /** clear finite elements
    */
  virtual PetscErrorCode clear_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** clear finite elements
    */
  virtual PetscErrorCode clear_finite_elements(const std::string &name,const Range &ents,int verb = -1) = 0;

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
  virtual PetscErrorCode build_adjacencies(const Range &ents,int verb = -1) = 0;

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
  virtual PetscErrorCode build_adjacencies(const BitRefLevel &bit,int verb = -1) = 0;

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
  virtual PetscErrorCode build_adjacencies(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

 /**@}*/

  /** \name Clear adjacencies */

 /**@{*/

  /** \brief clear adjacency map for finite elements on given bit level
    *
    * \param bit
    * \param mask
    */
  virtual PetscErrorCode clear_adjacencies_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** \brief clear adjacency map for entities on given bit level
    *
    * \param bit
    * \param mask
    */
  virtual PetscErrorCode clear_adjacencies_entities(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

 /**@}*/

  /** \name Build problems (DEPRECATED SHOULD NOT USE THIS)*/

 /**@{*/

  /** \brief build problem data structures
   *
   * \note If square_matrix is set to true, that indicate that problem is structurally
   * symmetric, i.e. rows and columns have the same dofs and are indexed in the same
   * way.

   \deprecated Use ProblemsManager to build and partition problems

   *
   * @param  problem pointer
   * @param  square_matrix structurally symmetric problem
   * @param  verb          verbosity level
   * @return               error code
   *
   */
  DEPRECATED virtual PetscErrorCode build_problem(const std::string &name,const bool square_matrix,int verb = -1) = 0;

  /** \brief build problem data structures
   *
   * \note If square_matrix is set to true, that indicate that problem is structurally
   * symmetric, i.e. rows and columns have the same dofs and are indexed in the same
   * way.

   \deprecated Use ProblemsManager to build and partition problems

   * @param  name          problem name
   * @param  square_matrix structurally symmetric problem
   * @param  verb          verbosity level
   * @return               error code
   *
   */
  DEPRECATED virtual PetscErrorCode build_problem(Problem *problem_ptr,const bool square_matrix,int verb = -1) = 0;

  /** \brief clear problem
   * \ingroup mofem_problems
   */
  virtual PetscErrorCode clear_problem(const std::string &name,int verb = -1) = 0;

  /** \brief build problem data structures

   \deprecated Use MoFEM::Interface::build_problem(const std::string &name,const
   bool square_matrix,int verb = -1) instead. This function not allows to Control
   if problem is structurally symmetric.

   */
  DEPRECATED virtual PetscErrorCode build_problems(int verb = -1) = 0;

  /** \brief clear problems
   * \ingroup mofem_problems
   */
  virtual PetscErrorCode clear_problems(int verb = -1) = 0;

  /** \brief build problem data structures, assuming that mesh is distributed (collective)

   \deprecated Use ProblemsManager to build and partition problems

   Mesh is distributed, that means that each processor keeps only own part of
   the mesh and shared entities.

   Collective - need to be run on all processors in communicator, i.e. each
   function has to call this function.

   */
  DEPRECATED virtual PetscErrorCode build_problem_on_distributed_mesh(
    const std::string &name,const bool square_matrix = true,int verb = -1
  ) = 0;

  /** \brief build problem data structures, assuming that mesh is distributed (collective)

   \deprecated Use ProblemsManager to build and partition problems

   Mesh is distributed, that means that each processor keeps only own part of
   the mesh and shared entities.

   Collective - need to be run on all processors in communicator, i.e. each
   function has to call this function.

   */
  DEPRECATED virtual PetscErrorCode build_problem_on_distributed_mesh(
    Problem *problem_ptr,
    const bool square_matrix = true,
    int verb = -1
  ) = 0;

  /** \brief build problem data structures, assuming that mesh is distributed (collective)

   \deprecated Use ProblemsManager to build and partition problems

   Mesh is distributed, that means that each processor keeps only own part of
   the mesh and shared entities.

   Collective - need to be run on all processors in communicator, i.e. each
   function has to call this function.

   */
  DEPRECATED virtual PetscErrorCode build_problem_on_distributed_mesh(int verb = -1) = 0;

  /**
   * \brief Set partition tag to each finite element in the problem
   *
   * This will use one of the mesh partitioning programs available from PETSc
   * See <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatPartitioningType.html>
   *
   * @param  ents        Entities to partition
   * @param  dim         Dimension of entities to partition
   * @param  adj_dim     Adjacency dimension
   * @param  n_parts     Number of partitions
   * @param  verb        Verbosity level
   * @return             Error code
   */
  DEPRECATED virtual PetscErrorCode partition_mesh(
    const Range &ents,const int dim,const int adj_dim,const int n_parts,int verb = -1
  ) = 0;

  /** \brief partition problem dofs

   \deprecated Use ProblemsManager to build and partition problems

   *
   * \param name problem name
   */
  DEPRECATED virtual PetscErrorCode partition_simple_problem(const std::string &name,int verb = -1) = 0;

  /** \brief partition problem dofs (collective)

   \deprecated Use ProblemsManager to build and partition problems

   *
   * \param name problem name
   */
  DEPRECATED virtual PetscErrorCode partition_problem(const std::string &name,int verb = -1) = 0;

  /**
    * \brief build indexing and partition problem inheriting indexing and partitioning from two other problems
    *
    \deprecated Use ProblemsManager to build and partition problems
    * \param name problem name
    * \param problem_for_rows problem used to index rows
    * \param copy_rows just copy rows dofs
    * \param problem_for_cols problem used to index cols
    * \param copy_cols just copy cols dofs
    *
    * If copy_rows/copy_cols is set to false only partition is copied between problems.
    *
    */
  DEPRECATED virtual PetscErrorCode partition_compose_problem(
    const std::string &name,
    const std::string &problem_for_rows,
    const bool copy_rows,
    const std::string &problem_for_cols,
    const bool copy_cols,
    int verb = -1
  ) = 0;

  /**
   * \brief build sub problem

   \deprecated Use ProblemsManager to build and partition problems

   * @param  out_name problem
   * @param  fields_row  vector of fields composing problem
   * @param  fields_col  vector of fields composing problem
   * @param  main_problem main problem
   * @return              error code
   */
  DEPRECATED virtual PetscErrorCode build_sub_problem(
    const std::string &out_name,
    const std::vector<std::string> &fields_row,
    const std::vector<std::string> &fields_col,
    const std::string &main_problem,
    const bool square_matrix = true,
    int verb = -1
  ) = 0;

  /** \brief determine ghost nodes
   * \ingroup mofem_field
   * \deprecated Use ProblemsManager to build and partition problems
   *
   * \param name problem name
   */
  DEPRECATED virtual PetscErrorCode partition_ghost_dofs(const std::string &name,int verb = -1) = 0;

  /** \brief partition finite elements
   *
   * Function which partition finite elements based on dofs partitioning.<br>
   * In addition it sets information about local row and cols dofs at given element on partition.
   *
   * \deprecated Use ProblemsManager to build and partition problems
   *
   * \param name problem name
   */
  DEPRECATED virtual PetscErrorCode partition_finite_elements(const std::string &name,
    bool part_from_moab = false,int low_proc = -1,int hi_proc = -1,int verb = -1
  ) = 0;

  /** \brief check if matrix fill in correspond to finite element indices

    This is used to check consistency of code. If problem is notices with
    additional non-zero elements in matrix, this function can help detect problem.
    Should be used as a part of atom tests

    \param problem_nanme
    \param row print info at particular row
    \param col print info at particular col

    */
  virtual PetscErrorCode partition_check_matrix_fill_in(
    const std::string &problem_name,int row,int col,int verb
  ) = 0;

  /**
   * \brief resolve shared entities for finite elements in the problem
   * \ingroup mofem_problems

   * @param  problem_ptr  problem pointer
   * @param  fe_name     finite element name
   * @param  verb        verbosity level
   * @return             error code
   *
   * This allows for tag reduction or tag exchange, f.e.

   \code
   Tag th;
   rval = mField.get_moab().tag_get_handle("ADAPT_ORDER",th); CHKERRQ_MOAB(rval);
   ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
   // rval = pcomm->reduce_tags(th,MPI_SUM,prisms);
   rval = pcomm->exchange_tags(th,prisms);
   \endcode

   *
   */
  virtual PetscErrorCode resolve_shared_ents(const Problem *problem_ptr,const std::string &fe_name,int verb = -1) = 0;

  /**
   * \brief resolve shared entities for finite elements in the problem
   * \ingroup mofem_problems

   * @param  name        problem name
   * @param  fe_name     finite element name
   * @param  verb        verbosity level
   * @return             error code
   *
   * This allows for tag reduction or tag exchange, f.e.

   \code
   ierr = m_field.resolve_shared_ents(problem_ptr,"SHELL_ELEMENT"); CHKERRQ(ierr);
   Tag th;
   rval = mField.get_moab().tag_get_handle("ADAPT_ORDER",th); CHKERRQ_MOAB(rval);
   ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
   // rval = pcomm->reduce_tags(th,MPI_SUM,prisms);
   rval = pcomm->exchange_tags(th,prisms);
   \endcode

   *
   */
  virtual PetscErrorCode resolve_shared_ents(const std::string &name,const std::string &fe_name,int verb = -1) = 0;

  /** \deprecated use ProblemsManager
   * \brief Get layout of elements in the problem
   * \ingroup mofem_problems
   *
   * In layout is stored information how many elements is on each processor, for
   * more information look int petsc documentation
   * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/IS/PetscLayoutCreate.html#PetscLayoutCreate>
   *
   * @param  name    problem name
   * @param  fe_name finite elements
   * @param  layout  layout
   * @param  verb    verbosity level
   * @return         error code
   */
  DEPRECATED virtual PetscErrorCode get_problem_elements_layout(
    const std::string &name,const std::string &fe_name,PetscLayout *layout,int verb = -1
  ) = 0;


  /**
    * \brief add finite elements to the meshset
    * \ingroup mofem_problems
    *
    * \param name is problem name
    * \param fe_name
    * \param meshset
    */
  virtual PetscErrorCode get_problem_finite_elements_entities(const std::string &name,const std::string &fe_name,const EntityHandle meshset) = 0;

 /**@}*/

  /** \name Create vectors */

 /**@{*/

  /** \deprecated use VecManager
   * \brief create local vector for problem
   * \ingroup mofem_vectors
   *
   * \param name problem name
   * \param RowColData specify what data is taken from Row, Col or Data
   * \param Vec the vector where data is stored
   */
  DEPRECATED virtual PetscErrorCode VecCreateSeq(const std::string &name,RowColData rc,Vec *V) const = 0;

  /** \deprecated use VecManager
   * \brief create ghost vector for problem (collective)
   * \ingroup mofem_vectors

  collective - need to be run on all processors in communicator

   * \param name problem name
   * \param RowColData specify what data is taken from Row, Col or Data
   * \param Vec the vector where data is stored
   */
  DEPRECATED virtual PetscErrorCode VecCreateGhost(const std::string &name,RowColData rc,Vec *V) const = 0;

 /**@}*/

  /** \name Create matrices */

 /**@{*/

  /**
    * \brief create Mat (MPIAIJ) for problem (collective)
    *
    * \param name of the problem
    */
  virtual PetscErrorCode MatCreateMPIAIJWithArrays(const std::string &name,Mat *Aij,int verb = -1) = 0;

  /**
   * \brief Create Adj matrix
   * @param  name [description]
   * @param  Adj  [description]
   * @param  verb [description]
   * @return      [description]
   */
  virtual PetscErrorCode MatCreateMPIAdj_with_Idx_mi_tag(const std::string &name,Mat *Adj,int verb = -1) = 0;



  /**
    * \brief create Mat (AIJ) for problem
    *
    * \param name of the problem
    */
  virtual PetscErrorCode MatCreateSeqAIJWithArrays(const std::string &name,Mat *Aij,PetscInt **i,PetscInt **j,PetscScalar **v,int verb = -1) = 0;

 /**@}*/

  /** \name Create IS */

 /**@{*/

  /** \deprecated Use ISManager
    * \brief create IS for give two problems and field

    Note that indices are ordered in ascending order of local indices in problem_y

    \param x_problem name of problem
    \param x_field_name name of field in problem_x
    \param x_rc that is ROW or COL
    \param y_problem name of problem
    \param y_field_name name of field in problem_y
    \param y_rc that is ROW or COL

    \retval idx indexes in problem_x
    \retval idy indexes in problem_y

    */
  DEPRECATED virtual PetscErrorCode ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    std::vector<int> &idx,std::vector<int> &idy,int verb = -1
  ) const = 0;

  /** \deprecated Use ISManager
    * \brief create IS for give two problems and field

    Indices are sorted by global PETSc index in problem_x.

    \param x_problem name of problem
    \param x_field_name name of field in problem_x
    \param x_rc that is ROW or COL
    \param y_problem name of problem
    \param y_field_name name of field in problem_y
    \param y_rc that is ROW or COL

    \retval ix IS indexes in problem_x (can be PETSC_NULL)
    \retval iy IS indexes in problem_y

    */
  DEPRECATED virtual PetscErrorCode ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    IS *ix,IS *iy,int verb = -1
  ) const = 0;

  /** \deprecated Use ISManager
    * \brief create IS for given order range (collective)

    * \param problem name
    * \param rc ROW or COL
    * \param min_order
    * \param max_order
    * \retval is out value

    */
  DEPRECATED virtual PetscErrorCode ISCreateProblemOrder(
    const std::string &problem,RowColData rc,int min_order,int max_order,IS *is,int verb = -1
  ) const = 0;

  /** \deprecated Use ISManager
    * \brief create IS for given problem, field and rank range (collective)
    * \ingroup mofem_vectors

    * \param problem name
    * \param rc ROW or COL
    * \param field name
    * \param min_coeff_idx
    * \param max_coeff_idx
    * \retval is out value

    */
  DEPRECATED virtual PetscErrorCode ISCreateProblemFieldAndRank(
    const std::string &problem,
    RowColData rc,
    const std::string &field,
    int min_coeff_idx,
    int max_coeff_idx,
    IS *is,
    int verb = -1
  ) const = 0;

 /**@}*/

  /** \name Scatter vectors */

 /**@{*/

  /** \deprecated use VecManager
    * \brief create scatter for vectors form one to another problem (collective)
    * \ingroup mofem_vectors
    *
    * User specify what name of field on one problem is scattered to another.
    *
    * \ingroup mofem_vectors
    *
    * \param xin vector
    * \param x_proble problem name
    * \param x_field name
    * \param yin vector
    * \param y_problem problem name
    * \param y_field_name
    * \retval newctx scatter

    */
  DEPRECATED virtual PetscErrorCode VecScatterCreate(
    Vec xin,
    const std::string &x_problem,
    const std::string &x_field_name,
    RowColData x_rc,
    Vec yin,
    const std::string &y_problem,
    const std::string &y_field_name,
    RowColData y_rc,
    VecScatter *newctx,
    int verb = -1
  ) const = 0;

  /** \deprecated use VecManager
    * \brief create scatter for vectors form one to another problem (collective)
    * \ingroup mofem_vectors
    *
    * \param xin vector
    * \param x_proble problem name
    * \param yin vector
    * \param y_problem problem name
    * \retval newctx scatter

    */
  DEPRECATED virtual PetscErrorCode VecScatterCreate(
    Vec xin,const std::string &x_problem,
    RowColData x_rc,
    Vec yin,
    const std::string &y_problem,
    RowColData y_rc,
    VecScatter *newctx,
    int verb = -1
  ) const = 0;

 /**@}*/

  /** \name Set vector and mesh values */

 /**@{*/

  /** \deprecated use VecManager
    * \brief set values of vector from/to meshdatabase
    * \ingroup mofem_vectors
    *
    * \param pointer to problem struture
    * \param RowColData for row or column:e (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities from V vector.
    *
    * SCATTER_FORWARD set vector V from data field entities
    *
    */
  DEPRECATED virtual PetscErrorCode set_local_ghost_vector(
    const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const = 0;

  /** \deprecated use VecManager
    * \brief set values of vector from/to meshdatabase
    * \ingroup mofem_vectors
    *
    * \param name of the problem
    * \param RowColData for row or column:e (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities from V vector.
    *
    * SCATTER_FORWARD set vector V from data field entities
    *
    */
  DEPRECATED virtual PetscErrorCode set_local_ghost_vector(
    const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const = 0;

  /** \deprecated use VecManager
    * \brief set values of vector from/to mesh database (collective)
    * \ingroup mofem_vectors

    collective - need tu be run on all processors in communicator

    * \param pointer to porblem struture
    * \param RowColData for row or column (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  DEPRECATED virtual PetscErrorCode set_global_ghost_vector(
    const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const = 0;

  /** \deprecated use VecManager
    * \brief set values of vector from/to mesh database (collective)
    * \ingroup mofem_vectors

    collective - need tu be run on all processors in communicator

    * \param name of the problem
    * \param RowColData for row or column (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  DEPRECATED virtual PetscErrorCode set_global_ghost_vector(
    const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const = 0;

  /** \deprecated use VecManager
    * \brief Copy vector to field which is not part of the problem
    * \ingroup mofem_vectors
    *
    * \param pointer to poroblem multi_index
    * \param field_name field name used for indexing petsc vectors used in the problem
    * \param cpy_field field name where data from vector are stored
    * \param RowColData for row or column
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  DEPRECATED virtual PetscErrorCode set_other_local_ghost_vector(
    const Problem *problem_ptr,
    const std::string& fiel_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb = -1
  ) = 0;

  /** \deprecated use VecManager
    * \brief Copy vector to field which is not part of the problem
    * \ingroup mofem_vectors
    *
    * \param name problem name
    * \param field_name field name used for indexing petsc vectors used in the problem
    * \param cpy_field field name where data from vector are stored
    * \param RowColData for row or column
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  DEPRECATED virtual PetscErrorCode set_other_local_ghost_vector(
    const std::string &name,
    const std::string& field_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb = -1
  ) = 0;

  /** \deprecated use VecManager
    * \brief Copy vector to field which is not part of the problem (collective)
    * \ingroup mofem_vectors

    collective - need tu be run on all processors in communicator

    * \param problem_ptr pointer to problem
    * \param field_name field name used for indexing petsc vectors used in the problem
    * \param cpy_field field name where data from vector are stored
    * \param RowColData for row or column
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  DEPRECATED virtual PetscErrorCode set_other_global_ghost_vector(
    const Problem *problem_ptr,
    const std::string& field_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb = -1
  ) = 0;

  /** \deprecated use VecManager
    * \brief Copy vector to field which is not part of the problem (collective)
    * \ingroup mofem_vectors

    collective - need tu be run on all processors in communicator

    * \param name problem name
    * \param field_name field name used for indexing petsc vectors used in the problem
    * \param cpy_field field name where data from vector are stored
    * \param RowColData for row or column
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  DEPRECATED virtual PetscErrorCode set_other_global_ghost_vector(
    const std::string &name,
    const std::string& field_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb = -1
  ) = 0;

 /**@}*/

  /** \name Field algebra */

 /**@{*/

  /** \deprecated use FieldBlas
    * \brief axpy fields
    * \ingroup mofem_field_algebra
    * \todo should be moved to independent interface, i.e. FieldAlgebra
    *
    * field_y = field_y + alpha*field_x
    *
    * \param alpha
    * \param field_name_x name of field_x
    * \param field_name_y name of field_y
    * \param error_if_missing throw error if entity/dof exist in field_x but not on field_y
    * \param create_if_missing creat dof in field_y from fiedl_x if it is not database
    *
    */
  DEPRECATED virtual PetscErrorCode field_axpy(const double alpha,const std::string& fiel_name_x,const std::string& field_name_y,bool error_if_missing = false,bool creat_if_missing = false) = 0;

  /** \deprecated use FieldBlas
    * \brief scale field
    * \ingroup mofem_field_algebra
    * \todo should be moved to independent interface, i.e. FieldAlgebra
    *
    * \param alpha is a scaling factor
    * \field_name  is a field name
    *
    */
  DEPRECATED virtual PetscErrorCode field_scale(const double alpha,const std::string& field_name) = 0;

  /** \brief use FieldBlas
    * \brief set field
    * \ingroup mofem_field_algebra
    * \todo should be moved to independent interface, i.e. FieldAlgebra
    *
    * field_y = val
    *
    * \param val
    * \param entity type
    * \param field_name
    *
    */
  DEPRECATED virtual PetscErrorCode set_field(const double val,const EntityType type,const std::string& field_name) = 0;

  /** \deprecated use FieldBlas
    * \brief set field
    * \ingroup mofem_field_algebra
    * \todo should be moved to independent interface, i.e. FieldAlgebra
    *
    * field_y = val
    *
    * \param val
    * \param entity type
    * \param on enties
    * \param field_name
    *
    */
  DEPRECATED virtual PetscErrorCode set_field(const double val,const EntityType type,const Range &ents,const std::string& field_name) = 0;

 /**@}*/

  /** \name Making loops on elements and entities */

 /**@{*/

  /** \brief Set data for BasicMethod
    *
    * This function set data about problem, adjacencies and other MultIindices
    * in database. This function can be used a special case when user need to
    * do some pre- and post-processing before matrix or vector is initiated, or
    * to assemble matrix for group of FEMethods. Is used by calsses classes
    * SnesCtx and TsCtx. Look for more details there.
    *
    * FIXME: Here we need example
    *
    * \param pointer to problem data structure
    * \param method user method derived from BasicMethod
    *
  **/
  virtual PetscErrorCode problem_basic_method_preProcess(const Problem *problem_ptr,BasicMethod &method,int verb = -1) = 0;

  /** \brief Set data for BasicMethod
    *
    * This function set data about problem, adjacencies and other MultIindices
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
  virtual PetscErrorCode problem_basic_method_preProcess(const std::string &problem_name,BasicMethod &method,int verb = -1) = 0;

  /** \brief Set data for BasicMethod
    * \ingroup mofem_loops
    *
    * This function set data about problem, adjacencies and other MultIindices
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
  virtual PetscErrorCode problem_basic_method_postProcess(const Problem *problem_ptr,BasicMethod &method,int verb = -1) = 0;

  /** \brief Set data for BasicMethod
    * \ingroup mofem_loops
    *
    * This function set data about problem, adjacencies and other MultIindices
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
  virtual PetscErrorCode problem_basic_method_postProcess(const std::string &problem_name,BasicMethod &method,int verb = -1) = 0;

  /** \brief Make a loop over finite elements.
   *
   * This function is like swiss knife, is can be used to post-processing or matrix
   * and vectors assembly. It makes loop over given finite element for given
   * problem. The particular methods executed on each element are given by
   * class derived form Interface::FEMethod. At beginning of each loop user defined
   * function (method)  preProcess() is called, for each element operator() is
   * executed, at the end loop finalizes with user defined function (method)
   * postProcess().
   *
   * Methods are executed only for local elements at given processor.
   *
   * For more details pleas look to examples.

   * @param  problem_name problem consisting set of elements
   * @param  fe_name      name of element in problem
   * @param  method       class derived form Interface::FEMethod
   * @param  bh           if bH = MF_EXIST, throws error if fe_name does not exist
   * @param  verb         verbosity level
   * @return              error code

   * \ingroup mofem_loops
  **/
  virtual PetscErrorCode loop_finite_elements(
    const std::string &problem_name,const std::string &fe_name,FEMethod &method,
    MoFEMTypes bh = MF_EXIST,int verb = -1
  ) = 0;

  /** \brief Make a loop over finite elements on partitions from upper to lower rank.
   *
   * This function is like swiss knife, is can be used to post-processing or matrix
   * and vectors assembly. It makes loop over given finite element for given
   * problem. The particular methods executed on each element are given by
   * class derived form Interface::FEMethod. At beginning of each loop user defined
   * function (method)  preProcess() is called, for each element operator() is
   * executed, at the end loop finalizes with user defined function (method)
   * postProcess().
   *
   * For more details please look to examples.
   *
   * Interface::FEMethod

   * @param  problem_ptr pointer to problem consisting set of elements
   * @param  fe_name      name of element in problem
   * @param  method       class derived form Interface::FEMethod
   * @param  lower_rank   lower rank of process owned by finite element
   * @param  upper_rank   lower rank of process owned by finite element
   * @param  bh           if bH = MF_EXIST, throws error if fe_name does not exist
   * @param  verb         verbosity level
   * @return              error code

   * \ingroup mofem_loops
  **/
  virtual PetscErrorCode loop_finite_elements(
    const Problem *problem_ptr,const std::string &fe_name,FEMethod &method,
    int lower_rank,int upper_rank,MoFEMTypes bh = MF_EXIST,int verb = -1
  ) = 0;

  /** \brief Make a loop over finite elements on partitions from upper to lower rank.
   *
   * This function is like swiss knife, is can be used to post-processing or matrix
   * and vectors assembly. It makes loop over given finite element for given
   * problem. The particular methods executed on each element are given by
   * class derived form Interface::FEMethod. At beginning of each loop user defined
   * function (method)  preProcess() is called, for each element operator() is
   * executed, at the end loop finalizes with user defined function (method)
   * postProcess().
   *
   * For more details please look to examples.
   *

   * @param  problem_name pointer to problem consisting set of elements
   * @param  fe_name      name of element in problem
   * @param  method       class derived form Interface::FEMethod
   * @param  lower_rank   lower rank of process owned by finite element
   * @param  upper_rank   lower rank of process owned by finite element
   * @param  bh           if bH = MF_EXIST, throws error if fe_name does not exist
   * @param  verb         verbosity level
   * @return              error code

   * \ingroup mofem_loops
  **/
  virtual PetscErrorCode loop_finite_elements(
    const std::string &problem_name,const std::string &fe_name,FEMethod &method,
    int lower_rank,int upper_rank,MoFEMTypes bh = MF_EXIST,int verb = -1
  ) = 0;

  /** \brief Make a loop over entities

    * \ingroup mofem_loops
    */
  virtual PetscErrorCode loop_dofs(
    const Problem *problem_ptr,const std::string &field_name,RowColData rc,
    EntMethod &method,int lower_rank,int upper_rank,int verb = -1
  ) = 0;

  /** \brief Make a loop over entities

    * \ingroup mofem_loops
    */
  virtual PetscErrorCode loop_dofs(
    const std::string &problem_name,const std::string &field_name,RowColData rc,
    EntMethod &method,int lower_rank,int upper_rank,int verb = -1
  ) = 0;


  /** \brief Make a loop over entities

    * \ingroup mofem_loops
    */
  virtual PetscErrorCode loop_dofs(
    const std::string &problem_name,const std::string &field_name,RowColData rc,
    EntMethod &method,int verb = -1
  ) = 0;

  /** \brief Make a loop over entities

    * \ingroup mofem_field
    */
  virtual PetscErrorCode loop_dofs(const std::string &field_name,EntMethod &method,int verb = -1) = 0;

 /**@}*/

  /** \name Get pointers to multi-index databses */

 /**@{*/

  /** \brief Get fields multi-index from database
    * \ingroup mofem_access
    */
  virtual PetscErrorCode get_fields(const Field_multiIndex **fields_ptr) const = 0;

  /** \brief Get ref entities multi-index from database
    * \ingroup mofem_access
    */
  virtual PetscErrorCode get_ref_ents(const RefEntity_multiIndex **refined_ents_ptr) const = 0;

  /** \brief Get ref finite elements multi-index form database
    * \ingroup mofem_access
    */
  virtual PetscErrorCode get_ref_finite_elements(const RefElement_multiIndex **refined_finite_elements_ptr) const = 0;

  /** \brief Get finite elements multi-index
    * \ingroup mofem_access
    */
  virtual PetscErrorCode get_finite_elements(const FiniteElement_multiIndex **fe_ptr) const = 0;

  /** \brief Get entities finite elements multi-index
    * \ingroup mofem_access
    */
  virtual PetscErrorCode get_ents_finite_elements(const EntFiniteElement_multiIndex **fe_ent_ptr) const = 0;

  /** \brief Get problem database (data structure)
    * \ingroup mofem_access
    */
  virtual PetscErrorCode get_problem(const std::string &problem_name,const Problem **problem_ptr) const = 0;

  /**
   * \brief Get pointer to problems multi-index
   * \ingroup mofem_access
   */
  virtual PetscErrorCode get_problems(const Problem_multiIndex **problems_ptr) const = 0;

  /** \brief Get field multi index
    *
    * \ingroup mofem_access
    *
    */
  virtual PetscErrorCode get_field_ents(const FieldEntity_multiIndex **field_ents) const = 0;

  /** \brief Get dofs multi index
    *
    * \ingroup mofem_access
    *
    */
  virtual PetscErrorCode get_dofs(const DofEntity_multiIndex **dofs_ptr) const = 0;

  /**
    * \brief get begin iterator of filed ents of given name (instead you can use _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    * \ingroup mofem_field
    *
    * for(_IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name
    */
  virtual FieldEntityByFieldName::iterator get_ent_moabfield_by_name_begin(const std::string &field_name) const = 0;

  /**
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    * \ingroup mofem_field
    *
    * for(_IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name
    */
  virtual FieldEntityByFieldName::iterator get_ent_moabfield_by_name_end(const std::string &field_name) const = 0;

  /** \brief loop over all dofs from a moFEM field and particular field
    * \ingroup mofem_field
    */
  #define _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    FieldEntityByFieldName::iterator IT = MFIELD.get_ent_moabfield_by_name_begin(NAME); \
      IT != MFIELD.get_ent_moabfield_by_name_end(NAME); IT++

  /**
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    * \ingroup mofem_field
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name
    */
  virtual DofEntityByFieldName::iterator get_dofs_by_name_begin(const std::string &field_name) const = 0;

  /**
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    * \ingroup mofem_field
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name
    */
  virtual DofEntityByFieldName::iterator get_dofs_by_name_end(const std::string &field_name) const = 0;

  /** loop over all dofs from a moFEM field and particular field
    * \ingroup mofem_field
    */
  #define _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    DofEntityByFieldName::iterator IT = MFIELD.get_dofs_by_name_begin(NAME); \
      IT != MFIELD.get_dofs_by_name_end(NAME); IT++

  /**
    * \brief get begin iterator of filed dofs of given name and ent(instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,ENT,IT)
    * \ingroup mofem_field
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(MFIELD,NAME,ENT,IT)) {
    * 	...
    * }
    *
    * \param field_name
    */
  virtual DofEntityByNameAndEnt::iterator
  get_dofs_by_name_and_ent_begin(const std::string &field_name,const EntityHandle ent) const = 0;

  /**
    * \brief get begin iterator of filed dofs of given name and ent (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,ENT,IT)
    * \ingroup mofem_field
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(MFIELD,NAME,ENT,IT)) {
    * 	...
    * }
    *
    * \param field_name
    */
  virtual DofEntityByNameAndEnt::iterator
  get_dofs_by_name_and_ent_end(const std::string &field_name,const EntityHandle ent) const = 0;

  /** \brief loop over all dofs from a moFEM field and particular field
    * \ingroup mofem_access
    */
  #define _IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(MFIELD,NAME,ENT,IT) \
    DofEntityByNameAndEnt::iterator IT = MFIELD.get_dofs_by_name_and_ent_begin(NAME,ENT); \
      IT != MFIELD.get_dofs_by_name_and_ent_end(NAME,ENT); IT++

  /**
    * \brief get field data from entity and field
    * \ingroup mofem_field
    *
    * this funciont is not recommended to be used in finite elemeny implementation
    *
    */
  template <typename DIT>
  PetscErrorCode get_field_dof_data(
    const std::string& name,const EntityHandle *ent,const int num_ents,DIT dit,int *count = NULL
  ) {
    PetscFunctionBegin;
    if(count!=NULL) *count = 0;
    for(int nn = 0;nn<num_ents;nn++) {
      for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_((*this),name,ent[nn],it)) {
        *(dit++) = (*it)->getFieldData();
        if(count!=NULL) (*count)++;
      }
    }
    PetscFunctionReturn(0);
  }

  /**
    * \brief get field data from entity and field
    * \ingroup mofem_field
    *
    * this function is not recommended to be used in finite element implementation
    *
    */
  template <typename DIT>
  PetscErrorCode get_field_dof_data(
    const std::string& name,const Range &ents,DIT dit,int *count = NULL
  ) {
    PetscFunctionBegin;
    if(count!=NULL) *count = 0;
    for(Range::const_iterator eit = ents.begin();eit!=ents.end();eit++) {
      for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_((*this),name,*eit,it)) {
        *(dit++) = (*it)->getFieldData();
        if(count!=NULL) (*count)++;
      }
    }
    PetscFunctionReturn(0);
  }

  /**
    * \brief get begin iterator of filed dofs of given name and ent type (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,TYPE,IT)
    * \ingroup mofem_field
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(MFIELD,NAME,TYPE,IT)) {
    * 	...
    * }
    *
    * \param field_name
    */
  virtual DofEntityByNameAndType::iterator
  get_dofs_by_name_and_type_begin(const std::string &field_name,const EntityType type) const = 0;

  /**
    * \brief get begin iterator of filed dofs of given name end ent type(instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,TYPE,IT)
    * \ingroup mofem_field
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(MFIELD,NAME,TYPE,IT)) {
    * 	...
    * }
    *
    * \param field_name
    */
  virtual DofEntityByNameAndType::iterator
  get_dofs_by_name_and_type_end(const std::string &field_name,const EntityType type) const = 0;

  /** \brief loop over all dofs from a moFEM field and particular field
    * \ingroup mofem_field
    */
  #define _IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(MFIELD,NAME,TYPE,IT) \
    DofEntityByNameAndType::iterator IT = MFIELD.get_dofs_by_name_and_type_begin(NAME,TYPE); \
      IT != MFIELD.get_dofs_by_name_and_type_end(NAME,TYPE); IT++

  /**
    * \brief get begin iterator of finite elements of given name (instead you can use _IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    * \ingroup mofem_access
    *
    * for(_IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param fe_name
    */
  virtual EntFiniteElementbyName::iterator
  get_fe_by_name_begin(const std::string &fe_name) const = 0;

  /**
    * \brief get end iterator of finite elements of given name (instead you can use _IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    * \ingroup mofem_access
    *
    * for(_IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param fe_name
    */
  virtual EntFiniteElementbyName::iterator
  get_fe_by_name_end(const std::string &fe_name) const = 0;

  /** \brief loop over all finite elements from a moFEM field and FE
    * \ingroup mofem_access
    */
  #define _IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    EntFiniteElementbyName::iterator IT = MFIELD.get_fe_by_name_begin(NAME); \
      IT != MFIELD.get_fe_by_name_end(NAME); IT++


};

}

#endif // __INTERFACE_HPP__

/***************************************************************************//**
 * \defgroup mofem_field Fields
 * \brief Data structure for adding and managing fields
 *
 * \ingroup mofem
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_ref_ents Get entities and adjacencies
 * \brief Get adjacencies/entities for given BitRefLevel (mesh refinement)
 *
 * \ingroup mofem
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_fe Finite elements
 * \brief Adding and managing finite elements
 *
 * \ingroup mofem
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_problems Problems
 * \brief Adding and managing problems
 *
 * \ingroup mofem
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_loops Loops
 * \brief Manages complexities for integrating over finite elements and dofs.
 *
 * \ingroup mofem
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_access Pointers to multi-indices
 * \brief Get direct access to multi-indices in database
 *
 * \ingroup mofem
 ******************************************************************************/
