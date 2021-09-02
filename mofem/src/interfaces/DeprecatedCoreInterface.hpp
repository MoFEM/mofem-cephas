/** \file DeprecatedCoreInterface.hpp
 * \brief Deprecated interface functions
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

#ifndef __INTERFACE_DEPRECATED_HPP__
#define __INTERFACE_DEPRECATED_HPP__

/** \brief name space of MoFEM library functions and classes
 */
namespace MoFEM {

/**
 * \brief Deprecated interface functions
 * \nosubgrouping
 */
struct DeprecatedCoreInterface : public CoreInterface {

  /** \name Interfaces */

  /**@}*/

  /** \name Seed entities */

  /**@{*/

  /** \deprecated use BitRefManager
   * \brief seed 2D entities (Triangles entities only) in the meshset and their
   * adjacencies (only TRIs adjacencies) in a particular BitRefLevel \todo
   * Should be outsourced to separate interface, i.e. BitLevelManager
   *
   * \param EntityHandle MeshSet
   * \param BitRefLevel bitLevel
   *
   */
  DEPRECATED virtual MoFEMErrorCode
  seed_ref_level_2D(const EntityHandle meshset, const BitRefLevel &bit,
                    int verb = -1);

  /** \deprecated use BitRefManager
  * \brief seed 2D entities in the meshset and their adjacencies (only TETs
  adjacencies) in a particular BitRefLevel
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
  DEPRECATED MoFEMErrorCode seed_ref_level_3D(const EntityHandle meshset,
                                              const BitRefLevel &bit,
                                              int verb = -1);

  /** \deprecated use BitRefManager
   * \brief seed entities in the range and their adjacencies in a particular
   * BitRefLevel \todo Should be outsourced to separate interface, i.e.
   * BitLevelManager
   */
  DEPRECATED MoFEMErrorCode seed_ref_level(const Range &ents,
                                           const BitRefLevel &bit,
                                           const bool only_tets = true,
                                           int verb = -1);

  /** \deprecated use BitRefManager
   * brief seed ref level by MESHSET that contains entities other than volumes
   *
   * \param EntityHandle MeshSet
   * \param BitRefLevel bitLevel
   */
  DEPRECATED MoFEMErrorCode seed_ref_level_MESHSET(const EntityHandle meshset,
                                                   const BitRefLevel &bit,
                                                   int verb = -1);

  /** \deprecated use BitRefManager and setElementsBitRefLevel
   * Create finite elements based from entities in meshsets. Throw error if
   * entity is not in database \todo Should be outsourced to separate
   * interface, i.e. BitLevelManager
   *
   * \param EntityHandle meshset
   *
   */
  DEPRECATED MoFEMErrorCode seed_finite_elements(const EntityHandle meshset,
                                                 int verb = -1);

  /** \deprecated use BitRefManager and setElementsBitRefLevel
   * Create finite elements based from entities in meshsets. Throw error if
   * entity is not in database \todo Should be outsourced to separate
   * interface, i.e. BitLevelManager
   *
   * \param Range entities
   *
   */
  DEPRECATED MoFEMErrorCode seed_finite_elements(const Range &entities,
                                                 int verb = -1);

  /**@}*/

  /** \name Getting entities by BitRefLevel */

  /**@{*/

  /**\brief add all ents from ref level given by bit to meshset
   * \deprecated Use MoFEM::BitRefManager interface instead
   * \ingroup mofem_ref_ents
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityType type of entities
   * \retval EntityHandle meshset
   *
   */
  DEPRECATED MoFEMErrorCode get_entities_by_type_and_ref_level(
      const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
      const EntityHandle meshset, int verb = -1);

  /**\brief add all ents from ref level given by bit to meshset
   * \deprecated Use MoFEM::BitRefManager interface instead
   * \ingroup mofem_ref_ents
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityType type of entities
   * \retval ents
   *
   */
  DEPRECATED MoFEMErrorCode get_entities_by_type_and_ref_level(
      const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
      Range &ents, int verb = -1);

  /**\brief add all ents from ref level given by bit to meshset
   * \deprecated Use MoFEM::BitRefManager interface instead
   * \ingroup mofem_ref_ents
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityHandle meshset
   *
   */
  DEPRECATED MoFEMErrorCode
  get_entities_by_ref_level(const BitRefLevel &bit, const BitRefLevel &mask,
                            const EntityHandle meshset);

  /**\brief add all ents from ref level given by bit to meshset
   * \ingroup mofem_ref_ents
   * \deprecated Use MoFEM::BitRefManager interface instead
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \retval ents
   */
  DEPRECATED MoFEMErrorCode get_entities_by_ref_level(const BitRefLevel &bit,
                                                      const BitRefLevel &mask,
                                                      Range &ents);

  /**@}*/

  /**@}*/

  /** \name Mange meshest */

  /**@{*/

  /**
  * \brief check for CUBIT Id and CUBIT type

  \deprecated use MeshsetsManager
  \todo All cubit interface functions should be outsourced to dedicated
  interface

  * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
  * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
  */
  DEPRECATED bool check_msId_meshset(const int msId,
                                     const CubitBCType cubit_bc_type);

  /**
  * \brief add cubit meshset

  \deprecated use MeshsetsManager

  * \param see CubitBC (NODESET, SIDESET or BLOCKSET and more)
  * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET
  * \param name of set

  */
  DEPRECATED MoFEMErrorCode add_cubit_msId(const CubitBCType cubit_bc_tyep,
                                           const int msId,
                                           const std::string name = "");

  /**
  * \brief set attributes to cubit meshset

  \deprecated use MeshsetsManager

  * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET, NODESET,
  SIDESET
  * @param  ms_id         id of meshset
  * @param  attributes    attributes
  * @return               error code
  */
  DEPRECATED MoFEMErrorCode set_cubit_msId_attribites(
      const CubitBCType cubit_bc_type, const int ms_id,
      const std::vector<double> &attributes, const std::string name = "");

  /**
  * \brief set (material) data structure to cubit meshset

  \deprecated use MeshsetsManager
  \todo All cubit interface functions should be outsourced to dedicated
  interface

  * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET, NODESET,
  SIDESET
  * @param  ms_id         id of meshset
  * @param  attributes    attributes
  * @return               error code
  */
  DEPRECATED MoFEMErrorCode set_cubit_msId_attribites_data_structure(
      const CubitBCType cubit_bc_type, const int ms_id,
      const GenericAttributeData &data, const std::string name = "");

  /**
  * \brief set boundary data structure to meshset

  \deprecated use MeshsetsManager
  \todo All cubit interface functions should be outsourced to dedicated
  interface

  * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET, NODESET,
  SIDESET
  * @param  ms_id         id of meshset
  * @param  data          data structure
  * @return               error code
  */
  DEPRECATED MoFEMErrorCode set_cubit_msId_bc_data_structure(
      const CubitBCType cubit_bc_type, const int ms_id,
      const GenericCubitBcData &data);

  /**
  * \brief delete cubit meshset
  * \ingroup mopfem_bc

  \deprecated use MeshsetsManager
  \todo All cubit interface functions should be outsourced to dedicated
  interface

  * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
  * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
  *
  */
  DEPRECATED MoFEMErrorCode delete_cubit_msId(const CubitBCType cubit_bc_type,
                                              const int msId);

  /**
  * \brief get cubit meshset

  \deprecated use MeshsetsManager
  \todo All cubit interface functions should be outsourced to dedicated
  interface

  */
  DEPRECATED MoFEMErrorCode
  get_cubit_msId(const int msId, const CubitBCType cubit_bc_type,
                 const CubitMeshSets **cubit_meshset_ptr);

  DEPRECATED MoFEMErrorCode get_cubit_msId_entities_by_dimension(
      const int ms_id, const CubitBCType cubit_bc_type, const int dimension,
      Range &entities, const bool recursive = false);
  DEPRECATED MoFEMErrorCode get_cubit_msId_entities_by_dimension(
      const int ms_id, const CubitBCType cubit_bc_type, Range &entities,
      const bool recursive = false);
  // DEPRECATED MoFEMErrorCode get_cubit_msId_entities_by_dimension(
  //   const int ms_id,const unsigned int cubit_bc_type, const int
  //   dimension,Range &entities,const bool recursive = false
  // );
  // DEPRECATED MoFEMErrorCode get_cubit_msId_entities_by_dimension(
  //   const int ms_id,const unsigned int cubit_bc_type, Range &entities,const
  //   bool recursive = false
  // );

  /**
  * \brief get entities from CUBIT/meshset of a particular entity dimension

  * Nodeset can contain nodes, edges, triangles and tets. This applies to other
  meshsets too.
  * The nodeset's meshset contain the nodes in the MIDDLE of the surface or
  volume which is done by default in Cubit,
  * Hence if all nodes on a particular nodeset are required,
  * one should get all triangles or tetrahedron for which the nodeset was create
  in Cubit,
  * and get all the connectivities of tris/tets.

  \deprecated use MeshsetsManager
  \todo All cubit interface functions should be outsourced to dedicated
  interface

  \code
  MeshsetsManager *meshset_manager_ptr;
  ierr = m_field.getInterface(meshset_manager_ptr); CHKERRG(ierr);
  ierr =
  meshset_manager_ptr->getEntitiesByDimension(ms_id,cubit_bc_type,dimension,entities,true);
  CHKERRG(ierr); \endcode

  * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
  * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
  * \param dimensions (0 - Nodes, 1 - Edges, 2 - Faces, 3 - Volume(tetrahedral))
  * \param Range containing the retreived entities
  * \param recursive If true, meshsets containing meshsets are queried
  recursively. Returns the contents of meshsets, but not the meshsets themselves
  if true.
  */
  DEPRECATED MoFEMErrorCode get_cubit_msId_entities_by_dimension(
      const int msId, const unsigned int cubit_bc_type, const int dimension,
      Range &entities, const bool recursive = false);

  /**
  * \brief get entities related to CUBIT/meshset,

  * NODESET will get Vertices only, even if the NODESET contains edges, tris and
  tets
  * SIDESET will get Tris, BLOCKSET will get Tets, DISPLACEMENTSET and FORCESET
  are stored in NODESET, PRESSURESET is stored in Sideset.

  \deprecated use MeshsetsManager
  \todo All cubit interface functions should be outsourced to dedicated
  interface

  * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
  * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
  * \param Range containing the retreived entities related to the
  * \param recursive If true, meshsets containing meshsets are queried
  recursively.  Returns the contents of meshsets, but not the meshsets
  themselves if true.
  */
  DEPRECATED MoFEMErrorCode get_cubit_msId_entities_by_dimension(
      const int msId, const unsigned int cubit_bc_type, Range &entities,
      const bool recursive = false);

  /**
  * \brief get meshset from CUBIT Id and CUBIT type

  \deprecated use MeshsetsManager
  \todo All cubit interface functions should be outsourced to dedicated
  interface

  * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
  * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
  * \param meshset where to store the retrieved entities
  */
  DEPRECATED MoFEMErrorCode get_cubit_msId_meshset(
      const int msId, const unsigned int cubit_bc_type, EntityHandle &meshset);

  /**
  * \brief get all CUBIT meshsets by CUBIT type

  \deprecated use MeshsetsManager
  \todo All cubit interface functions should be outsourced to dedicated
  interface

  * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more).
  * \param meshsets is range of meshsets
  */
  DEPRECATED MoFEMErrorCode get_cubit_meshsets(const unsigned int cubit_bc_type,
                                               Range &meshsets);

  /** \deprecated use MeshsetsManager interface instead
   */
  DEPRECATED MoFEMErrorCode print_cubit_displacement_set() const;

  /** \deprecated use MeshsetsManager interface instead
   */
  DEPRECATED MoFEMErrorCode print_cubit_pressure_set() const;

  /** \deprecated use MeshsetsManager interface instead
   */
  DEPRECATED MoFEMErrorCode print_cubit_force_set() const;

  /** \deprecated use MeshsetsManager interface instead
   */
  DEPRECATED MoFEMErrorCode print_cubit_materials_set() const;

  /**@}*/

  /** \name Updating entities (DO NOT USE THIS USE BitRefManager interface)*/

  /**@{*/

  /** \brief Get child entities form meshset containing parent entities
   * \deprecated do not us this use BitRefManager interface
   *
   * Search for refined entities of given type whose parent are entities in the
   * parent meshset. It can be used for example to transfer information about
   * boundary conditions to refined mesh or split mesh by interface
   * elements. It is used by function refine_MESHSET, to update MESHSET finite
   *elements.
   *
   * \param parent meshset
   * \param child_bit refinement level
   * \param type of refined entity
   * \param child_type meshset where child entities are stored (if the child
   *meshset is set to be the parent meshset, the parent would be updated with
   *the refined entities) \param recursive if true parent meshset is searched
   *recursively
   *
   **/
  DEPRECATED MoFEMErrorCode update_meshset_by_entities_children(
      const EntityHandle parent, const BitRefLevel &child_bit,
      const EntityHandle child, EntityType child_type,
      const bool recursive = false, int verb = -1);

  /** \brief update fields meshesets by child entities
   * \deprecated do not us this use BitRefManager interface
   * \ingroup mofem_field
   */
  DEPRECATED MoFEMErrorCode update_field_meshset_by_entities_children(
      const BitRefLevel &child_bit, int verb = -1);

  /** \brief update field mesheset by child entities
   * \deprecated do not us this use BitRefManager interface
   * \ingroup mofem_field
   */
  DEPRECATED MoFEMErrorCode update_field_meshset_by_entities_children(
      const std::string name, const BitRefLevel &child_bit, int verb = -1);

  /** \brief update finite element mesheset by child entities
   * \deprecated do not us this use BitRefManager interface
   */
  DEPRECATED MoFEMErrorCode update_finite_element_meshset_by_entities_children(
      const std::string name, const BitRefLevel &child_bit,
      const EntityType fe_ent_type, int verb = -1);

  /**@}*/

  /**@}*/

  /** \name Shift BitRefLevl */

  /**@{*/

  /** \deprecated use BitRefManager to do this operation
   * \brief right shift bit ref level
   * \todo Should be outsourced to separate interface, i.e. BitLevelManage
   */
  DEPRECATED MoFEMErrorCode shift_right_bit_ref(const int shift, int verb = -1);

  /**@}*/


  /** \brief add TET elements from given refinement level to finite element
   * database given by name \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_bit_ref with mask explicitly
   * given
   *
   * \param BitRefLevel bit
   * \param finite element name
   * \param finite elenent type
   * \param verrbose level
   */
  DEPRECATED virtual MoFEMErrorCode
  add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,
                                                const std::string &name,
                                                EntityType type,
                                                int verb = -1) = 0;

  /** \brief add TET entities from given refinement level to finite element
   * database given by name \ingroup mofem_fe
   *
   * \deprecated use add_ents_to_finite_element_by_bit_ref with mask explicitly
   * given
   *
   * \param BitRefLevel bit
   * \param BitRefLevel mask
   * \param finite element name
   * \param finite element type
   * \param verrbose level
   */
  DEPRECATED virtual MoFEMErrorCode
  add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,
                                                const BitRefLevel &mask,
                                                const std::string &name,
                                                EntityType type,
                                                int verb = -1) = 0;

  /** \name Build problems (DEPRECATED SHOULD NOT USE THIS)*/

  /**@{*/

  /** \brief build problem data structures
  *
  * \note If square_matrix is set to true, that indicate that problem is
  structurally
  * symmetric, i.e. rows and columns have the same dofs and are indexed in the
  same
  * way.

  \deprecated Use ProblemsManager to build and partition problems

  *
  * @param  problem pointer
  * @param  square_matrix structurally symmetric problem
  * @param  verb          verbosity level
  * @return               error code
  *
  */
  DEPRECATED MoFEMErrorCode build_problem(const std::string &name,
                                          const bool square_matrix,
                                          int verb = -1);

  /** \brief build problem data structures
  *
  * \note If square_matrix is set to true, that indicate that problem is
  structurally
  * symmetric, i.e. rows and columns have the same dofs and are indexed in the
  same
  * way.

  \deprecated Use ProblemsManager to build and partition problems

  * @param  name          problem name
  * @param  square_matrix structurally symmetric problem
  * @param  verb          verbosity level
  * @return               error code
  *
  */
  DEPRECATED MoFEMErrorCode build_problem(Problem *problem_ptr,
                                          const bool square_matrix,
                                          int verb = -1);

  /** \brief build problem data structures

  \deprecated Use MoFEM::Interface::build_problem(const std::string &name,const
  bool square_matrix,int verb = -1) instead. This function not allows to Control
  if problem is structurally symmetric.

  */
  DEPRECATED virtual MoFEMErrorCode build_problems(int verb = -1) = 0;

  /** \brief build problem data structures, assuming that mesh is distributed
  (collective)

  \deprecated Use ProblemsManager to build and partition problems

  Mesh is distributed, that means that each processor keeps only own part of
  the mesh and shared entities.

  Collective - need to be run on all processors in communicator, i.e. each
  function has to call this function.

  */
  DEPRECATED MoFEMErrorCode build_problem_on_distributed_mesh(
      const std::string &name, const bool square_matrix = true, int verb = -1);

  /** \brief build problem data structures, assuming that mesh is distributed
  (collective)

  \deprecated Use ProblemsManager to build and partition problems

  Mesh is distributed, that means that each processor keeps only own part of
  the mesh and shared entities.

  Collective - need to be run on all processors in communicator, i.e. each
  function has to call this function.

  */
  DEPRECATED MoFEMErrorCode build_problem_on_distributed_mesh(
      Problem *problem_ptr, const bool square_matrix = true, int verb = -1);

  /** \brief build problem data structures, assuming that mesh is distributed
  (collective)

  \deprecated Use ProblemsManager to build and partition problems

  Mesh is distributed, that means that each processor keeps only own part of
  the mesh and shared entities.

  Collective - need to be run on all processors in communicator, i.e. each
  function has to call this function.

  */
  DEPRECATED virtual MoFEMErrorCode
  build_problem_on_distributed_mesh(int verb = -1) = 0;

  /**
   * \brief Set partition tag to each finite element in the problem
   *
   * This will use one of the mesh partitioning programs available from PETSc
   * See
   * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatPartitioningType.html>
   *
   * @param  ents        Entities to partition
   * @param  dim         Dimension of entities to partition
   * @param  adj_dim     Adjacency dimension
   * @param  n_parts     Number of partitions
   * @param  verb        Verbosity level
   * @return             Error code
   */
  DEPRECATED MoFEMErrorCode partition_mesh(const Range &ents, const int dim,
                                           const int adj_dim, const int n_parts,
                                           int verb = -1);

  /** \brief partition problem dofs

  \deprecated Use ProblemsManager to build and partition problems

  *
  * \param name problem name
  */
  DEPRECATED MoFEMErrorCode partition_simple_problem(const std::string &name,
                                                     int verb = -1);

  /** \brief partition problem dofs (collective)

  \deprecated Use ProblemsManager to build and partition problems

  *
  * \param name problem name
  */
  DEPRECATED MoFEMErrorCode partition_problem(const std::string &name,
                                              int verb = -1);

  /**
  * \brief build indexing and partition problem inheriting indexing and
  partitioning from two other problems
  *
  \deprecated Use ProblemsManager to build and partition problems
  * \param name problem name
  * \param problem_for_rows problem used to index rows
  * \param copy_rows just copy rows dofs
  * \param problem_for_cols problem used to index cols
  * \param copy_cols just copy cols dofs
  *
  * If copy_rows/copy_cols is set to false only partition is copied between
  problems.
  *
  */
  DEPRECATED MoFEMErrorCode partition_compose_problem(
      const std::string &name, const std::string &problem_for_rows,
      const bool copy_rows, const std::string &problem_for_cols,
      const bool copy_cols, int verb = -1);

  /**
  * \brief build sub problem

  \deprecated Use ProblemsManager to build and partition problems

  * @param  out_name problem
  * @param  fields_row  vector of fields composing problem
  * @param  fields_col  vector of fields composing problem
  * @param  main_problem main problem
  * @return              error code
  */
  DEPRECATED MoFEMErrorCode build_sub_problem(
      const std::string &out_name, const std::vector<std::string> &fields_row,
      const std::vector<std::string> &fields_col,
      const std::string &main_problem, const bool square_matrix = true,
      int verb = -1);

  /** \brief determine ghost nodes
   * \ingroup mofem_field
   * \deprecated Use ProblemsManager to build and partition problems
   *
   * \param name problem name
   */
  DEPRECATED MoFEMErrorCode partition_ghost_dofs(const std::string &name,
                                                 int verb = -1);

  /** \brief partition finite elements
   *
   * Function which partition finite elements based on dofs partitioning.<br>
   * In addition it sets information about local row and cols dofs at given
   * element on partition.
   *
   * \deprecated Use ProblemsManager to build and partition problems
   *
   * \param name problem name
   */
  DEPRECATED MoFEMErrorCode partition_finite_elements(
      const std::string &name, bool part_from_moab = false, int low_proc = -1,
      int hi_proc = -1, int verb = -1);

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
  DEPRECATED MoFEMErrorCode get_problem_elements_layout(
      const std::string &name, const std::string &fe_name, PetscLayout *layout,
      int verb = -1);

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
  DEPRECATED MoFEMErrorCode VecCreateSeq(const std::string &name, RowColData rc,
                                         Vec *V) const;

  /** \deprecated use VecManager
  * \brief create ghost vector for problem (collective)
  * \ingroup mofem_vectors

  collective - need to be run on all processors in communicator

  * \param name problem name
  * \param RowColData specify what data is taken from Row, Col or Data
  * \param Vec the vector where data is stored
  */
  DEPRECATED MoFEMErrorCode VecCreateGhost(const std::string &name,
                                           RowColData rc, Vec *V) const;

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
  DEPRECATED MoFEMErrorCode ISCreateFromProblemFieldToOtherProblemField(
      const std::string &x_problem, const std::string &x_field_name,
      RowColData x_rc, const std::string &y_problem,
      const std::string &y_field_name, RowColData y_rc, std::vector<int> &idx,
      std::vector<int> &idy, int verb = -1) const;

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
  DEPRECATED MoFEMErrorCode ISCreateFromProblemFieldToOtherProblemField(
      const std::string &x_problem, const std::string &x_field_name,
      RowColData x_rc, const std::string &y_problem,
      const std::string &y_field_name, RowColData y_rc, IS *ix, IS *iy,
      int verb = -1) const;

  /** \deprecated Use ISManager
  * \brief create IS for given order range (collective)

  * \param problem name
  * \param rc ROW or COL
  * \param min_order
  * \param max_order
  * \retval is out value

  */
  DEPRECATED MoFEMErrorCode ISCreateProblemOrder(const std::string &problem,
                                                 RowColData rc, int min_order,
                                                 int max_order, IS *is,
                                                 int verb = -1) const;

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
  DEPRECATED MoFEMErrorCode ISCreateProblemFieldAndRank(
      const std::string &problem, RowColData rc, const std::string &field,
      int min_coeff_idx, int max_coeff_idx, IS *is, int verb = -1) const;

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
  DEPRECATED MoFEMErrorCode VecScatterCreate(
      Vec xin, const std::string &x_problem, const std::string &x_field_name,
      RowColData x_rc, Vec yin, const std::string &y_problem,
      const std::string &y_field_name, RowColData y_rc, VecScatter *newctx,
      int verb = -1) const;

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
  DEPRECATED MoFEMErrorCode
  VecScatterCreate(Vec xin, const std::string &x_problem, RowColData x_rc,
                   Vec yin, const std::string &y_problem, RowColData y_rc,
                   VecScatter *newctx, int verb = -1) const;

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
   * \param scatter_mode see petsc manual for ScatterMode (The available modes
   * are: SCATTER_FORWARD or SCATTER_REVERSE)
   *
   * SCATTER_REVERSE set data to field entities from V vector.
   *
   * SCATTER_FORWARD set vector V from data field entities
   *
   */
  DEPRECATED MoFEMErrorCode
  set_local_ghost_vector(const Problem *problem_ptr, RowColData rc, Vec V,
                         InsertMode mode, ScatterMode scatter_mode) const;

  /** \deprecated use VecManager
   * \brief set values of vector from/to meshdatabase
   * \ingroup mofem_vectors
   *
   * \param name of the problem
   * \param RowColData for row or column:e (i.e. Row,Col)
   * \param V vector
   * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
   * \param scatter_mode see petsc manual for ScatterMode (The available modes
   * are: SCATTER_FORWARD or SCATTER_REVERSE)
   *
   * SCATTER_REVERSE set data to field entities from V vector.
   *
   * SCATTER_FORWARD set vector V from data field entities
   *
   */
  DEPRECATED MoFEMErrorCode
  set_local_ghost_vector(const std::string &name, RowColData rc, Vec V,
                         InsertMode mode, ScatterMode scatter_mode) const;

  /** \deprecated use VecManager
  * \brief set values of vector from/to mesh database (collective)
  * \ingroup mofem_vectors

  collective - need tu be run on all processors in communicator

  * \param pointer to porblem struture
  * \param RowColData for row or column (i.e. Row,Col)
  * \param V vector
  * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
  * \param scatter_mode see petsc manual for ScatterMode (The available modes
  are: SCATTER_FORWARD or SCATTER_REVERSE)
  *
  * SCATTER_REVERSE set data to field entities form V vector.
  *
  */
  DEPRECATED MoFEMErrorCode
  set_global_ghost_vector(const Problem *problem_ptr, RowColData rc, Vec V,
                          InsertMode mode, ScatterMode scatter_mode) const;

  /** \deprecated use VecManager
  * \brief set values of vector from/to mesh database (collective)
  * \ingroup mofem_vectors

  collective - need tu be run on all processors in communicator

  * \param name of the problem
  * \param RowColData for row or column (i.e. Row,Col)
  * \param V vector
  * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
  * \param scatter_mode see petsc manual for ScatterMode (The available modes
  are: SCATTER_FORWARD or SCATTER_REVERSE)
  *
  * SCATTER_REVERSE set data to field entities form V vector.
  *
  */
  DEPRECATED MoFEMErrorCode
  set_global_ghost_vector(const std::string &name, RowColData rc, Vec V,
                          InsertMode mode, ScatterMode scatter_mode) const;

  /** \deprecated use VecManager
   * \brief Copy vector to field which is not part of the problem
   * \ingroup mofem_vectors
   *
   * \param pointer to poroblem multi_index
   * \param field_name field name used for indexing petsc vectors used in the
   * problem \param cpy_field field name where data from vector are stored
   * \param RowColData for row or column
   * \param V vector
   * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
   * \param scatter_mode see petsc manual for ScatterMode (The available modes
   * are: SCATTER_FORWARD or SCATTER_REVERSE)
   *
   * SCATTER_REVERSE set data to field entities form V vector.
   *
   */
  DEPRECATED MoFEMErrorCode set_other_local_ghost_vector(
      const Problem *problem_ptr, const std::string &fiel_name,
      const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
      ScatterMode scatter_mode, int verb = -1);

  /** \deprecated use VecManager
   * \brief Copy vector to field which is not part of the problem
   * \ingroup mofem_vectors
   *
   * \param name problem name
   * \param field_name field name used for indexing petsc vectors used in the
   * problem \param cpy_field field name where data from vector are stored
   * \param RowColData for row or column
   * \param V vector
   * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
   * \param scatter_mode see petsc manual for ScatterMode (The available modes
   * are: SCATTER_FORWARD or SCATTER_REVERSE)
   *
   * SCATTER_REVERSE set data to field entities form V vector.
   *
   */
  DEPRECATED MoFEMErrorCode set_other_local_ghost_vector(
      const std::string &name, const std::string &field_name,
      const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
      ScatterMode scatter_mode, int verb = -1);

  /** \deprecated use VecManager
  * \brief Copy vector to field which is not part of the problem (collective)
  * \ingroup mofem_vectors

  collective - need tu be run on all processors in communicator

  * \param problem_ptr pointer to problem
  * \param field_name field name used for indexing petsc vectors used in the
  problem
  * \param cpy_field field name where data from vector are stored
  * \param RowColData for row or column
  * \param V vector
  * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
  * \param scatter_mode see petsc manual for ScatterMode (The available modes
  are: SCATTER_FORWARD or SCATTER_REVERSE)
  *
  * SCATTER_REVERSE set data to field entities form V vector.
  *
  */
  DEPRECATED MoFEMErrorCode set_other_global_ghost_vector(
      const Problem *problem_ptr, const std::string &field_name,
      const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
      ScatterMode scatter_mode, int verb = -1);

  /** \deprecated use VecManager
  * \brief Copy vector to field which is not part of the problem (collective)
  * \ingroup mofem_vectors

  collective - need tu be run on all processors in communicator

  * \param name problem name
  * \param field_name field name used for indexing petsc vectors used in the
  problem
  * \param cpy_field field name where data from vector are stored
  * \param RowColData for row or column
  * \param V vector
  * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
  * \param scatter_mode see petsc manual for ScatterMode (The available modes
  are: SCATTER_FORWARD or SCATTER_REVERSE)
  *
  * SCATTER_REVERSE set data to field entities form V vector.
  *
  */
  DEPRECATED MoFEMErrorCode set_other_global_ghost_vector(
      const std::string &name, const std::string &field_name,
      const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
      ScatterMode scatter_mode, int verb = -1);

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
   * \param error_if_missing throw error if entity/dof exist in field_x but not
   * on field_y \param create_if_missing creat dof in field_y from fiedl_x if it
   * is not database
   *
   */
  DEPRECATED MoFEMErrorCode field_axpy(const double alpha,
                                       const std::string &fiel_name_x,
                                       const std::string &field_name_y,
                                       bool error_if_missing = false,
                                       bool creat_if_missing = false);

  /** \deprecated use FieldBlas
   * \brief scale field
   * \ingroup mofem_field_algebra
   * \todo should be moved to independent interface, i.e. FieldAlgebra
   *
   * \param alpha is a scaling factor
   * \field_name  is a field name
   *
   */
  DEPRECATED MoFEMErrorCode field_scale(const double alpha,
                                        const std::string &field_name);

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
  DEPRECATED MoFEMErrorCode set_field(const double val, const EntityType type,
                                      const std::string &field_name);

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
  DEPRECATED MoFEMErrorCode set_field(const double val, const EntityType type,
                                      const Range &ents,
                                      const std::string &field_name);

  /**@}*/

  /** \name Get adjacencies */

  /**@{*/

  /** \brief Get the adjacencies associated with a entity to entities of a
   * specified dimension. \ingroup mofem_ref_ents \todo Should be outsourced to
   * separate interface, i.e. BitLevelManager
   *
   * bit ref level of adjacent entities is equal to bit ref level of adjacent
   * entities
   */
  DEPRECATED MoFEMErrorCode
  get_adjacencies_equality(const EntityHandle from_entiti,
                           const int to_dimension, Range &adj_entities) const;

  /** \brief Get the adjacencies associated with a entity to entities of a
   * specified dimension. \ingroup mofem_ref_ents
   *
   * bit ref level of adjacent entities is any of bit ref level of adjacent
   * entities
   */
  DEPRECATED MoFEMErrorCode get_adjacencies_any(const EntityHandle from_entiti,
                                                const int to_dimension,
                                                Range &adj_entities) const;

  /** \brief Get the adjacencies associated with a entity to entities of a
   * specified dimension. \ingroup mofem_ref_ents \todo Should be outsourced to
   * separate interface, i.e. BitLevelManage
   *
   * bit ref level of adjacent entities is equal to bit ref level of adjacent
   * entities
   */
  DEPRECATED MoFEMErrorCode get_adjacencies(
      const Problem *problem_ptr, const EntityHandle *from_entities,
      const int num_netities, const int to_dimension, Range &adj_entities,
      const int operation_type = moab::Interface::INTERSECT,
      const int verb = 0) const;

  /** \brief Get the adjacencies associated with a entity to entities of a
   * specified dimension. \ingroup mofem_ref_ents \todo Should be outsourced to
   * separate interface, i.e. BitLevelManage
   *
   * bit ref level of adjacent entities is equal to bit ref level of adjacent
   * entities
   */
  DEPRECATED MoFEMErrorCode get_adjacencies(
      const BitRefLevel &bit, const EntityHandle *from_entities,
      const int num_netities, const int to_dimension, Range &adj_entities,
      const int operation_type = moab::Interface::INTERSECT,
      const int verb = 0) const;

  /**@}*/

  /**@}*/

  /** \name Clear and remove */

  /**@{*/

  /** \deprecated Clear dofs by bit level
   */
  DEPRECATED MoFEMErrorCode clear_dofs_fields(const BitRefLevel &bit,
                                              const BitRefLevel &mask,
                                              int verb = -1);

  /** \deprecated use clear_ents_fields_by_bit_ref
   */
  DEPRECATED MoFEMErrorCode clear_ents_fields(const BitRefLevel &bit,
                                              const BitRefLevel &mask,
                                              int verb = -1);

  /** \deprecated use clear_finite_elements_by_bit_ref
   */
  DEPRECATED MoFEMErrorCode clear_finite_elements(const BitRefLevel &bit,
                                                  const BitRefLevel &mask,
                                                  int verb = -1);

  /**@}*/

  /**@}*/

  /** \name Clear and remove */

  /**@{*/

  /**
   * \brief set field entities on vertices
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param nodes contains set vertices
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_VERTICEs(
      const Range &nodes, const std::string &name, int verb = -1);

  /**
   * \brief set field entities on vertices
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param meshset contains set vertices
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_VERTICEs(
      const EntityHandle meshset, const std::string &name, int verb = -1);

  /**
   * \brief set field entities form adjacencies of edges
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param range contains set edges
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_EDGEs(const Range &edges,
                                                       const std::string &name,
                                                       int verb = -1);

  /**
   * \brief set field entities form adjacencies of edges
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param meshset contains set edges
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_EDGEs(
      const EntityHandle meshset, const std::string &name, int verb = -1);

  /**
   * \brief set field entities form adjacencies of triangles
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param meshset contains set triangles
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_TRIs(
      const EntityHandle meshset, const std::string &name, int verb = -1);

  /**
   * \brief set field entities form adjacencies of triangles
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param range triangles
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_TRIs(const Range &tris,
                                                      const std::string &name,
                                                      int verb = -1);

  /**
   * \brief set field entities from adjacencies of tetrahedron
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param meshset contains set tetrahedron
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_TETs(
      const EntityHandle meshset, const std::string &name, int verb = -1);

  /**
   * \brief set field entities from adjacencies of tetrahedron
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param range contains set tetrahedron
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_TETs(const Range &tets,
                                                      const std::string &name,
                                                      int verb = -1);

  /**
   * \deprecated use add_ents_to_field_by_type
   * \brief set field entities from adjacencies of quads
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param quads range contains set quads
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_QUADs(const Range &quads,
                                                       const std::string &name,
                                                       int verb = -1);

  /**
   * \deprecated use add_ents_to_field_by_type
   * \brief set field entities from adjacencies of quads
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param meshset contains set quads
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_QUADs(EntityHandle meshset,
                                                       const std::string &name,
                                                       int verb = -1);

  /**
   * \deprecated use add_ents_to_field_by_type
   * \brief set field entities from adjacencies of prisms
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param prisms range contains set prisms
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_PRISMs(const Range &prisms,
                                                        const std::string &name,
                                                        int verb = -1);

  /**
   * \deprecated use add_ents_to_field_by_type
   * \brief set field entities from adjacencies of prisms
   * \ingroup mofem_field
   *
   * The lower dimension entities are added depending on the space type
   * \param meshset contains set prisms
   * \param name of the field
   */
  DEPRECATED MoFEMErrorCode add_ents_to_field_by_PRISMs(EntityHandle meshset,
                                                        const std::string &name,
                                                        int verb = -1);

  /**@}*/

  /** \name Create matrices (will be moved to independent interface) */

  /**@{*/

  /**
   * @deprecated Use MatrixInterface
   *
   * \code
   * CHKERR m_field.getInterface<MatrixManager>()
   *  ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(name, Aij, verb);
   * \endcode
   */
  DEPRECATED MoFEMErrorCode MatCreateMPIAIJWithArrays(
      const std::string &name, Mat *Aij, int verb = DEFAULT_VERBOSITY);

  /**
   * @deprecated Use MatrixInterface
   *
   * \code
   * CHKERR m_field.getInterface<MatrixManager>()
   *  ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(name, Aij, verb);
   * \endcode
   */
  DEPRECATED MoFEMErrorCode MatCreateMPIAdj_with_Idx_mi_tag(
      const std::string &name, Mat *Adj, int verb = DEFAULT_VERBOSITY);

  /**
   * @deprecated Use MatrixInterface
   *
   * \code
   * CHKERR m_field.getInterface<MatrixManager>()
   *  ->createSeqAIJWithArrays<PetscLocalIdx_mi_tag>(name, Adj, verb);
   * \endcode
   *
   */
  DEPRECATED MoFEMErrorCode MatCreateSeqAIJWithArrays(
      const std::string &name, Mat *Aij, PetscInt **i, PetscInt **j,
      PetscScalar **v, int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /** \name Making loops on elements and entities */

  /**@{*/

  using CoreInterface::loop_finite_elements;

  /**
   * @deprecated Use version from core interface
   */
  DEPRECATED MoFEMErrorCode loop_finite_elements(const Problem *problem_ptr,
                                                 const std::string &fe_name,
                                                 FEMethod &method,
                                                 int lower_rank, int upper_rank,
                                                 MoFEMTypes bh,
                                                 int verb = DEFAULT_VERBOSITY);

  /**
   * @deprecated Use version from core interface
   */
  DEPRECATED MoFEMErrorCode loop_finite_elements(
      const std::string &problem_name, const std::string &fe_name,
      FEMethod &method, int lower_rank, int upper_rank, MoFEMTypes bh,
      int verb = DEFAULT_VERBOSITY);

  /**
   * @deprecated Use version from core interface
   */
  DEPRECATED MoFEMErrorCode loop_finite_elements(
      const std::string &problem_name, const std::string &fe_name,
      FEMethod &method, MoFEMTypes bh, int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /** \name Comm */

  /**@{*/

  /**
   * @brief make entities from proc 0 shared on all proc
   * @deprecated Use CommInterface
   *
   * \note collective - need tu be run on all processors in communicator
   *
   * @param entities
   * @param num_entities
   * @param my_proc default proc id to share from
   * @param verb
   * @return MoFEMErrorCode
   */
  DEPRECATED MoFEMErrorCode make_entities_multishared(
      const EntityHandle *entities, const int num_entities,
      const int my_proc = 0, int verb = DEFAULT_VERBOSITY);
  /**
   * @brief make entities from proc 0 shared on all proc
   * @deprecated Use CommInterface
   *
   * \note collective - need tu be run on all processors in communicator
   *
   * @param entities
   * @param my_proc default proc id to share from
   * @param verb
   * @return MoFEMErrorCode
   */
  DEPRECATED MoFEMErrorCode make_entities_multishared(
      Range &entities, const int my_proc = 0, int verb = DEFAULT_VERBOSITY);

  /**
   * @brief make field entities multi shared
   * @deprecated Use CommInterface
   *
   * \note collective - need tu be run on all processors in communicator
   *
   * @param field_name
   * @param owner_proc
   * @param verb
   * @return MoFEMErrorCode
   */
  DEPRECATED MoFEMErrorCode make_field_entities_multishared(
      const std::string field_name, const int owner_proc = 0,
      int verb = DEFAULT_VERBOSITY);

  /**
   * @brief Exchange field data
   * @deprecated Use CommInterface
   *
   * Exchange field for all shared and ghosted entities. This function should be
   * called collectively over the communicator for this ParallelComm. If the
   * entities vector is empty, all shared entities participate in the exchange.
   * If a proc has no owned entities this function must still be called since it
   * is collective.
   *
   * \note collective - need tu be run on all processors in communicator
   *
   * \todo It is not working if field has entities diffrent than vertices.
   *
   * @param verb
   * @param field_name @return MoFEMErrorCode
   */
  DEPRECATED MoFEMErrorCode exchange_field_data(const std::string field_name,
                                                int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /** \name Synchronize */

  /**@{*/

  /** synchronize entity range on processors (collective)

  collective - need tu be run on all processors in communicator

  @deprecated Use Comm Interface
  */
  DEPRECATED MoFEMErrorCode synchronise_entities(Range &ent,
                                                 int verb = DEFAULT_VERBOSITY);

  /** synchronize entity range on processors (collective)
  * \ingroup mofem_field

  collective - need tu be run on all processors in communicator

  \param name field
  \param verbose level

  \deprecated Use CommInterface
  */
  DEPRECATED MoFEMErrorCode synchronise_field_entities(
      const std::string &name, int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /** \name Delete and remove */

  /**@{*/

  /**
   * \deprecated use  remove_parents_by_bit_ref
   */
  DEPRECATED MoFEMErrorCode
  remove_parents_by_by_bit_ref(const BitRefLevel bit, const BitRefLevel mask,
                               int verb = DEFAULT_VERBOSITY);

  /**@}*/

};

} // namespace MoFEM

#endif // __INTERFACE_DEPRECATED_HPP__
