/** \file Core.hpp
 * \brief Core Interface class for user interface
 *
 * Low level data structures not used directly by user
 *
 * FIXME It is a mess with names of core cpp files need better organization
 *
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

#ifndef __MOABFIELD_CORE_HPP__
#define __MOABFIELD_CORE_HPP__

namespace MoFEM {

struct MeshsetsManager;

/** \brief Core (interface) class
 *  \ingroup mofem

 This is the implementation of abstract MoFEM::Interface class. Similarly to the
 convention used in MoAB, we use small letters to name function of purely
 abstract classes. This is an exception used only here. For more details about
 naming functions see \ref coding_practice

 This class is not used directly by the user. For internal use only.   It is
 database with basic functions to access data. Abstraction of this is MoFEM
 Interface structure.

 It is deign to hide come complexities for users and allow low development
 without interfering with users modules programmer work.

 */
struct Core: public Interface {

  moab::Interface& moab;

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);
  PetscErrorCode query_interface_type(const std::type_info& iface_type,void*& ptr) const;

  /**
   * Contruct core database
   */
  Core(
    moab::Interface& moab,              ///< MoAB interface
    MPI_Comm comm = PETSC_COMM_WORLD,   ///< MPI communicator
    int verbose = 1                     ///< Verbosity level
  );
  ~Core();

  inline Tag get_th_RefParentHandle() const { return th_RefParentHandle; }
  inline Tag get_th_RefBitLevel() const { return th_RefBitLevel; }
  inline Tag get_th_RefBitEdge() const { return th_RefBitEdge; }
  inline Tag get_th_RefType() const { return th_RefType; }

  /**
   * Is used to check consistency. I n future properly this will be removed and
   * replaced by other solution. It is only for internal use.
   */
  enum SemaphoresBuildMofem {
    BUILD_FIELD = 1<<0,
    BUILD_FE = 1<<1,
    BUILD_ADJ = 1<<2,
    BUILD_PROBLEM = 1<<3,
    PARTITION_PROBLEM = 1<<4,
    PARTITION_FE = 1<<5,
    PARTITION_GHOST_DOFS = 1<<6,
    PARTITION_MESH = 1<<7
  };

  /**
   * \brief Get flags/semaphores for different stages
   */
  inline int& getBuildMoFEM() const { return *buildMoFEM; }

  //add prims element FIXME This is wrong solution
  PetscErrorCode addPrismToDatabase(const EntityHandle prism,int verb = -1);

  protected:

  mutable boost::ptr_map<unsigned long,UnknownInterface> iFaces;
  mutable int *buildMoFEM; ///< keeps flags/semaphores for different stages

  //Database
  MoABErrorCode rval;
  PetscErrorCode ierr;

  //Data and low level methods
  Tag th_Part;  ///< Tag for partition number
  Tag th_RefParentHandle,th_RefBitLevel,th_RefBitLevel_Mask,th_RefBitEdge,th_RefFEMeshset;
  Tag th_RefType;
  Tag th_FieldId,th_FieldName,th_FieldName_DataNamePrefix,th_FieldSpace,th_FieldBase;
  Tag th_FEId,th_FEName;
  Tag th_FEIdCol,th_FEIdRow,th_FEIdData;
  Tag th_ProblemId,th_ProblemName,th_ProblemFEId;
  Tag th_ProblemNbDofsRow,th_ProblemNbDofsCol;
  Tag th_ProblemLocalNbDofRow,th_ProblemGhostNbDofRow;
  Tag th_ProblemLocalNbDofCol,th_ProblemGhostNbDofCol;
  Tag th_ProblemShift,th_FieldShift,th_FEShift;
  Tag th_ElemType;                    ///< Needed for VTK files

  boost::shared_ptr<BasicEntityData> basicEntityDataPtr;

  /**
   * \brief Get pointer to basic entity data.
   *
   * This structure keeps data like tags handlers and other data used to construct
   * mofem entities, dofs and finite elements.
   *
   */
  boost::shared_ptr<BasicEntityData> get_basic_entity_data_ptr() {
    return basicEntityDataPtr;
  }

  int *fShift,*feShift,*pShift;
  int verbose;

  // Managing and storing basic entities
  RefEntity_multiIndex refinedEntities;		       ///< refined entities
  RefElement_multiIndex refinedFiniteElements;	 ///< refined elements

  // Managing and storung DOFs
  Field_multiIndex fIelds;			           ///< field
  FieldEntity_multiIndex entsFields;			 ///< entities on field
  DofEntity_multiIndex dofsField;		       ///< dofs on fiels

  // Managing and storing finite elements
  FiniteElement_multiIndex finiteElements;		        ///< finite elements
  EntFiniteElement_multiIndex entsFiniteElements;			///< finite element entities

  // Managing and storing adjacencies
  FieldEntityEntFiniteElementAdjacencyMap_multiIndex entFEAdjacencies;	///< adjacencies of elements to dofs

  //pRoblems
  Problem_multiIndex pRoblems;					 ///< problems

  //safety nets
  Tag th_MoFEMBuild;

  //core methods
  PetscErrorCode getTags(int verb = -1);
  PetscErrorCode clearMap();
  BitFieldId getFieldShift();
  BitFEId getFEShift();
  BitProblemId getProblemShift();
  PetscErrorCode initialiseDatabseInformationFromMesh(int verb = -1);

  //moab interface
  inline moab::Interface& get_moab() { return moab; }
  inline const moab::Interface& get_moab() const { return moab; }

  //FiedlInterface

  //check consistency
  PetscErrorCode check_number_of_ents_in_ents_field(const std::string& name) const;
  PetscErrorCode check_number_of_ents_in_ents_field() const;
  PetscErrorCode check_number_of_ents_in_ents_finite_element(const std::string& name) const;
  PetscErrorCode check_number_of_ents_in_ents_finite_element() const;

  PetscErrorCode clear_database(int verb  = -1);
  PetscErrorCode rebuild_database(int verb = -1);

  //cubit meshsets

  MeshsetsManager* get_meshsets_manager_ptr();
  const MeshsetsManager* get_meshsets_manager_ptr() const;
  inline MeshsetsManager& get_meshsets_manager() {
    return *get_meshsets_manager_ptr();
  }
  inline const MeshsetsManager& get_meshsets_manager() const {
    return *get_meshsets_manager_ptr();
  }

  // Should not be used, access data meshesets by MeshsetsManager interface
  DEPRECATED bool check_msId_meshset(const int ms_id,const CubitBCType cubit_bc_type);
  DEPRECATED PetscErrorCode add_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id,const std::string name = "");
  DEPRECATED PetscErrorCode set_cubit_msId_attribites(
    const CubitBCType cubit_bc_type,const int ms_id,const std::vector<double> &attributes,const std::string name = ""
  );
  DEPRECATED PetscErrorCode set_cubit_msId_attribites_data_structure(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericAttributeData &data,const std::string name = ""
  );
  DEPRECATED PetscErrorCode set_cubit_msId_bc_data_structure(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericCubitBcData &data
  );
  DEPRECATED PetscErrorCode delete_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id);
  DEPRECATED PetscErrorCode get_cubit_msId(
    const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr
  );
  DEPRECATED PetscErrorCode get_cubit_msId_entities_by_dimension(
    const int ms_id,const CubitBCType cubit_bc_type, const int dimension,Range &entities,const bool recursive = false
  );
  DEPRECATED PetscErrorCode get_cubit_msId_entities_by_dimension(
    const int ms_id,const CubitBCType cubit_bc_type, Range &entities,const bool recursive = false
  );
  DEPRECATED PetscErrorCode get_cubit_msId_entities_by_dimension(
    const int ms_id,const unsigned int cubit_bc_type, const int dimension,Range &entities,const bool recursive = false
  );
  DEPRECATED PetscErrorCode get_cubit_msId_entities_by_dimension(
    const int ms_id,const unsigned int cubit_bc_type, Range &entities,const bool recursive = false
  );
  DEPRECATED PetscErrorCode get_cubit_msId_meshset(
    const int ms_id,const unsigned int cubit_bc_type,EntityHandle &meshset
  );
  DEPRECATED PetscErrorCode get_cubit_meshsets(
    const unsigned int cubit_bc_type,Range &meshsets
  );
  DEPRECATED PetscErrorCode print_cubit_displacement_set() const;
  DEPRECATED PetscErrorCode print_cubit_pressure_set() const;
  DEPRECATED PetscErrorCode print_cubit_force_set() const;
  DEPRECATED PetscErrorCode print_cubit_temperature() const;
  DEPRECATED PetscErrorCode print_cubit_heat_flux_set() const;
  DEPRECATED PetscErrorCode print_cubit_materials_set() const;

  //refine
  PetscErrorCode seed_finite_elements(const Range &entities,int verb = -1);
  PetscErrorCode seed_finite_elements(const EntityHandle meshset,int verb = -1);
  PetscErrorCode seed_ref_level(const Range &ents,const BitRefLevel &bit,const bool only_tets = true,int verb = -1);
  PetscErrorCode seed_ref_level_2D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb = -1);
  PetscErrorCode get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb = -1);
  PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset);
  PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents);
  // PetscErrorCode add_ref_level_to_entities(const BitRefLevel &bit,Range &ents);
  // PetscErrorCode set_ref_level_to_entities(const BitRefLevel &bit,Range &ents);
  PetscErrorCode update_meshset_by_entities_children(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,
    const bool recursive = false, int verb = -1
  );
  PetscErrorCode update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb = -1);
  PetscErrorCode update_field_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,int verb = -1);
  PetscErrorCode update_finite_element_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb = -1);

  //remove entities
  PetscErrorCode delete_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const bool remove_parent = false,int verb = -1);
  PetscErrorCode remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode delete_finite_elements_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode shift_left_bit_ref(const int shif,int verb = -1);
  PetscErrorCode shift_right_bit_ref(const int shift,int verb = -1);

  //synchronize entities
  PetscErrorCode synchronise_entities(Range &ent,int verb = -1);
  PetscErrorCode synchronise_field_entities(const BitFieldId id,int verb = -1);
  PetscErrorCode synchronise_field_entities(const std::string& name,int verb = -1);

  /**
   * \brief Add filed
   * @param  name           Field name
   * @param  space          Space L2,H1,Hdiv,Hcurl
   * @param  base           Approximation base AINSWORTH_LEGENDRE_BASE, AINSWORTH_BERNSTEIN_BEZIER_BASE ...
   * @param  nb_cooficients Number of field coefficients
   * @param  tag_type       Tag type, MB_TAG_DENSE or MB_TAG_SPARSE (default)
   * @param  bh             Control behavior, if MF_EXCL throws error if exist
   * @param  verb           Verbosity level
   * @return                Return error code

   TODO: \todo MB_TAG_DENSE will not work properly in general case. It is need to separate field
   tags for each entity separately. That will allow for HO orders but homogenous approx. order
   on each entity. Need some discussion what is optimal solution. MB_TAG_SPARSE gives flexibility,
   but it not memory efficient. MB_TAG_DENSE uses memory more efficient and in principle allow
   for better efficiency if properly utilized.

   FIXME: \bug Need to resolve problem of dense tags at this stage of development will make only problems

   */
  PetscErrorCode add_field(
    const std::string& name,
    const FieldSpace space,
    const FieldApproximationBase base,
    const FieldCoefficientsNumber nb_cooficients,
    const TagType tag_type = MB_TAG_SPARSE,
    const enum MoFEMTypes bh = MF_EXCL,
    int verb = -1
  );

  PetscErrorCode add_ents_to_field_by_VERTICEs(const Range &nodes,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_VERTICEs(const Range &nodes,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_EDGEs(const Range &edges,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_EDGEs(const Range &edges,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_EDGEs(const EntityHandle meshset,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_EDGEs(const EntityHandle meshset,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TRIs(const EntityHandle meshset,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TRIs(const EntityHandle meshset,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TRIs(const Range &tris,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TRIs(const Range &tris,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TETs(const Range &tets,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TETs(const Range &tets,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_QUADs(const Range &prisms,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_QUADs(const Range &prisms,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_QUADs(EntityHandle meshset,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_PRISMs(const Range &prisms,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_PRISMs(const Range &prisms,const std::string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_PRISMs(EntityHandle meshset,const std::string& name,int verb = -1);
  PetscErrorCode remove_ents_from_field_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode remove_ents_from_field(const std::string& name,const EntityHandle meshset,const EntityType type,int verb = -1);
  PetscErrorCode remove_ents_from_field(const std::string& name,const Range &ents,int verb = -1);

  //set apprix oorder
  PetscErrorCode set_field_order(const Range &ents,const BitFieldId id,const ApproximationOrder order,int verb = -1);
  PetscErrorCode set_field_order(
    const EntityHandle meshset,
    const EntityType type,
    const BitFieldId id,
    const ApproximationOrder order,
    int verb = -1
  );
  PetscErrorCode set_field_order(const Range &ents,const std::string& name,const ApproximationOrder order,int verb = -1);
  PetscErrorCode set_field_order(
    const EntityHandle meshset,
    const EntityType type,
    const std::string& name,
    const ApproximationOrder order,
    int verb = -1
  );
  PetscErrorCode set_field_order_by_entity_type_and_bit_ref(
    const BitRefLevel &bit,
    const BitRefLevel &mask,
    const EntityType type,
    const BitFieldId id,
    const ApproximationOrder order,
    int verb = -1
  );
  PetscErrorCode set_field_order_by_entity_type_and_bit_ref(
    const BitRefLevel &bit,
    const BitRefLevel &mask,
    const EntityType type,
    const std::string& name,
    const ApproximationOrder order,
    int verb = -1
  );

  //build fiels
  PetscErrorCode buildFieldForNoField(const BitFieldId id,std::map<EntityType,int> &dof_counter,int verb = -1);
  PetscErrorCode buildFieldForL2H1HcurlHdiv(
    const BitFieldId id,std::map<EntityType,int> &dof_counter,
    std::map<EntityType,
    int> &inactive_dof_counter,
    int verb = -1
  );
  PetscErrorCode build_fields(int verb = -1);
  PetscErrorCode clear_inactive_dofs(int verb = -1);
  PetscErrorCode clear_dofs_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_ents_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_dofs_fields(const std::string &name,const Range ents,int verb = -1);
  PetscErrorCode clear_ents_fields(const std::string &name,const Range enst,int verb = -1);

  //other auxiliary functions for fields
  PetscErrorCode list_dofs_by_field_name(const std::string &name) const;
  PetscErrorCode list_fields() const;

  BitFieldId get_BitFieldId(const std::string& name) const;
  std::string get_BitFieldId_name(const BitFieldId id) const;
  EntityHandle get_field_meshset(const BitFieldId id) const;
  EntityHandle get_field_meshset(const std::string& name) const;
  PetscErrorCode get_field_entities_by_dimension(const std::string name,int dim,Range &ents) const;
  PetscErrorCode get_field_entities_by_type(const std::string name,EntityType type,Range &ents) const;
  PetscErrorCode get_field_entities_by_handle(const std::string name,Range &ents) const;



  bool check_field(const std::string& name) const;
  const Field* get_field_structure(const std::string& name);

  //FiniteElement
  bool check_finite_element(const std::string& name) const;
  PetscErrorCode add_finite_element(const std::string &fe_name,enum MoFEMTypes bh = MF_EXCL);
  PetscErrorCode modify_finite_element_adjacency_table(const std::string &fe_name,const EntityType type,ElementAdjacencyFunct function);
  PetscErrorCode modify_finite_element_add_field_data(const std::string &fe_name,const std::string &name_filed);
  PetscErrorCode modify_finite_element_add_field_row(const std::string &fe_name,const std::string &name_row);
  PetscErrorCode modify_finite_element_add_field_col(const std::string &fe_name,const std::string &name_col);
  PetscErrorCode modify_finite_element_off_field_data(const std::string &fe_name,const std::string &name_filed);
  PetscErrorCode modify_finite_element_off_field_row(const std::string &fe_name,const std::string &name_row);
  PetscErrorCode modify_finite_element_off_field_col(const std::string &fe_name,const std::string &name_col);
  PetscErrorCode add_ents_to_finite_element_by_type(
    const EntityHandle meshset,const EntityType type,const std::string &name,const bool recursive = true
  );
  PetscErrorCode add_ents_to_finite_element_by_dimension(
    const EntityHandle meshset,const int dim,const std::string &name,const bool recursive = true
  );
  PetscErrorCode add_ents_to_finite_element_by_type(const Range& ents,const EntityType type,const std::string &name);
  PetscErrorCode add_ents_to_finite_element_by_dim(const Range& ents,const int dim,const std::string &name);
  PetscErrorCode add_ents_to_finite_element_by_bit_ref(
    const BitRefLevel &bit,const BitRefLevel &mask,const std::string &name,EntityType type,int verb = -1
  );
  PetscErrorCode add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const std::string& name,const bool recursive = false);

  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_VERTICEs(const Range& vert,const std::string &name);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_EDGEs(const Range& vert,const std::string &name);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_EDGEs(const EntityHandle meshset,const std::string &name,const bool recursive = false);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_TRIs(const Range& tris,const std::string &name);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_TRIs(const EntityHandle meshset,const std::string &name,const bool recursive = false);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_TETs(const Range& tets,const std::string &name);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const std::string &name,const bool recursive = false);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_PRISMs(const Range& prims,const BitFEId id);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_PRISMs(const Range& prims,const std::string &name);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const std::string &name,const bool recursive = false);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const std::string &name,EntityType type,int verb = -1);
  DEPRECATED PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(
    const BitRefLevel &bit,const BitRefLevel &mask,const std::string &name,EntityType type,int verb = -1
  );
  PetscErrorCode remove_ents_from_finite_element_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode remove_ents_from_finite_element(const std::string &name,const EntityHandle meshset,const EntityType type,int verb = -1);
  PetscErrorCode remove_ents_from_finite_element(const std::string &name,const Range &ents,int verb = -1);
  PetscErrorCode delete_finite_element(const std::string name,int verb = -1);

  //other auxiliary functions for finite element
  BitFEId getBitFEId(const std::string& name) const;
  std::string getBitFEId_name(const BitFEId id) const;
  EntityHandle get_finite_element_meshset(const BitFEId id) const;
  EntityHandle get_finite_element_meshset(const std::string& name) const;
  PetscErrorCode get_finite_element_entities_by_dimension(
    const std::string name,int dim,Range &ents
  ) const;
  PetscErrorCode get_finite_element_entities_by_type(
    const std::string name,EntityType type,Range &ents
  ) const;
  PetscErrorCode get_finite_element_entities_by_handle(
    const std::string name,Range &ents
  ) const;
  PetscErrorCode list_finite_elements() const;

  //problem
  PetscErrorCode add_problem(const BitProblemId id,const std::string& name);
  PetscErrorCode add_problem(const std::string& name,enum MoFEMTypes bh = MF_EXCL,int verb = -1);
  bool check_problem(const std::string name);
  PetscErrorCode delete_problem(const std::string name);
  PetscErrorCode modify_problem_add_finite_element(const std::string &name_problem,const std::string &MoFEMFiniteElement_name);
  PetscErrorCode modify_problem_unset_finite_element(const std::string &name_problem,const std::string &MoFEMFiniteElement_name);
  PetscErrorCode modify_problem_ref_level_add_bit(const std::string &name_problem,const BitRefLevel &bit);
  PetscErrorCode modify_problem_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit);
  PetscErrorCode modify_problem_mask_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit);
  BitProblemId get_BitProblemId(const std::string& name) const;
  PetscErrorCode list_problem() const;
  PetscErrorCode clear_problem(const std::string &name,int verb = -1);
  PetscErrorCode clear_problems(int verb = -1);

  ///add entity EntFe to finite element data databse and resolve dofs on that entity
  //loop over all finite elements, resolve its meshsets, and resolve dofs on that entitie
  PetscErrorCode build_finite_elements(int verb = -1);
  PetscErrorCode build_finite_elements(const BitRefLevel &bit,int verb = -1);
  PetscErrorCode build_finite_elements(const boost::shared_ptr<FiniteElement> fe,const Range *ents_ptr = NULL,int verb = -1);
  PetscErrorCode build_finite_elements(const string fe_name,const Range *ents_ptr = NULL,int verb = -1);
  PetscErrorCode clear_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_finite_elements(const std::string &name,const Range &ents,int verb = -1);
  PetscErrorCode resolve_shared_ents(const Problem *problem_ptr,const std::string &fe_name,int verb = -1);
  PetscErrorCode resolve_shared_ents(const std::string &name,const std::string &fe_name,int verb = -1);
  PetscErrorCode get_problem_elements_layout(
    const std::string &name,const std::string &fe_name,PetscLayout *layout,int verb = -1
  );


  //entFEAdjacencies
  PetscErrorCode build_adjacencies(const Range &ents,int verb = -1);
  PetscErrorCode build_adjacencies(const BitRefLevel &bit,int verb = -1);
  PetscErrorCode build_adjacencies(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_adjacencies_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_adjacencies_entities(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_adjacencies_finite_elements(const std::string &name,const Range &ents,int verb = -1);
  PetscErrorCode clear_adjacencies_entities(const std::string &name,const Range &ents,int verb = -1);

  PetscErrorCode list_adjacencies() const;

  //problem building
  DEPRECATED PetscErrorCode build_problem_on_distributed_mesh(int verb = -1);
  DEPRECATED PetscErrorCode build_problem_on_distributed_mesh(
    const std::string &name,const bool square_matrix = true,int verb = -1
  );
  DEPRECATED PetscErrorCode build_problem_on_distributed_mesh(
    Problem *problem_ptr,const bool square_matrix = true,int verb = -1
  );

  DEPRECATED PetscErrorCode partition_mesh(
    const Range &ents,const int dim,const int adj_dim,const int n_parts,int verb = -1
  );
  DEPRECATED PetscErrorCode build_problem(const std::string &name,const bool square_matrix,int verb = -1);
  DEPRECATED PetscErrorCode build_problem(Problem *problem_ptr,const bool square_matrix,int verb = -1);
  DEPRECATED PetscErrorCode build_problems(int verb = -1);
  DEPRECATED PetscErrorCode partition_simple_problem(const std::string &name,int verb = -1);
  DEPRECATED PetscErrorCode partition_problem(const std::string &name,int verb = -1);
  DEPRECATED PetscErrorCode partition_compose_problem(
    const std::string &name,
    const std::string &problem_for_rows,
    const bool copy_rows,
    const std::string &problem_for_cols,
    const bool copy_cols,
    int verb = -1
  );
  DEPRECATED PetscErrorCode build_sub_problem(
    const std::string &out_name,
    const std::vector<std::string> &fields_row,
    const std::vector<std::string> &fields_col,
    const std::string &main_problem,
    const bool square_matrix = true,
    int verb = -1
  );
  DEPRECATED PetscErrorCode printPartitionedProblem(const Problem *problem_ptr,int verb = -1);
  DEPRECATED PetscErrorCode debugPartitionedProblem(const Problem *problem_ptr,int verb = -1);
  DEPRECATED PetscErrorCode partition_ghost_dofs(const std::string &name,int verb = -1);
  DEPRECATED PetscErrorCode partition_finite_elements(
    const std::string &name,
    bool part_from_moab = false,
    int low_proc = -1,
    int hi_proc = -1,
    int verb = -1
  );
  PetscErrorCode partition_check_matrix_fill_in(
    const std::string &problem_neme,int row,int col,int verb
  );

  ///save meshsets
  PetscErrorCode get_problem_finite_elements_entities(const std::string &name,const std::string &fe_name,const EntityHandle meshset);

  //vector and matrices
  PetscErrorCode MatCreateMPIAIJWithArrays(const std::string &name,Mat *Aij,int verb = -1);
  PetscErrorCode MatCreateMPIAdj_with_Idx_mi_tag(const std::string &name,Mat *Adj,int verb = -1);

  PetscErrorCode MatCreateSeqAIJWithArrays(const std::string &name,Mat *Aij,PetscInt **i,PetscInt **j,PetscScalar **v,int verb = -1);

  PetscErrorCode VecCreateSeq(const std::string &name,RowColData rc,Vec *V) const;
  PetscErrorCode VecCreateGhost(const std::string &name,RowColData rc,Vec *V) const;

  PetscErrorCode set_local_ghost_vector(const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) const;
  PetscErrorCode set_local_ghost_vector(const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) const;
  PetscErrorCode set_global_ghost_vector(const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) const;
  PetscErrorCode set_global_ghost_vector(const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) const;

  /// get IS for order
  PetscErrorCode ISCreateProblemOrder(
    const std::string &problem,RowColData rc,int min_order,int max_order,IS *is,int verb = -1
  ) const;
  /// get IS for field and rank
  PetscErrorCode ISCreateProblemFieldAndRank(
    const std::string &problem,
    RowColData rc,
    const std::string &field,
    int min_coeff_idx,
    int max_coeff_idx,
    IS *is,
    int verb = -1
  ) const;

  //scatter from problem filed to other problem field
  PetscErrorCode ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    std::vector<int> &idx,std::vector<int> &idy,int verb = -1
  ) const;
  PetscErrorCode ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    IS *ix,IS *iy,int verb = -1
  ) const;
  PetscErrorCode VecScatterCreate(
    Vec xin,const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    Vec yin,const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    VecScatter *newctx,int verb = -1
  ) const;

  //scatter from problem to other problem
  PetscErrorCode ISCreateFromProblemToOtherProblem(
    const std::string &x_problem,
    RowColData x_rc,
    const std::string &y_problem,
    RowColData y_rc,
    std::vector<int> &idx,
    std::vector<int> &idy,
    int verb = -1
  ) const;
  PetscErrorCode ISCreateFromProblemToOtherProblem(
    const std::string &x_problem,
    RowColData x_rc,
    const std::string &y_problem,
    RowColData y_rc,
    IS *ix,
    IS *iy,
    int verb = -1
  ) const;
  PetscErrorCode VecScatterCreate(
    Vec xin,
    const std::string &x_problem,
    RowColData x_rc,
    Vec yin,
    const std::string &y_problem,
    RowColData y_rc,
    VecScatter *newctx,
    int verb = -1
  ) const;

  //local
  PetscErrorCode set_other_local_ghost_vector(
    const Problem *problem_ptr,
    const std::string& fiel_name,
    const std::string& cpy_field_name,
    RowColData rc,Vec V,
    InsertMode mode,ScatterMode scatter_mode,int verb = -1
  );
  PetscErrorCode set_other_local_ghost_vector(
    const std::string &name,
    const std::string& fiel_name,
    const std::string& cpy_field_name,
    RowColData rc,Vec V,
    InsertMode mode,ScatterMode scatter_mode,int verb = -1
  );
  //global
  PetscErrorCode set_other_global_ghost_vector(
    const Problem *problem_ptr,
    const std::string& fiel_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb = -1
  );
  PetscErrorCode set_other_global_ghost_vector(
    const std::string &name,
    const std::string& fiel_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb = -1
  );

  //loops
  PetscErrorCode problem_basic_method_preProcess(const Problem *problem_ptr,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_preProcess(const std::string &problem_name,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_postProcess(const Problem *problem_ptr,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_postProcess(const std::string &problem_name,BasicMethod &method,int verb = -1);

  PetscErrorCode loop_finite_elements(
    const Problem *problem_ptr,const std::string &fe_name,FEMethod &method,
    int lower_rank,int upper_rank,MoFEMTypes bh = MF_EXIST,int verb = -1
  );
  PetscErrorCode loop_finite_elements(
    const std::string &problem_name,const std::string &fe_name,FEMethod &method,
    int lower_rank,int upper_rank,MoFEMTypes bh = MF_EXIST,int verb = -1
  );
  PetscErrorCode loop_finite_elements(
    const std::string &problem_name,const std::string &fe_name,FEMethod &method,
    MoFEMTypes bh = MF_EXIST,int verb = -1
  );

  PetscErrorCode loop_dofs(
    const Problem *problem_ptr,const std::string &field_name,RowColData rc,
    EntMethod &method,int lower_rank,int upper_rank,int verb = -1
  );
  PetscErrorCode loop_dofs(
    const std::string &problem_name,const std::string &field_name,RowColData rc,
    EntMethod &method,int lower_rank,int upper_rank,int verb = -1
  );
  PetscErrorCode loop_dofs(
    const std::string &problem_name,const std::string &field_name,RowColData rc,
    EntMethod &method,int verb = -1
  );
  PetscErrorCode loop_dofs(
    const std::string &field_name,EntMethod &method,int verb = -1
  );

  //get multi_index form database
  PetscErrorCode get_fields(const Field_multiIndex **fields_ptr) const;
  PetscErrorCode get_ref_ents(const RefEntity_multiIndex **refined_entities_ptr) const;
  PetscErrorCode get_ref_finite_elements(const RefElement_multiIndex **refined_finite_elements_ptr) const;
  PetscErrorCode get_finite_elements(const FiniteElement_multiIndex **fe_ptr) const;
  PetscErrorCode get_ents_finite_elements(const EntFiniteElement_multiIndex **fe_ent_ptr) const;
  PetscErrorCode get_field_ents(const FieldEntity_multiIndex **field_ents) const;
  PetscErrorCode get_dofs(const DofEntity_multiIndex **dofs_ptr) const ;
  PetscErrorCode get_problem(const std::string &problem_name,const Problem **problem_ptr) const;
  PetscErrorCode get_problems(const Problem_multiIndex **problems_ptr) const;


  FieldEntityByFieldName::iterator get_ent_moabfield_by_name_begin(const std::string &field_name) const;
  FieldEntityByFieldName::iterator get_ent_moabfield_by_name_end(const std::string &field_name) const;

  DofEntityByFieldName::iterator get_dofs_by_name_begin(const std::string &field_name) const;
  DofEntityByFieldName::iterator get_dofs_by_name_end(const std::string &field_name) const;
  DofEntityByNameAndEnt::iterator get_dofs_by_name_and_ent_begin(const std::string &field_name,const EntityHandle ent) const;
  DofEntityByNameAndEnt::iterator get_dofs_by_name_and_ent_end(const std::string &field_name,const EntityHandle ent) const;
  DofEntityByNameAndType::iterator get_dofs_by_name_and_type_begin(const std::string &field_name,const EntityType type) const;
  DofEntityByNameAndType::iterator get_dofs_by_name_and_type_end(const std::string &field_name,const EntityType ent) const;

  EntFiniteElementbyName::iterator get_fe_by_name_begin(const std::string &fe_name) const;
  EntFiniteElementbyName::iterator get_fe_by_name_end(const std::string &fe_name) const;

  //Copy field values to another field
  PetscErrorCode field_axpy(
    const double alpha,const std::string& fiel_name_x,const std::string& field_name_y,
    bool error_if_missing = false,bool creat_if_missing = false
  );
  PetscErrorCode field_scale(const double alpha,const std::string& fiel_name);
  PetscErrorCode set_field(const double val,const EntityType type,const std::string& fiel_name);
  PetscErrorCode set_field(const double val,const EntityType type,const Range &ents,const std::string& field_name);

  //Get adjacencies
  PetscErrorCode get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const;
  PetscErrorCode get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const;
  PetscErrorCode get_adjacencies(
    const Problem *problem_ptr,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type = moab::Interface::INTERSECT,
    const int verb = 0
  ) const;
  PetscErrorCode get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type = moab::Interface::INTERSECT,
    const int verb = 0
  ) const;

  // //Coordinate systems
  // DEPRECATED PetscErrorCode add_coordinate_system(const int cs_dim[],const std::string name);
  // DEPRECATED PetscErrorCode set_field_coordinate_system(const std::string field_name,const std::string cs_name);

  //Petsc Logs
  PetscLogEvent USER_EVENT_preProcess;
  PetscLogEvent USER_EVENT_operator;
  PetscLogEvent USER_EVENT_postProcess;
  PetscLogEvent USER_EVENT_createMat;
  PetscLogEvent USER_EVENT_buildProblem;

  // size and rank of communicator

  mutable MPI_Comm cOmm;      ///< MoFEM communicator
  mutable ParallelComm *pComm; ///< MOAB communicator struture

  int sIze; ///< MoFEM communicator size
  int rAnk; ///< MOFEM communicator rank

  inline MPI_Comm& get_comm() const { return cOmm; }
  inline int get_comm_size() const { return sIze; }
  inline int get_comm_rank() const { return rAnk; }

  DEPRECATED inline int getCommSize() const { return sIze; }
  DEPRECATED inline int getCommRank() const { return rAnk; }

  private:

  static bool isGloballyInitialised;

};

}

#endif // __MOABFIELD_CORE_HPP__
