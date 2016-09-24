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

/** \brief Core Interface class
 *  \ingroup mofem

  This class is not used directly by the user. For internal use only.   It is
  database with basic functions to access data. Abstraction of this is MoFEM
  Interface structure.

  It is deign to hide come complexities for users and allow low development
  without interfering with users modules programmer work.

 */
struct Core: public Interface, MeshRefinment, PrismInterface, SeriesRecorder {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);
  PetscErrorCode query_interface_type(const std::type_info& iface_type,void*& ptr) const;

  moab::Interface& moab;
  MPI_Comm comm;

  Core(moab::Interface& _moab,MPI_Comm _comm = PETSC_COMM_WORLD,int _verbose = 1);
  ~Core();

  Tag get_th_RefParentHandle() { return th_RefParentHandle; }
  Tag get_th_RefBitLevel() { return th_RefBitLevel; }

  protected:

  mutable boost::ptr_map<unsigned long,UnknownInterface *> iFaces;

  //Database
  ErrorCode rval;
  PetscErrorCode ierr;

  //Data and low level methods
  Tag th_Part;  ///< Tag for partition number
  Tag th_RefParentHandle,th_RefBitLevel,th_RefBitLevel_Mask,th_RefBitEdge,th_RefFEMeshset;
  Tag th_FieldId,th_FieldName,th_FieldName_DataNamePrefix,th_FieldSpace,th_FieldBase;
  Tag th_FEId,th_FEName;
  Tag th_FEIdCol,th_FEIdRow,th_FEIdData;
  Tag th_ProblemId,th_ProblemName,th_ProblemFEId;
  Tag th_ProblemNbDofsRow,th_ProblemNbDofsCol;
  Tag th_ProblemLocalNbDofRow,th_ProblemGhostNbDofRow;
  Tag th_ProblemLocalNbDofCol,th_ProblemGhostNbDofCol;
  Tag th_ProblemShift,th_FieldShift,th_FEShift;
  Tag th_ElemType;                    ///< Needed for VTK files
  Tag th_SeriesName;                  ///< Recorded series name

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

  //ref
  RefEntity_multiIndex refinedEntities;		       ///< refined entities
  RefElement_multiIndex refinedFiniteElements;	 ///< refined elements
  //field
  Field_multiIndex fIelds;			           ///< field
  MoFEMEntity_multiIndex entsFields;			 ///< entities on field
  DofEntity_multiIndex dofsField;		       ///< dofs on fiels
  //finite element
  FiniteElement_multiIndex finiteElements;		        ///< finite elements
  EntFiniteElement_multiIndex entsFiniteElements;			///< finite element entities
  //entFEAdjacencies
  MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex entFEAdjacencies;	///< adjacencies of elements to dofs
  //pRoblems
  MoFEMProblem_multiIndex pRoblems;					 ///< problems
  //cubit
  // CubitMeshSet_multiIndex cubitMeshsets;	   ///< cubit meshsets
  //series
  Series_multiIndex sEries;							///< recorded series
  SeriesStep_multiIndex seriesSteps;						///< recorded series steps

  //safety nets
  Tag th_MoFEMBuild;
  int *buildMoFEM;

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

  //core methods
  PetscErrorCode getTags(int verb = -1);
  PetscErrorCode clearMap();
  BitFieldId getFieldShift();
  BitFEId getFEShift();
  BitProblemId getProblemShift();
  PetscErrorCode initialiseDatabseInformationFromMesh(int verb = -1);

  //moab interface
  moab::Interface& get_moab();
  const moab::Interface& get_moab() const;

  //communicator MoFEM
  MPI_Comm get_comm() const;

  //add prims element
  PetscErrorCode addPrismToDatabase(const EntityHandle prism,int verb = -1);

  //MeshRefinemnt
  PetscErrorCode add_verices_in_the_middel_of_edges(
    const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1);
  PetscErrorCode add_verices_in_the_middel_of_edges(const Range &edges,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface = false);
  PetscErrorCode refine_TET(const Range &test,const BitRefLevel &bit,const bool respect_interface = false);
  PetscErrorCode refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1);

  //SeriesRecorder
  //add/delete series
  PetscErrorCode add_series_recorder(const std::string& series_name);
  PetscErrorCode delete_recorder_series(const std::string& series_name);
  //initialize/finalize recording
  PetscErrorCode initialize_series_recorder(const std::string& serie_name);
  PetscErrorCode finalize_series_recorder(const std::string& serie_name);
  //start recording
  PetscErrorCode record_begin(const std::string& serie_name);
  //recording functions
  PetscErrorCode record_problem(const std::string& serie_name,const MoFEMProblem *problemPtr,RowColData rc);
  PetscErrorCode record_problem(const std::string& serie_name,const std::string& problem_name,RowColData rc);
  PetscErrorCode record_field(const std::string& serie_name,const std::string& field_name,const BitRefLevel &bit,const BitRefLevel &mask);
  //end recording
  PetscErrorCode record_end(const std::string& serie_name,double time = 0);
  PetscErrorCode print_series_steps();
  bool check_series(const std::string& name) const;
  //get data back
  PetscErrorCode load_series_data(const std::string& serie_name,const int step_number);

  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator get_series_steps_byName_begin(const std::string& name);
  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator get_series_steps_byName_end(const std::string& name);

  //FiedlInterface

  //check consistency
  PetscErrorCode check_number_of_ents_in_ents_field(const std::string& name) const;
  PetscErrorCode check_number_of_ents_in_ents_field() const;
  PetscErrorCode check_number_of_ents_in_ents_finite_element(const std::string& name) const;
  PetscErrorCode check_number_of_ents_in_ents_finite_element() const;

  PetscErrorCode clear_database(int verb  = -1);
  PetscErrorCode rebuild_database(int verb = -1);

  PetscErrorCode get_msId_3dENTS_sides(
    const int msId,
    const CubitBCType cubit_bc_type,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1
  );
  PetscErrorCode get_msId_3dENTS_sides(
    const EntityHandle SIDESET,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1
  );
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const int msId,const CubitBCType cubit_bc_type,
    const bool add_iterfece_entities,const bool recursive = false,int verb = -1
  );
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const EntityHandle SIDESET,const bool add_iterfece_entities,const bool recursive = false,int verb = -1
  );
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const BitRefLevel &inheret_from_bit_level,const BitRefLevel &inheret_from_bit_level_mask,
    const EntityHandle SIDESET,const bool add_iterfece_entities,const bool recursive = false,int verb = -1
  );

  //cubit meshsets

  MeshsetsManager* meshsetsManagerPtr;
  MeshsetsManager* get_meshsets_manager_ptr() { return meshsetsManagerPtr; }
  const MeshsetsManager* get_meshsets_manager_ptr() const { return meshsetsManagerPtr; }
  MeshsetsManager& get_meshsets_manager() { return *meshsetsManagerPtr; }
  const MeshsetsManager& get_meshsets_manager() const { return *meshsetsManagerPtr; }

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
  DEPRECATED PetscErrorCode get_cubit_msId(const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr);
  DEPRECATED PetscErrorCode get_cubit_msId_entities_by_dimension(const int ms_id,const CubitBCType cubit_bc_type, const int dimension,Range &entities,const bool recursive = false);
  DEPRECATED PetscErrorCode get_cubit_msId_entities_by_dimension(const int ms_id,const CubitBCType cubit_bc_type, Range &entities,const bool recursive = false);
  DEPRECATED PetscErrorCode get_cubit_msId_entities_by_dimension(const int ms_id,const unsigned int cubit_bc_type, const int dimension,Range &entities,const bool recursive = false);
  DEPRECATED PetscErrorCode get_cubit_msId_entities_by_dimension(const int ms_id,const unsigned int cubit_bc_type, Range &entities,const bool recursive = false);
  DEPRECATED PetscErrorCode get_cubit_msId_meshset(const int ms_id,const unsigned int cubit_bc_type,EntityHandle &meshset);
  DEPRECATED PetscErrorCode get_cubit_meshsets(const unsigned int cubit_bc_type,Range &meshsets);

  DEPRECATED PetscErrorCode print_cubit_displacement_set() const;
  DEPRECATED PetscErrorCode print_cubit_pressure_set() const;
  DEPRECATED PetscErrorCode print_cubit_force_set() const;
  DEPRECATED PetscErrorCode print_cubit_temperature() const;
  DEPRECATED PetscErrorCode print_cubit_heat_flux_set() const;
  DEPRECATED PetscErrorCode print_cubit_materials_set() const;

  //refine
  PetscErrorCode seed_finite_elements(const Range &entities,int verb = -1);
  PetscErrorCode seed_finite_elements(const EntityHandle meshset,int verb = -1);
  PetscErrorCode seed_ref_level(const Range &ents,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode seed_ref_level_2D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb = -1);
  PetscErrorCode get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb = -1);
  PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset);
  PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents);
  PetscErrorCode add_ref_level_to_entities(const BitRefLevel &bit,Range &ents);
  PetscErrorCode set_ref_level_to_entities(const BitRefLevel &bit,Range &ents);
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
   * @param  base           Approximation base AINSWORTH_COLE_BASE, BERNSTEIN_BEZIER_BASE ...
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
  PetscErrorCode dofs_NoField(const BitFieldId id,std::map<EntityType,int> &dof_counter,int verb = -1);
  PetscErrorCode dofs_L2H1HcurlHdiv(
    const BitFieldId id,std::map<EntityType,int> &dof_counter,std::map<EntityType,int> &inactive_dof_counter,int verb = -1
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
  PetscErrorCode add_ents_to_finite_element_by_VERTICEs(const Range& vert,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_VERTICEs(const Range& vert,const std::string &name);
  PetscErrorCode add_ents_to_finite_element_by_EDGEs(const Range& vert,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_EDGEs(const Range& vert,const std::string &name);
  PetscErrorCode add_ents_to_finite_element_by_TRIs(const Range& tris,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_TRIs(const Range& tris,const std::string &name);
  PetscErrorCode add_ents_to_finite_element_by_TRIs(const EntityHandle meshset,const std::string &name,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const Range& tets,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const Range& tets,const std::string &name);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const BitFEId id,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const std::string &name,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_by_PRISMs(const Range& prims,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_PRISMs(const Range& prims,const std::string &name);
  PetscErrorCode add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const BitFEId id,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const std::string &name,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const std::string& name,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const std::string &name,EntityType type,int verb = -1);
  PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const std::string &name,EntityType type,int verb = -1);
  PetscErrorCode remove_ents_from_finite_element_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode remove_ents_from_finite_element(const std::string &name,const EntityHandle meshset,const EntityType type,int verb = -1);
  PetscErrorCode remove_ents_from_finite_element(const std::string &name,const Range &ents,int verb = -1);
  PetscErrorCode delete_finite_element(const std::string name,int verb = -1);

  //other auxiliary functions for finite element
  BitFEId getBitFEId(const std::string& name) const;
  std::string getBitFEId_name(const BitFEId id) const;
  EntityHandle get_finite_element_meshset(const BitFEId id) const;
  EntityHandle get_finite_element_meshset(const std::string& name) const;
  PetscErrorCode list_finite_elements() const;

  //problem
  PetscErrorCode add_problem(const BitProblemId id,const std::string& name);
  PetscErrorCode add_problem(const std::string& name,enum MoFEMTypes bh = MF_EXCL,int verb = -1);
  PetscErrorCode delete_problem(const std::string name);
  PetscErrorCode modify_problem_add_finite_element(const std::string &name_problem,const std::string &MoFEMFiniteElement_name);
  PetscErrorCode modify_problem_unset_finite_element(const std::string &name_problem,const std::string &MoFEMFiniteElement_name);
  PetscErrorCode modify_problem_ref_level_add_bit(const std::string &name_problem,const BitRefLevel &bit);
  PetscErrorCode modify_problem_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit);
  PetscErrorCode modify_problem_dof_mask_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit);
  BitProblemId get_BitProblemId(const std::string& name) const;
  PetscErrorCode list_problem() const;

  ///add entity EntFe to finite element data databse and resolve dofs on that entity
  //loop over all finite elements, resolve its meshsets, and resolve dofs on that entitie
  PetscErrorCode build_finite_elements(int verb = -1);
  PetscErrorCode build_finite_elements(const BitRefLevel &bit,int verb = -1);
  PetscErrorCode build_finite_elements(const boost::shared_ptr<FiniteElement> fe,const Range *ents_ptr = NULL,int verb = -1);
  PetscErrorCode build_finite_elements(const string fe_name,const Range *ents_ptr = NULL,int verb = -1);
  PetscErrorCode clear_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_finite_elements(const std::string &name,const Range &ents,int verb = -1);

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
  PetscErrorCode build_problem_on_partitioned_mesh(MoFEMProblem *problem_ptr,bool square_matrix = true,int verb = -1);
  PetscErrorCode build_problem_on_distributed_mesh(int verb = -1);
  PetscErrorCode build_problem_on_distributed_mesh(const std::string &name,bool square_matrix = true,int verb = -1);
  PetscErrorCode build_problem_on_distributed_mesh(MoFEMProblem *problem_ptr,bool square_matrix = true,int verb = -1);
  PetscErrorCode partition_mesh(Range &ents,int dim,int adj_dim,int n_parts,int verb = -1);
  PetscErrorCode build_problem(const std::string &name,int verb = -1);
  PetscErrorCode clear_problem(const std::string &name,int verb = -1);
  PetscErrorCode build_problem(MoFEMProblem *problem_ptr,int verb = -1);
  PetscErrorCode build_problems(int verb = -1);
  PetscErrorCode clear_problems(int verb = -1);
  PetscErrorCode partition_simple_problem(const std::string &name,int verb = -1);
  PetscErrorCode partition_problem(const std::string &name,int verb = -1);
  PetscErrorCode partition_compose_problem(const std::string &name,const std::string &problem_for_rows,bool copy_rows,const std::string &problem_for_cols,bool copy_cols,int verb = -1);
  PetscErrorCode partition_ghost_dofs(const std::string &name,int verb = -1);
  PetscErrorCode partition_finite_elements(const std::string &name,bool part_from_moab = false,int low_proc = -1,int hi_proc = -1,int verb = -1);
  PetscErrorCode partition_check_matrix_fill_in(const std::string &problem_neme,int row,int col,int verb);
  PetscErrorCode printPartitionedProblem(const MoFEMProblem *problem_ptr,int verb = -1);
  PetscErrorCode debugPartitionedProblem(const MoFEMProblem *problem_ptr,int verb = -1);
  PetscErrorCode resolve_shared_ents(const MoFEMProblem *problem_ptr,const std::string &fe_name,int verb = -1);
  PetscErrorCode resolve_shared_ents(const std::string &name,const std::string &fe_name,int verb = -1);
  PetscErrorCode get_problem_elements_layout(
    const std::string &name,const std::string &fe_name,PetscLayout *layout,int verb = -1
  );

  ///save meshsets
  PetscErrorCode get_problem_finite_elements_entities(const std::string &name,const std::string &fe_name,const EntityHandle meshset);

  //vector and matrices
  PetscErrorCode MatCreateMPIAIJWithArrays(const std::string &name,Mat *Aij,int verb = -1);
  PetscErrorCode MatCreateSeqAIJWithArrays(const std::string &name,Mat *Aij,PetscInt **i,PetscInt **j,PetscScalar **v,int verb = -1);

  PetscErrorCode VecCreateSeq(const std::string &name,RowColData rc,Vec *V) const;
  PetscErrorCode VecCreateGhost(const std::string &name,RowColData rc,Vec *V) const;

  PetscErrorCode set_local_ghost_vector(const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) const;
  PetscErrorCode set_local_ghost_vector(const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) const;
  PetscErrorCode set_global_ghost_vector(const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) const;
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
    const MoFEMProblem *problem_ptr,const std::string& fiel_name,const std::string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1);
  PetscErrorCode set_other_local_ghost_vector(
    const std::string &name,const std::string& fiel_name,const std::string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1);
  //global
  PetscErrorCode set_other_global_ghost_vector(
    const MoFEMProblem *problem_ptr,
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
  PetscErrorCode problem_basic_method_preProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_preProcess(const std::string &problem_name,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_postProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_postProcess(const std::string &problem_name,BasicMethod &method,int verb = -1);

  PetscErrorCode loop_finite_elements(
    const MoFEMProblem *problem_ptr,const std::string &fe_name,FEMethod &method,
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
    const MoFEMProblem *problem_ptr,const std::string &field_name,RowColData rc,
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
  PetscErrorCode get_problem(const std::string &problem_name,const MoFEMProblem **problem_ptr) const;
  PetscErrorCode get_field_ents(const MoFEMEntity_multiIndex **field_ents) const;
  PetscErrorCode get_dofs(const DofEntity_multiIndex **dofs_ptr) const ;
  PetscErrorCode get_finite_elements(const FiniteElement_multiIndex **finiteElements_ptr) const;

  MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_ent_moabfield_by_name_begin(const std::string &field_name) const;
  MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_ent_moabfield_by_name_end(const std::string &field_name) const;

  DofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_begin(const std::string &field_name) const;
  DofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_end(const std::string &field_name) const;
  DofEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator get_dofs_by_name_and_ent_begin(const std::string &field_name,const EntityHandle ent) const;
  DofEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator get_dofs_by_name_and_ent_end(const std::string &field_name,const EntityHandle ent) const;
  DofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator get_dofs_by_name_and_type_begin(const std::string &field_name,const EntityType type) const;
  DofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator get_dofs_by_name_and_type_end(const std::string &field_name,const EntityType ent) const;

  EntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator get_fe_by_name_begin(const std::string &fe_name) const;
  EntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator get_fe_by_name_end(const std::string &fe_name) const;

  //Copy field values to another field
  PetscErrorCode field_axpy(const double alpha,const std::string& fiel_name_x,const std::string& field_name_y,bool error_if_missing = false,bool creat_if_missing = false);
  PetscErrorCode field_scale(const double alpha,const std::string& fiel_name);
  PetscErrorCode set_field(const double val,const EntityType type,const std::string& fiel_name);
  PetscErrorCode set_field(const double val,const EntityType type,const Range &ents,const std::string& field_name);

  //Get adjacencies
  PetscErrorCode get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const;
  PetscErrorCode get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const;
  PetscErrorCode get_adjacencies(
    const MoFEMProblem *problem_ptr,
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
  int sIze,rAnk;
  int getCommSize() const { return sIze; }
  int getCommRank() const { return rAnk; }

  private:

  static bool isGloballyInitialised;

};

}

#endif // __MOABFIELD_CORE_HPP__
