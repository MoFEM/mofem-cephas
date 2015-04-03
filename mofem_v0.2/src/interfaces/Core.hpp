/** \file Core.hpp
 * \brief Core FieldInterface class for user interface
 * 
 * Low level data structures not used directly by user
 *
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

#include "FieldInterface.hpp"
#include "MeshRefinment.hpp"
#include "PrismInterface.hpp"
#include "SeriesRecorder.hpp"

namespace MoFEM {

/** \brief Core FieldInterface class
 *  \ingroup mofem
 
  This class is not used directly by the user. It is database with basic
  functions to access data. Abstraction of this is MoFEM Interface structure.

  It is deign to hide come complexities for users and allow low development
  without interfering with users modules programmer work.

 */
struct Core: 
  public FieldInterface, MeshRefinment, PrismInterface, SeriesRecorder {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface);
  PetscErrorCode query_interface_type(const std::type_info& iface_type, void*& ptr);

  Interface& moab;
  MPI_Comm comm;

  Core(Interface& _moab,MPI_Comm _comm = PETSC_COMM_WORLD,int _verbose = 1);
  ~Core();

  Tag get_th_RefParentHandle() { return th_RefParentHandle; }
  Tag get_th_RefBitLevel() { return th_RefBitLevel; }

  protected:

  boost::ptr_map<unsigned long,FieldUnknownInterface *> iFaces;

  //Database
  ErrorCode rval;
  PetscErrorCode ierr;
  //Data and low level methods 
  Tag th_Part;
  Tag th_RefType,th_RefParentHandle,th_RefBitLevel,th_RefBitLevel_Mask,th_RefBitEdge,th_RefFEMeshset;
  Tag th_FieldId,th_FieldName,th_FieldName_DataNamePrefix,th_FieldSpace;
  Tag th_FEId,th_FEName;
  Tag th_FEIdCol,th_FEIdRow,th_FEIdData;
  Tag th_ProblemId,th_ProblemName,th_ProblemFEId;
  Tag th_ProblemNbDofsRow,th_ProblemNbDofsCol;
  Tag th_ProblemLocalNbDofRow,th_ProblemGhostNbDofRow;
  Tag th_ProblemLocalNbDofCol,th_ProblemGhostNbDofCol;
  Tag th_ProblemShift,th_FieldShift,th_FEShift;
  Tag nsTag,ssTag,nsTag_data,ssTag_data,bhTag,bhTag_header;
  Tag th_ElemType;
  Tag th_SeriesName;

  int *fShift,*feShift,*pShift;
  int verbose;

  //ref
  RefMoFEMEntity_multiIndex refinedEntities;		///< refined entities
  RefMoFEMElement_multiIndex refinedFiniteElements;	///< refined elements
  //field
  MoFEMField_multiIndex moabFields;			///< field
  MoFEMEntity_multiIndex entsMoabField;			///< entities on field
  DofMoFEMEntity_multiIndex dofsMoabField;		///< dofs on fiels
  //finite element
  MoFEMFiniteElement_multiIndex finiteElements;		///< finite elements
  EntMoFEMFiniteElement_multiIndex finiteElementsMoFEMEnts;			///< finite element entities
  //entFEAdjacencies
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex entFEAdjacencies;	///< adjacencies of elements to dofs
  //moFEMProblems
  MoFEMProblem_multiIndex moFEMProblems;					///< problems	
  //cubit
  CubitMeshSet_multiIndex cubitMeshsets;					///< cubit meshsets
  //series
  Series_multiIndex sEries;							///< recorded series
  SeriesStep_multiIndex seriesSteps;						///< recorded series steps

  //safety nets
  Tag th_MoFEMBuild;
  int *build_MoFEM;

  //core methods
  PetscErrorCode clearMap();
  BitFieldId getFieldShift();
  BitFEId getFEShift();
  BitProblemId getProblemShift();
  PetscErrorCode initialiseDatabseInformationFromMesh(int verb = -1);

  //moab interface
  Interface& get_moab();

  //communicator MoFEM
  MPI_Comm get_comm();

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
  PetscErrorCode add_series_recorder(const string& series_name);
  PetscErrorCode delete_recorder_series(const string& series_name);
  //initialize/finalize recording
  PetscErrorCode initialize_series_recorder(const string& serie_name);
  PetscErrorCode finalize_series_recorder(const string& serie_name);
  //start recording
  PetscErrorCode record_begin(const string& serie_name);
  //recording functions
  PetscErrorCode record_problem(const string& serie_name,const MoFEMProblem *problemPtr,RowColData rc);
  PetscErrorCode record_problem(const string& serie_name,const string& problem_name,RowColData rc);
  PetscErrorCode record_field(const string& serie_name,const string& field_name,const BitRefLevel &bit,const BitRefLevel &mask);
  //end recording
  PetscErrorCode record_end(const string& serie_name,double time = 0);
  PetscErrorCode print_series_steps();
  bool check_series(const string& name) const;
  //get data back
  PetscErrorCode load_series_data(const string& serie_name,const int step_number);

  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator get_series_steps_byName_begin(const string& name);
  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator get_series_steps_byName_end(const string& name);

  //PrismInrerface

  PetscErrorCode get_msId_3dENTS_sides(
    const int msId,
    const CubitBCType cubit_bc_type,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1);
  PetscErrorCode get_msId_3dENTS_sides(
    const EntityHandle SIDESET,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1);
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const int msId,const CubitBCType cubit_bc_type,
    const bool add_iterfece_entities,const bool recursive = false,int verb = -1);
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const EntityHandle SIDESET,const bool add_iterfece_entities,const bool recursive = false,int verb = -1);
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const BitRefLevel &inheret_from_bit_level,const BitRefLevel &inheret_from_bit_level_mask,
    const EntityHandle SIDESET,const bool add_iterfece_entities,const bool recursive = false,int verb = -1);


  //FiedlInterface

  //check consistency
  PetscErrorCode check_number_of_ents_in_ents_field(const string& name);
  PetscErrorCode check_number_of_ents_in_ents_field();
  PetscErrorCode check_number_of_ents_in_ents_finite_element(const string& name);
  PetscErrorCode check_number_of_ents_in_ents_finite_element();
  PetscErrorCode rebuild_database(int verb = -1);


  //cubit meshsets
  bool check_msId_meshset(const int msId,const CubitBCType cubit_bc_type);
  PetscErrorCode add_cubit_msId(const CubitBCType cubit_bc_type,const int msId);
  PetscErrorCode delete_cubit_msId(const CubitBCType cubit_bc_type,const int msId);
  PetscErrorCode get_cubit_msId(const int msId,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr);
  PetscErrorCode get_cubit_msId_entities_by_dimension(const int msId,const CubitBCType cubit_bc_type, const int dimension,Range &entities,const bool recursive = false);
  PetscErrorCode get_cubit_msId_entities_by_dimension(const int msId,const CubitBCType cubit_bc_type, Range &entities,const bool recursive = false);
  PetscErrorCode get_cubit_msId_entities_by_dimension(const int msId,const unsigned int cubit_bc_type, const int dimension,Range &entities,const bool recursive = false);
  PetscErrorCode get_cubit_msId_entities_by_dimension(const int msId,const unsigned int cubit_bc_type, Range &entities,const bool recursive = false);
  PetscErrorCode get_cubit_msId_meshset(const int msId,const unsigned int cubit_bc_type,EntityHandle &meshset);
  PetscErrorCode get_cubit_meshsets(const unsigned int cubit_bc_type,Range &meshsets);
  CubitMeshSet_multiIndex::iterator get_cubit_meshsets_begin() { return cubitMeshsets.begin(); }
  CubitMeshSet_multiIndex::iterator get_cubit_mesh_sets_end() { return cubitMeshsets.end(); }
  CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_cubit_meshsets_begin(const unsigned int cubit_bc_type) { 
    return cubitMeshsets.get<CubitMeshSets_mi_tag>().lower_bound(cubit_bc_type); 
  }
  CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_cubit_mesh_sets_end(const unsigned int cubit_bc_type) { 
    return cubitMeshsets.get<CubitMeshSets_mi_tag>().upper_bound(cubit_bc_type); 
  }
  CubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_begin(const unsigned int cubit_bc_type) { 
    return cubitMeshsets.get<CubitMeshSets_mask_meshset_mi_tag>().lower_bound(cubit_bc_type); 
  }
  CubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_end(const unsigned int cubit_bc_type) { 
    return cubitMeshsets.get<CubitMeshSets_mask_meshset_mi_tag>().upper_bound(cubit_bc_type); 
  }
  CubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator get_CubitMeshSets_byName_begin(const string& name) { 
    return cubitMeshsets.get<CubitMeshSets_name>().lower_bound(name); 
  }
  CubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator get_CubitMeshSets_byName_end(const string& name) { 
    return cubitMeshsets.get<CubitMeshSets_name>().upper_bound(name); 
  }

  template<class _CUBIT_BC_DATA_TYPE_>
  PetscErrorCode printCubitSet(_CUBIT_BC_DATA_TYPE_& data,unsigned long int type) {
    PetscFunctionBegin;
    try {
      FieldInterface& this_mField = *this;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(this_mField,type,it)) {
	ierr = it->get_cubit_bc_data_structure(data); CHKERRQ(ierr);
	ostringstream ss;
	ss << *it << endl;
	ss << data << endl;
	Range tets,tris,edges,nodes;
	rval = moab.get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
	rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
	rval = moab.get_entities_by_type(it->meshset,MBEDGE,edges,true); CHKERR_PETSC(rval);
	rval = moab.get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
	ss << "name "<< it->get_Cubit_name() << endl;
	ss << "msId "<< it->get_msId() << " nb. tets " << tets.size() << endl;
	ss << "msId "<< it->get_msId() << " nb. tris " << tris.size() << endl;
	ss << "msId "<< it->get_msId() << " nb. edges " << edges.size() << endl;
	ss << "msId "<< it->get_msId() << " nb. nodes " << nodes.size() << endl;
	ss << endl;
	PetscPrintf(comm,ss.str().c_str());
      }
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,1,msg);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode print_cubit_displacement_set() {
    PetscFunctionBegin;
    DisplacementCubitBcData mydata;
    ierr = printCubitSet(mydata,NODESET|mydata.type.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  PetscErrorCode print_cubit_pressure_set() {
    PetscFunctionBegin;
    PressureCubitBcData mydata;
    ierr = printCubitSet(mydata,SIDESET|mydata.type.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
   
  PetscErrorCode print_cubit_force_set() {
    PetscFunctionBegin;
    ForceCubitBcData mydata;
    ierr = printCubitSet(mydata,NODESET|mydata.type.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode print_cubit_temperature() {
    PetscFunctionBegin;
    TemperatureCubitBcData mydata;
    ierr = printCubitSet(mydata,NODESET|mydata.type.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
    
  PetscErrorCode print_cubit_heat_flux_set() {
    PetscFunctionBegin;
    HeatfluxCubitBcData mydata;
    ierr = printCubitSet(mydata,SIDESET|mydata.type.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode print_cubit_materials_set() {
    PetscFunctionBegin;
    FieldInterface& this_mField = *this;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(this_mField,BLOCKSET|MAT_ELASTICSET,it)) {
      Mat_Elastic data;
      ierr = it->get_attribute_data_structure(data); CHKERRQ(ierr);
      ostringstream ss;
      ss << *it << endl;
      ss << data;
      Range tets;
      rval = moab.get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
      ss << "MAT_ELATIC msId "<< it->get_msId() << " nb. tets " << tets.size() << endl;
      ss << endl;
      PetscPrintf(comm,ss.str().c_str());
    }
      
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(this_mField,BLOCKSET|MAT_THERMALSET,it)) {
        Mat_Thermal data;
        ierr = it->get_attribute_data_structure(data); CHKERRQ(ierr);
        ostringstream ss;
        ss << *it << endl;
        ss << data;
        PetscPrintf(comm,ss.str().c_str());
    }
    
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(this_mField,BLOCKSET|MAT_MOISTURESET,it)) {
      Mat_Moisture data;
      ierr = it->get_attribute_data_structure(data); CHKERRQ(ierr);
      ostringstream ss;
      ss << *it << endl;
      ss << data;
      PetscPrintf(comm,ss.str().c_str());
    }

    
    PetscFunctionReturn(0);
  }

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
    const bool recursive = false, int verb = -1);
  PetscErrorCode update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb = -1);
  PetscErrorCode update_field_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,int verb = -1);
  PetscErrorCode update_finite_element_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb = -1);

  //remove entities
  PetscErrorCode delete_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const bool remove_parent = false,int verb = -1);
  PetscErrorCode remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode delete_finite_elements_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode shift_left_bit_ref(const int shif,int verb = -1);
  PetscErrorCode shift_right_bit_ref(const int shift,int verb = -1);

  //synchronize entities
  PetscErrorCode synchronise_entities(Range &ent,int verb = -1);

  //field
  PetscErrorCode add_field(
    const string& name,const FieldSpace space,const ApproximationRank rank,enum MoFEMTypes bh = MF_EXCL,int verb = -1);
  PetscErrorCode add_ents_to_field_by_VERTICEs(const Range &nodes,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_VERTICEs(const Range &nodes,const string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_EDGEs(const EntityHandle meshset,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_EDGEs(const EntityHandle meshset,const string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TRIs(const EntityHandle meshset,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TRIs(const EntityHandle meshset,const string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TRIs(const Range &tris,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TRIs(const Range &tris,const string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TETs(const Range &tets,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const BitFieldId id,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TETs(const Range &tets,const string& name,int verb = -1);
  PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name,int verb = -1);
  PetscErrorCode remove_ents_from_field_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode remove_ents_from_field(const string& name,const EntityHandle meshset,const EntityType type,int verb = -1);
  PetscErrorCode remove_ents_from_field(const string& name,const Range &ents,int verb = -1);

  //set apprix oorder
  PetscErrorCode set_field_order(const Range &ents,const BitFieldId id,const ApproximationOrder order,int verb = -1);
  PetscErrorCode set_field_order(const EntityHandle meshset,const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb = -1);
  PetscErrorCode set_field_order(const Range &ents,const string& name,const ApproximationOrder order,int verb = -1);
  PetscErrorCode set_field_order(const EntityHandle meshset,const EntityType type,const string& name,const ApproximationOrder order,int verb = -1);
  PetscErrorCode set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb = -1);
  PetscErrorCode set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const string& name,const ApproximationOrder order,int verb = -1);

  //build fiels
  PetscErrorCode dofs_NoField(const BitFieldId id,map<EntityType,int> &dof_counter,int verb = -1);
  PetscErrorCode dofs_L2H1HcurlHdiv(const BitFieldId id,map<EntityType,int> &dof_counter,int verb = -1);
  PetscErrorCode build_fields(int verb = -1);
  PetscErrorCode clear_dofs_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_ents_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_dofs_fields(const string &name,const Range ents,int verb = -1);
  PetscErrorCode clear_ents_fields(const string &name,const Range enst,int verb = -1);

  //other auxiliary functions for fields
  PetscErrorCode list_dofs_by_field_name(const string &name) const;
  PetscErrorCode list_fields() const;

  BitFieldId get_BitFieldId(const string& name) const;
  string get_BitFieldId_name(const BitFieldId id) const;
  EntityHandle get_field_meshset(const BitFieldId id) const;
  EntityHandle get_field_meshset(const string& name) const;
    bool check_field(const string& name) const;
  const MoFEMField* get_field_structure(const string& name);

  //MoFEMFiniteElement
  PetscErrorCode add_finite_element(const string &MoFEMFiniteElement_name,enum MoFEMTypes bh = MF_EXCL);
  PetscErrorCode modify_finite_element_adjacency_table(const string &MoFEMFiniteElement_name,const EntityType type,ElementAdjacencyFunct function);
  PetscErrorCode modify_finite_element_add_field_data(const string &MoFEMFiniteElement_name,const string &name_filed);
  PetscErrorCode modify_finite_element_add_field_row(const string &MoFEMFiniteElement_name,const string &name_row);
  PetscErrorCode modify_finite_element_add_field_col(const string &MoFEMFiniteElement_name,const string &name_col);
  PetscErrorCode modify_finite_element_off_field_data(const string &MoFEMFiniteElement_name,const string &name_filed);
  PetscErrorCode modify_finite_element_off_field_row(const string &MoFEMFiniteElement_name,const string &name_row);
  PetscErrorCode modify_finite_element_off_field_col(const string &MoFEMFiniteElement_name,const string &name_col);
  PetscErrorCode add_ents_to_finite_element_by_VERTICEs(const Range& vert,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_VERTICEs(const Range& vert,const string &name);
  PetscErrorCode add_ents_to_finite_element_by_EDGEs(const Range& vert,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_EDGEs(const Range& vert,const string &name);
  PetscErrorCode add_ents_to_finite_element_by_TRIs(const Range& tris,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_TRIs(const Range& tris,const string &name);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const Range& tets,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const Range& tets,const string &name);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const BitFEId id,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const string &name,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_by_PRISMs(const Range& prims,const BitFEId id);
  PetscErrorCode add_ents_to_finite_element_by_PRISMs(const Range& prims,const string &name);
  PetscErrorCode add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const BitFEId id,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const string &name,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const string& name,const bool recursive = false);
  PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const string &name,EntityType type,int verb = -1);
  PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const string &name,EntityType type,int verb = -1);
  PetscErrorCode remove_ents_from_finite_element_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode remove_ents_from_finite_element(const string &name,const EntityHandle meshset,const EntityType type,int verb = -1);
  PetscErrorCode remove_ents_from_finite_element(const string &name,const Range &ents,int verb = -1);

  //other auxiliary functions for finite element
  BitFEId get_BitFEId(const string& name) const;
  string get_BitFEId_name(const BitFEId id) const;
  EntityHandle get_finite_element_meshset(const BitFEId id) const;
  EntityHandle get_finite_element_meshset(const string& name) const;
  PetscErrorCode list_finite_elements() const;

  //problem
  PetscErrorCode add_problem(const BitProblemId id,const string& name);
  PetscErrorCode add_problem(const string& name,enum MoFEMTypes bh = MF_EXCL,int verb = -1);
  PetscErrorCode modify_problem_add_finite_element(const string &name_problem,const string &MoFEMFiniteElement_name);
  PetscErrorCode modify_problem_ref_level_add_bit(const string &name_problem,const BitRefLevel &bit);
  PetscErrorCode modify_problem_ref_level_set_bit(const string &name_problem,const BitRefLevel &bit);
  PetscErrorCode modify_problem_dof_mask_ref_level_set_bit(const string &name_problem,const BitRefLevel &bit);
  BitProblemId get_BitProblemId(const string& name) const;
  PetscErrorCode list_problem() const;

  ///add entity EntFe to finite element data databse and resolve dofs on that entity
  //loop over all finite elements, resolve its meshsets, and resolve dofs on that entitie
  PetscErrorCode build_finite_element_data_dofs(EntMoFEMFiniteElement &ent_fe,int verb = -1);
  PetscErrorCode build_finite_element_uids_view(EntMoFEMFiniteElement &ent_fe,int verb = -1);
  PetscErrorCode build_finite_elements(int verb = -1);
  PetscErrorCode clear_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_finite_elements(const string &name,const Range &ents,int verb = -1);

  //entFEAdjacencies
  PetscErrorCode build_adjacencies(const Range &ents,int verb = -1);
  PetscErrorCode build_adjacencies(const BitRefLevel &bit,int verb = -1);
  PetscErrorCode build_adjacencies(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_adjacencies_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_adjacencies_entities(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_adjacencies_finite_elements(const string &name,const Range &ents,int verb = -1);
  PetscErrorCode clear_adjacencies_entities(const string &name,const Range &ents,int verb = -1);

  PetscErrorCode list_adjacencies() const;

  //problem building
  PetscErrorCode build_problem(const string &name,int verb = -1);
  PetscErrorCode build_problem(MoFEMProblem *problem_ptr,int verb = -1);
  PetscErrorCode build_problems(int verb = -1);
  PetscErrorCode clear_problems(int verb = -1);
  PetscErrorCode build_partitioned_problems(int verb = -1);
  PetscErrorCode build_partitioned_problem(const string &name,bool square_matrix = true,int verb = -1);
  PetscErrorCode build_partitioned_problem(MoFEMProblem *problem_ptr,bool square_matrix = true,int verb = -1);
  PetscErrorCode simple_partition_problem(const string &name,int verb = -1);
  PetscErrorCode partition_problem(const string &name,int verb = -1);
  PetscErrorCode compose_problem(const string &name,const string &problem_for_rows,bool copy_rows,const string &problem_for_cols,bool copy_cols,int verb = -1);
  PetscErrorCode block_problem(const string &name,const vector<string> block_problems,int verb = -1);
  PetscErrorCode partition_ghost_dofs(const string &name,int verb = -1);
  PetscErrorCode partition_finite_elements(
    const string &name,bool part_from_moab = false,int low_proc = -1,int hi_proc = -1,int verb = -1);
  PetscErrorCode partition_check_matrix_fill_in(const string &problem_neme,int verb);
  PetscErrorCode print_partitioned_problem(const MoFEMProblem *problem_ptr,int verb = -1);
  PetscErrorCode debug_partitioned_problem(const MoFEMProblem *problem_ptr,int verb = -1);

  ///save meshsets
  PetscErrorCode get_problem_finite_elements_entities(const string &name,const string &fe_name,const EntityHandle meshset);

  //vector and matrices 
  PetscErrorCode MatCreateMPIAIJWithArrays(const string &name,Mat *Aij,int verb = -1);
  PetscErrorCode MatCreateSeqAIJWithArrays(const string &name,Mat *Aij,PetscInt **i,PetscInt **j,PetscScalar **v,int verb = -1);

  PetscErrorCode VecCreateSeq(const string &name,RowColData rc,Vec *V);
  PetscErrorCode VecCreateGhost(const string &name,RowColData rc,Vec *V);
  PetscErrorCode set_local_ghost_vector(const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode);
  PetscErrorCode set_local_ghost_vector(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode);
  PetscErrorCode set_global_ghost_vector(const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode); 
  PetscErrorCode set_global_ghost_vector(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode);

  /// get IS for order
  PetscErrorCode ISCreateProblemOrder(const string &problem,RowColData rc,int min_order,int max_order,IS *is,int verb = -1);
  /// get IS for field and rank
  PetscErrorCode ISCreateProblemFieldAndRank(const string &problem,RowColData rc,const string &field,int min_rank,int max_rank,IS *is,int verb = -1);

  //scatter from problem filed to other problem field
  PetscErrorCode ISCreateFromProblemFieldToOtherProblemField(
    const string &x_problem,const string &x_field_name,RowColData x_rc,
    const string &y_problem,const string &y_field_name,RowColData y_rc,
    vector<int> &idx,vector<int> &idy,int verb = -1);
  PetscErrorCode ISCreateFromProblemFieldToOtherProblemField(
    const string &x_problem,const string &x_field_name,RowColData x_rc,
    const string &y_problem,const string &y_field_name,RowColData y_rc,
    IS *ix,IS *iy,int verb = -1);
  PetscErrorCode VecScatterCreate(
    Vec xin,const string &x_problem,const string &x_field_name,RowColData x_rc,
    Vec yin,const string &y_problem,const string &y_field_name,RowColData y_rc,VecScatter *newctx,int verb = -1);

  //scatter from problem to other problem
  PetscErrorCode ISCreateFromProblemToOtherProblem(
    const string &x_problem,RowColData x_rc,const string &y_problem,RowColData y_rc,vector<int> &idx,vector<int> &idy,int verb = -1);
  PetscErrorCode ISCreateFromProblemToOtherProblem(
    const string &x_problem,RowColData x_rc,const string &y_problem,RowColData y_rc,IS *ix,IS *iy,int verb = -1);
    PetscErrorCode VecScatterCreate(Vec xin,const string &x_problem,RowColData x_rc,Vec yin,const string &y_problem,RowColData y_rc,VecScatter *newctx,int verb = -1);
  //local
  PetscErrorCode set_other_local_ghost_vector(
    const MoFEMProblem *problem_ptr,const string& fiel_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1);
  PetscErrorCode set_other_local_ghost_vector(
    const string &name,const string& fiel_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1);
  //global
  PetscErrorCode set_other_global_ghost_vector(
    const MoFEMProblem *problem_ptr,const string& fiel_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1);
  PetscErrorCode set_other_global_ghost_vector(
    const string &name,const string& fiel_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1);

  //loops
  PetscErrorCode problem_basic_method_preProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_preProcess(const string &problem_name,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_postProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_postProcess(const string &problem_name,BasicMethod &method,int verb = -1);
  PetscErrorCode loop_finite_elements(const MoFEMProblem *problem_ptr,const string &fe_name,FEMethod &method,int lower_rank,int upper_rank,int verb = -1);
  PetscErrorCode loop_finite_elements(const string &problem_name,const string &fe_name,FEMethod &method,int lower_rank,int upper_rank,int verb = -1);
  PetscErrorCode loop_finite_elements(const string &problem_name,const string &fe_name,FEMethod &method,int verb = -1);
  PetscErrorCode loop_dofs(const MoFEMProblem *problem_ptr,const string &field_name,RowColData rc,EntMethod &method,int lower_rank,int upper_rank,int verb = -1);
  PetscErrorCode loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int lower_rank,int upper_rank,int verb = -1);
  PetscErrorCode loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int verb = -1);
  PetscErrorCode loop_dofs(const string &field_name,EntMethod &method,int verb = -1);

  //get multi_index form database
  PetscErrorCode get_ref_ents(const RefMoFEMEntity_multiIndex **refined_entities_ptr);
  PetscErrorCode get_ref_finite_elements(const RefMoFEMElement_multiIndex **refined_finite_elements_ptr);
  PetscErrorCode get_problem(const string &problem_name,const MoFEMProblem **problem_ptr);
  PetscErrorCode get_dofs(const DofMoFEMEntity_multiIndex **dofs_ptr);
  PetscErrorCode get_finite_elements(const MoFEMFiniteElement_multiIndex **finiteElements_ptr);

  MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_ent_moabfield_by_name_begin(const string &field_name);
  MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_ent_moabfield_by_name_end(const string &field_name);

  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_begin(const string &field_name) const;
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_end(const string &field_name) const;
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator get_dofs_by_name_and_ent_begin(const string &field_name,const EntityHandle ent);
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator get_dofs_by_name_and_ent_end(const string &field_name,const EntityHandle ent);
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator get_dofs_by_name_and_type_begin(const string &field_name,const EntityType type);
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator get_dofs_by_name_and_type_end(const string &field_name,const EntityType ent);

  EntMoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator get_fes_moabfield_by_name_begin(const string &fe_name);
  EntMoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator get_fes_moabfield_by_name_end(const string &fe_name);

  //Copy field values to another field
  PetscErrorCode field_axpy(const double alpha,const string& fiel_name_x,const string& field_name_y,bool error_if_missing = false,bool creat_if_missing = false);
  PetscErrorCode field_scale(const double alpha,const string& fiel_name);
  PetscErrorCode set_field(const double val,const EntityType type,const string& fiel_name);
  PetscErrorCode set_field(const double val,const EntityType type,const Range &ents,const string& field_name);

  //Get adjacencies
  PetscErrorCode get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities);
  PetscErrorCode get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities);
  PetscErrorCode get_adjacencies(
    const MoFEMProblem *problem_ptr,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type = Interface::INTERSECT,const int verb = 0);
  PetscErrorCode get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type = Interface::INTERSECT,const int verb = 0);

  
  //Petsc Logs
  PetscLogEvent USER_EVENT_preProcess;
  PetscLogEvent USER_EVENT_operator;
  PetscLogEvent USER_EVENT_postProcess;
  PetscLogEvent USER_EVENT_createMat;

  //
  int sIze,rAnk;

  private:

  static bool isGloballyInitialised;

};

}

#endif // __MOABFIELD_CORE_HPP__
