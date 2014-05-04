/** \file FieldCore.hpp
 * \brief Core FieldInterface class for user interface
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __MOABFIELD_CORE_HPP__
#define __MOABFIELD_CORE_HPP__

#include "FieldInterface.hpp"
#include "CoreDataStructures.hpp"

namespace MoFEM {

/** \brief Core FieldInterface class
 *
 * This class is not used directly by the user
 */
struct FieldCore: public FieldInterface {
  ErrorCode rval;
  PetscErrorCode ierr;
  //Data and low level methods 
  Tag th_Part;
  Tag th_RefType,th_RefParentHandle,th_RefBitLevel,th_RefBitEdge,th_RefFEMeshset;
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

  Interface& moab;
  int *f_shift,*MoFEMFiniteElement_shift,*p_shift;
  int verbose;

  //database

  //ref
  RefMoFEMEntity_multiIndex refinedMoFemEntities;
  RefMoFEMElement_multiIndex refinedMoFemElements;
  //field
  MoFEMField_multiIndex moabFields;
  MoFEMEntity_multiIndex entsMoabField;
  DofMoFEMEntity_multiIndex dofsMoabField;
  //finite element
  MoFEMFiniteElement_multiIndex finiteElements;
  EntMoFEMFiniteElement_multiIndex finiteElementsMoFEMEnts;
  //entFEAdjacencies
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex entFEAdjacencies;
  //moFEMProblems
  MoFEMProblem_multiIndex moFEMProblems;
  //cubit
  moabCubitMeshSet_multiIndex cubit_meshsets;

  //safty nets
  Tag th_MoFEMBuild;
  int *build_MoFEM;

  //core methods
  PetscErrorCode clear_map();
  BitFieldId get_field_shift();
  BitFEId get_BitFEId();
  BitProblemId get_problem_shift();
  PetscErrorCode initialiseDatabseInformationFromMesh(int verb = -1);
  Interface& get_moab();

  //check consistency
  PetscErrorCode check_number_of_ents_in_ents_field(const string& name);
  PetscErrorCode check_number_of_ents_in_ents_field();
  PetscErrorCode check_number_of_ents_in_ents_finite_element(const string& name);
  PetscErrorCode check_number_of_ents_in_ents_finite_element();
  PetscErrorCode rebuild_database(int verb = -1);

  //cubit meshsets
  bool check_msId_meshset(const int msId,const Cubit_BC_bitset CubitBCType);
  PetscErrorCode add_Cubit_msId(const Cubit_BC_bitset CubitBCType,const int msId);
  PetscErrorCode delete_Cubit_msId(const Cubit_BC_bitset CubitBCType,const int msId);
  PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType, const int dimension,Range &entities,const bool recursive = false);
  PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType, Range &entities,const bool recursive = false);
  PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, const int dimension,Range &entities,const bool recursive = false);
  PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, Range &entities,const bool recursive = false);
  PetscErrorCode get_Cubit_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset);
  PetscErrorCode get_Cubit_meshsets(const unsigned int CubitBCType,Range &meshsets);
  moabCubitMeshSet_multiIndex::iterator get_CubitMeshSets_begin() { return cubit_meshsets.begin(); }
  moabCubitMeshSet_multiIndex::iterator get_CubitMeshSets_end() { return cubit_meshsets.end(); }
  moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_CubitMeshSets_begin(const unsigned int CubitBCType) { 
    return cubit_meshsets.get<CubitMeshSets_mi_tag>().lower_bound(CubitBCType); 
  }
  moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_CubitMeshSets_end(const unsigned int CubitBCType) { 
    return cubit_meshsets.get<CubitMeshSets_mi_tag>().upper_bound(CubitBCType); 
  }
  moabCubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_begin(const unsigned int CubitBCType) { 
    return cubit_meshsets.get<CubitMeshSets_mask_meshset_mi_tag>().lower_bound(CubitBCType); 
  }
  moabCubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_end(const unsigned int CubitBCType) { 
    return cubit_meshsets.get<CubitMeshSets_mask_meshset_mi_tag>().upper_bound(CubitBCType); 
  }
  moabCubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator get_CubitMeshSets_byName_begin(const string& name) { 
    return cubit_meshsets.get<CubitMeshSets_name>().lower_bound(name); 
  }
  moabCubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator get_CubitMeshSets_byName_end(const string& name) { 
    return cubit_meshsets.get<CubitMeshSets_name>().upper_bound(name); 
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
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,1,msg);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode printCubitDisplacementSet() {
    PetscFunctionBegin;
    displacement_cubit_bc_data mydata;
    ierr = printCubitSet(mydata,NodeSet|mydata.type.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  PetscErrorCode printCubitPressureSet() {
    PetscFunctionBegin;
    pressure_cubit_bc_data mydata;
    ierr = printCubitSet(mydata,SideSet|mydata.type.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
   
  PetscErrorCode printCubitForceSet() {
    PetscFunctionBegin;
    force_cubit_bc_data mydata;
    ierr = printCubitSet(mydata,NodeSet|mydata.type.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode printCubitTemperatureSet() {
        PetscFunctionBegin;
        temperature_cubit_bc_data mydata;
        ierr = printCubitSet(mydata,NodeSet|mydata.type.to_ulong()); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
    
  PetscErrorCode printCubitHeatFluxSet() {
        PetscFunctionBegin;
        heatflux_cubit_bc_data mydata;
        ierr = printCubitSet(mydata,SideSet|mydata.type.to_ulong()); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

  PetscErrorCode printCubitMaterials() {
    PetscFunctionBegin;
    FieldInterface& this_mField = *this;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(this_mField,BlockSet|Mat_ElasticSet,it)) {
      Mat_Elastic data;
      ierr = it->get_attribute_data_structure(data); CHKERRQ(ierr);
      ostringstream ss;
      ss << *it << endl;
      ss << data;
      Range tets;
      rval = moab.get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
      ss << "MAT_ELATIC msId "<< it->get_msId() << " nb. tets " << tets.size() << endl;
      ss << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
      
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(this_mField,BlockSet|Mat_ThermalSet,it)) {
        Mat_Thermal data;
        ierr = it->get_attribute_data_structure(data); CHKERRQ(ierr);
        ostringstream ss;
        ss << *it << endl;
        ss << data;
        PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(this_mField,BlockSet|Mat_TransIsoSet,it)) {
        Mat_TransIso data;
        ierr = it->get_attribute_data_structure(data); CHKERRQ(ierr);
        ostringstream ss;
        ss << *it << endl;
        ss << data;
	Range tets;
	rval = moab.get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
	ss << "MAT_TRANSISO msId "<< it->get_msId() << " nb. tets " << tets.size() << endl;
	ss << endl;
        PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

  //refine
  PetscErrorCode seed_finite_elements(const Range &entities,int verb = -1);
  PetscErrorCode seed_finite_elements(const EntityHandle meshset,int verb = -1);
  PetscErrorCode seed_ref_level_2D(const Range &ents2d,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode seed_ref_level_2D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode seed_ref_level_3D(const Range &ents3d,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit);
  PetscErrorCode get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb = -1);
  PetscErrorCode get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb = -1);
  PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset);
  PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents);
  PetscErrorCode update_meshset_by_entities_children(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,
    const bool recursive = false, int verb = -1);
  PetscErrorCode update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb = -1);
  PetscErrorCode update_field_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,int verb = -1);
  PetscErrorCode update_finite_element_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb = -1);

  //remove entities
  PetscErrorCode delete_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode delete_finite_elements_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode shift_left_bit_ref(const int shif,int verb = -1);
  PetscErrorCode shift_right_bit_ref(const int shift,int verb = -1);

  //field
  PetscErrorCode add_field(
    const string& name,const BitFieldId id,const FieldSpace space,const ApproximationRank rank,enum MoFEMTypes bh = MF_EXCL,int verb = -1);
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
  PetscErrorCode dofs_NoField(const BitFieldId id,map<EntityType,int> &dof_counter);
  PetscErrorCode dofs_L2H1HcurlHdiv(const BitFieldId id,map<EntityType,int> &dof_counter,int verb = -1);
  PetscErrorCode build_fields(int verb = -1);
  PetscErrorCode clear_dofs_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_ents_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1);
  PetscErrorCode clear_dofs_fields(const string &name,const Range ents,int verb = -1);
  PetscErrorCode clear_ents_fields(const string &name,const Range enst,int verb = -1);

  //other auxiliary functions for fields
  PetscErrorCode list_dof_by_field_name(const string &name) const;
  PetscErrorCode list_ent_by_field_name(const string &name) const;
  PetscErrorCode list_dof_by_field_id(const BitFieldId id) const;
  PetscErrorCode list_ent_by_field_id(const BitFieldId id) const;
  PetscErrorCode list_field() const;
  BitFieldId get_BitFieldId(const string& name) const;
  string get_BitFieldId_name(const BitFieldId id) const;
  EntityHandle get_field_meshset(const BitFieldId id) const;
  EntityHandle get_field_meshset(const string& name) const;
  bool check_field(const string& name) const;
  const MoFEMField* get_field_structure(const string& name);

  //MoFEMFiniteElement
  PetscErrorCode add_finite_element(const string &MoFEMFiniteElement_name,enum MoFEMTypes bh = MF_EXCL);
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
  BitProblemId get_BitProblemId(const string& name) const;
  PetscErrorCode list_problem() const;

  ///add entity EntFe to finite element data databse and resolve dofs on that entity
  //loop over all finite elements, resolve its meshsets, and resolve dofs on that entitie
  PetscErrorCode build_finite_element_data_dofs(EntMoFEMFiniteElement &EntFe,int verb = -1);
  PetscErrorCode build_finite_element_uids_view(EntMoFEMFiniteElement &EntFe,int verb = -1);
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
  PetscErrorCode build_problems(int verb = -1);
  PetscErrorCode clear_problems(int verb = -1);
  PetscErrorCode simple_partition_problem(const string &name,int verb = -1);
  PetscErrorCode partition_problem(const string &name,int verb = -1);
  PetscErrorCode compose_problem(const string &name,const string &problem_for_rows,const string &problem_for_cols,int var = -1);
  PetscErrorCode compose_problem(const string &name,const string &problem_for_rows,bool copy_rows,const string &problem_for_cols,bool copy_cols,int verb = -1);
  PetscErrorCode partition_ghost_dofs(const string &name,int verb = -1);
  PetscErrorCode partition_finite_elements(const string &name,bool do_skip = true,int verb = -1);
  PetscErrorCode partition_check_matrix_fill_in(const string &problem_neme,int verb);

  //save meshsets
  PetscErrorCode problem_get_FE(const string &name,const string &fe_name,const EntityHandle meshset);

  //vector and matrices 
  PetscErrorCode VecCreateGhost(const string &name,RowColData rc,Vec *V);
  PetscErrorCode set_local_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode);
  PetscErrorCode set_global_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode);
  PetscErrorCode MatCreateMPIAIJWithArrays(const string &name,Mat *Aij,int verb = -1);
  PetscErrorCode MatCreateSeqAIJWithArrays(const string &name,Mat *Aij,PetscInt **i,PetscInt **j,PetscScalar **v,int verb = -1);
  PetscErrorCode VecScatterCreate(Vec xin,string &x_problem,RowColData x_rc,Vec yin,string &y_problem,RowColData y_rc,VecScatter *newctx,int verb = -1);

  //Mesh refine and interfaces
  PetscErrorCode get_msId_3dENTS_sides(
    const int msId,
    const Cubit_BC_bitset CubitBCType,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1);
  PetscErrorCode get_msId_3dENTS_sides(
    const EntityHandle SideSet,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1);
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const int msId,const Cubit_BC_bitset CubitBCType,
    const bool add_iterfece_entities,const bool recursive = false,int verb = -1);
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const EntityHandle SideSet,const bool add_iterfece_entities,const bool recursive = false,int verb = -1);
  PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const BitRefLevel &inheret_from_bit_level,const BitRefLevel &inheret_from_bit_level_mask,
    const EntityHandle SideSet,const bool add_iterfece_entities,const bool recursive = false,int verb = -1);

  PetscErrorCode add_prism_to_mofem_database(const EntityHandle prism,int verb = -1);

  PetscErrorCode add_verices_in_the_middel_of_edges(
    const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1);
  PetscErrorCode add_verices_in_the_middel_of_edges(const Range &edges,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface = false);
  PetscErrorCode refine_TET(const Range &test,const BitRefLevel &bit,const bool respect_interface = false);
  PetscErrorCode refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1);
  PetscErrorCode refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1);

  //loops
  PetscErrorCode problem_basic_method_preProcess(const string &problem_name,BasicMethod &method,int verb = -1);
  PetscErrorCode problem_basic_method_postProcess(const string &problem_name,BasicMethod &method,int verb = -1);
  PetscErrorCode loop_finite_elements(
    const string &problem_name,const string &fe_name,FEMethod &method,
    int lower_rank,int upper_rank,int verb = -1);
  PetscErrorCode loop_finite_elements(const string &problem_name,const string &fe_name,FEMethod &method,int verb = -1);
  PetscErrorCode loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int verb = -1);
  PetscErrorCode loop_dofs(const string &field_name,EntMethod &method,int verb = -1);

  //get multi_index form database
  PetscErrorCode get_ref_ents(const RefMoFEMEntity_multiIndex **refinedMoFemEntities_ptr);
  PetscErrorCode get_problem(const string &problem_name,const MoFEMProblem **problem_ptr);
  PetscErrorCode get_dofs(const DofMoFEMEntity_multiIndex **dofsMoabField_ptr);
  PetscErrorCode get_finite_elements(const MoFEMFiniteElement_multiIndex **finiteElements_ptr);

  MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_ent_moabfield_by_name_begin(const string &field_name);
  MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_ent_moabfield_by_name_end(const string &field_name);

  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_begin(const string &field_name) const;
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_end(const string &field_name) const;
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator get_dofs_by_name_and_ent_begin(const string &field_name,const EntityHandle ent);
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator get_dofs_by_name_and_ent_end(const string &field_name,const EntityHandle ent);
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator get_dofs_by_name_and_type_begin(const string &field_name,const EntityType type);
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator get_dofs_by_name_and_type_end(const string &field_name,const EntityType ent);

  EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_fes_moabfield_by_name_begin(const string &fe_name);
  EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_fes_moabfield_by_name_end(const string &fe_name);

  //Copy Vector of Field to Another
  PetscErrorCode set_other_global_VecCreateGhost(
    const string &name,const string& fiel_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1);
  PetscErrorCode field_axpy(const double alpha,const string& fiel_name_x,const string& field_name_y,bool error_if_missing = false,bool creat_if_missing = false);
  PetscErrorCode field_scale(const double alpha,const string& fiel_name);
  PetscErrorCode set_field(const double val,const EntityType type,const string& fiel_name);

  //Get adjacencies
  PetscErrorCode get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities);
  PetscErrorCode get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities);
  PetscErrorCode get_adjacencies(
    const MoFEMProblem *problem_ptr,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type = Interface::INTERSECT,const int verb = 0);
  PetscErrorCode get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type = Interface::INTERSECT,const int verb = 0);

  //constructor
  FieldCore(Interface& _moab,int _verbose = 1);
  ~FieldCore();


  //templates
  template<typename Tag> 
  PetscErrorCode partition_create_Mat(
    const string &name,Mat *M,const MatType type,PetscInt **_i,PetscInt **_j,PetscScalar **_v,const bool no_diagonals = true,int verb = -1);
  
  //low level finite element data
  double diffN_TET[12]; 

  //Petsc Logs
  PetscLogEvent USER_EVENT_preProcess;
  PetscLogEvent USER_EVENT_operator;
  PetscLogEvent USER_EVENT_postProcess;
};

//templates

#define PARALLEL_PARTITIONING 0
#if PARALLEL_PARTITIONING
  #define PARTITIONING_MPIADJ_COMM PETSC_COMM_WORLD
#else 
  #define PARTITIONING_MPIADJ_COMM PETSC_COMM_SELF
#endif

template<typename Tag> 
PetscErrorCode FieldCore::partition_create_Mat(
    const string &name,Mat *M,const MatType type,PetscInt **_i,PetscInt **_j,PetscScalar **_v,
    const bool no_diagonals,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    typedef typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Tag>::type NumeredDofMoFEMEntitys_by_idx;
    typedef MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type adj_by_ent;
    //find p_miit
    typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
    moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
    moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
    if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > is not found (top tip: check spelling)",name.c_str());
    //
    const NumeredDofMoFEMEntitys_by_idx &dofs_row_by_idx = p_miit->numered_dofs_rows.get<Tag>();
    const NumeredDofMoFEMEntitys_by_idx &dofs_col_by_idx = p_miit->numered_dofs_cols.get<Tag>();
    DofIdx nb_dofs_row = dofs_row_by_idx.size();
    if(p_miit->get_nb_dofs_row()!=nb_dofs_row) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    if((unsigned int)p_miit->get_nb_dofs_col()!=p_miit->numered_dofs_cols.size()) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    if(nb_dofs_row == 0) {
      SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> has zero rows",name.c_str());
    }
    typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Tag>::type::iterator miit_row,hi_miit_row;
    if(Tag::IamNotPartitioned) {
      //get range of local indices
      #if PARALLEL_PARTITIONING
      PetscLayout layout;
      ierr = PetscLayoutCreate(PETSC_COMM_WORLD,&layout); CHKERRQ(ierr);
      ierr = PetscLayoutSetBlockSize(layout,1); CHKERRQ(ierr);
      ierr = PetscLayoutSetSize(layout,nb_dofs_row); CHKERRQ(ierr);
      ierr = PetscLayoutSetUp(layout); CHKERRQ(ierr);
      PetscInt rstart,rend;
      ierr = PetscLayoutGetRange(layout,&rstart,&rend); CHKERRQ(ierr);
      ierr = PetscLayoutDestroy(&layout); CHKERRQ(ierr);
      if(verb > 0) {
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\tpartition_create_Mat: row lower %d row upper %d\n",rstart,rend);
	PetscSynchronizedFlush(PETSC_COMM_WORLD); 
      }
      miit_row = dofs_row_by_idx.lower_bound(rstart);
      hi_miit_row = dofs_row_by_idx.lower_bound(rend);
      if(distance(miit_row,hi_miit_row) != rend-rstart) {
	SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
	  "data inconsistency, distance(miit_row,hi_miit_row) != rend - rstart (%d != %d - %d = %d) ",
	  distance(miit_row,hi_miit_row),rend,rstart,rend-rstart);
      }
      #else
      miit_row = dofs_row_by_idx.begin();
      hi_miit_row = dofs_row_by_idx.end();
      #endif
    } else {
      miit_row = dofs_row_by_idx.lower_bound(pcomm->rank());
      hi_miit_row = dofs_row_by_idx.upper_bound(pcomm->rank());
    }
    int nb_loc_row_from_iterators = distance(miit_row,hi_miit_row);
    MoFEMEntity *MoFEMEntity_ptr = NULL;
    vector<PetscInt> i,j;
    vector<DofIdx> dofs_vec;
    // loop local rows
    for(;miit_row!=hi_miit_row;miit_row++) {
      i.push_back(j.size());
      if(strcmp(type,MATMPIADJ)==0) {
	DofIdx idx = Tag::get_index(miit_row);
	if(dofs_col_by_idx.find(idx)->get_unique_id()!=miit_row->get_unique_id()) {
	  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data insonsistency");
	}
      }
      if( (MoFEMEntity_ptr == NULL) ? 1 : (MoFEMEntity_ptr->get_unique_id() != miit_row->get_MoFEMEntity_ptr()->get_unique_id()) ) {
	// get field ptr
	MoFEMEntity_ptr = const_cast<MoFEMEntity*>(miit_row->get_MoFEMEntity_ptr());
	adj_by_ent::iterator adj_miit = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(MoFEMEntity_ptr->get_unique_id());
	adj_by_ent::iterator hi_adj_miit = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(MoFEMEntity_ptr->get_unique_id());
	NumeredDofMoFEMEntity_multiIndex_uid_view dofs_col_view;
	for(;adj_miit!=hi_adj_miit;adj_miit++) {
	  if(adj_miit->by_other&by_row) {
	    if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) {
	      // if element is not part of problem
	      continue; 
	    }
	    if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&miit_row->get_BitRefLevel()).none()) {
	      // if entity is not problem refinment level
	      continue; 
	    }
	    ierr = adj_miit->EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_col_dof_view( 
	      p_miit->numered_dofs_cols,dofs_col_view,Interface::UNION); CHKERRQ(ierr);
	  }
	}
	dofs_vec.resize(dofs_col_view.size());
	vector<DofIdx>::iterator vvit = dofs_vec.begin();
	NumeredDofMoFEMEntity_multiIndex_uid_view::iterator cvit;
	cvit = dofs_col_view.begin();
	for(;cvit!=dofs_col_view.end();cvit++,vvit++) {
	  int idx = Tag::get_index(*cvit);
	  if(no_diagonals) {
	    if(idx == Tag::get_index(miit_row)) {
	      continue;
	    }
	  }
	  *vvit = idx;
	}
	sort(dofs_vec.begin(),dofs_vec.end());
      }
      j.insert(j.end(),dofs_vec.begin(),dofs_vec.end());
    }
    //build adj matrix
    i.push_back(j.size());
    ierr = PetscMalloc(i.size()*sizeof(PetscInt),_i); CHKERRQ(ierr);
    ierr = PetscMalloc(j.size()*sizeof(PetscInt),_j); CHKERRQ(ierr);
    copy(i.begin(),i.end(),*_i);
    copy(j.begin(),j.end(),*_j);
    PetscInt nb_row_dofs = p_miit->get_nb_dofs_row();
    PetscInt nb_col_dofs = p_miit->get_nb_dofs_col();
    if(strcmp(type,MATMPIADJ)==0) { 
      if(i.size()-1 != (unsigned int)nb_loc_row_from_iterators) {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
      }
      ierr = MatCreateMPIAdj(PARTITIONING_MPIADJ_COMM,i.size()-1,nb_col_dofs,*_i,*_j,PETSC_NULL,M); CHKERRQ(ierr);
      ierr = MatSetOption(*M,MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);
    } else if(strcmp(type,MATMPIAIJ)==0) {
      if(i.size()-1 != (unsigned int)nb_loc_row_from_iterators) {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
      }
      PetscInt nb_local_dofs_row = p_miit->get_nb_local_dofs_row();
      if((unsigned int)nb_local_dofs_row!=i.size()-1) {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
      }
      PetscInt nb_local_dofs_col = p_miit->get_nb_local_dofs_col();
      ierr = ::MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD,nb_local_dofs_row,nb_local_dofs_col,nb_row_dofs,nb_col_dofs,*_i,*_j,PETSC_NULL,M); CHKERRQ(ierr);
    } else if(strcmp(type,MATAIJ)==0) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"not implemented");
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"not implemented");
    }
    PetscFunctionReturn(0);
  }

}

#endif // __MOABFIELD_CORE_HPP__
