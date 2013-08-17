/** \file moabField_Core.cpp
 * \brief Myltindex containes, data structures and other low-level functions 
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
 *
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

#include<moabField_Core.hpp>
#include<FEM.h>

namespace MoFEM {

const static int debug = 1;

moabField_Core::moabField_Core(Interface& _moab,int _verbose): 
  moab(_moab),verbose(_verbose) {
  const EntityHandle root_meshset = moab.get_root_set();
  // Version
  stringstream strs_version;
  strs_version << "MoFEM_version_" << MoFEM_VERSION_MAJOR << "." << MoFEM_VERSION_MINOR;
  Tag th_version;
  string version = strs_version.str();
  rval = moab.tag_get_handle("_MoFEM_VERSION",version.size()*sizeof(char),MB_TYPE_OPAQUE,
    th_version,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,NULL); CHKERR(rval);
  const char *ptr_version = version.c_str();
  rval = moab.tag_set_data(th_version,&root_meshset,1,ptr_version); CHKERR_THROW(rval);
  //tags saved in vtk-files
  const int def_part = -1;
  rval = moab.tag_get_handle("PARTITION",1,MB_TYPE_INTEGER,th_Part,MB_TAG_CREAT|MB_TAG_SPARSE,&def_part); CHKERR(rval);
  //Tags Ref
  rval = moab.tag_get_handle("_RefParentHandle",1,MB_TYPE_HANDLE,th_RefParentHandle,MB_TAG_CREAT|MB_TAG_SPARSE,&root_meshset); CHKERR(rval); CHKERR_THROW(rval);
  const int def_type[] = {0,0};
  rval = moab.tag_get_handle("_RefType",2,MB_TYPE_INTEGER,th_RefType,MB_TAG_CREAT|MB_TAG_SPARSE,def_type); CHKERR(rval);
  BitRefLevel def_bit_level = 0;
  rval = moab.tag_get_handle("_RefBitLevel",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_RefBitLevel,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); CHKERR(rval);
  BitRefEdges def_bit_egde = 0;
  rval = moab.tag_get_handle("_RefBitEdge",sizeof(BitRefEdges),MB_TYPE_OPAQUE,
    th_RefBitEdge,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_egde); CHKERR(rval);
  //Tags Field
  const unsigned long int def_id = 0;
  rval = moab.tag_get_handle("_FieldId",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FieldId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR(rval);
  FieldSpace def_space = LastSpace;
  rval = moab.tag_get_handle("_FieldSpace",sizeof(FieldSpace),MB_TYPE_OPAQUE,
    th_FieldSpace,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_space); CHKERR(rval);
  const int def_val_len = 0;
  rval = moab.tag_get_handle("_FieldName",def_val_len,MB_TYPE_OPAQUE,
    th_FieldName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR(rval);
  //Tags FE
  rval = moab.tag_get_handle("_FEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_FEId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR(rval);
  rval = moab.tag_get_handle("_FEName",def_val_len,MB_TYPE_OPAQUE,
    th_FEName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR(rval);
  rval = moab.tag_get_handle("_FEIdCol",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR(rval);
  rval = moab.tag_get_handle("_FEIdRow",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR(rval);
  rval = moab.tag_get_handle("_FEIdData",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR(rval);
  //Tags Problem
  rval = moab.tag_get_handle("_ProblemId",sizeof(BitProblemId),MB_TYPE_OPAQUE,
    th_ProblemId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR(rval);
  rval = moab.tag_get_handle("_ProblemFEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_ProblemFEId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR(rval);
  rval = moab.tag_get_handle("_ProblemName",def_val_len,MB_TYPE_OPAQUE,
    th_ProblemName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR(rval);
  DofIdx def_nbdofs = 0;
  rval = moab.tag_get_handle("_ProblemNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR(rval);
  rval = moab.tag_get_handle("_ProblemNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR(rval);
  rval = moab.tag_get_handle("_ProblemLocalNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR(rval);
  rval = moab.tag_get_handle("_ProblemGhostNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR(rval);
  rval = moab.tag_get_handle("_ProblemLocalNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR(rval);
  rval = moab.tag_get_handle("_ProblemGhostNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR(rval);
  //Global Variables
  //Fields
  int def_shift = 1;
  rval = moab.tag_get_handle("_FieldShift",1,MB_TYPE_INTEGER,th_FieldShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR(rval);
  const void* tag_data[1];
  rval = moab.tag_get_by_ptr(th_FieldShift,&root_meshset,1,tag_data); CHKERR(rval);
  f_shift = (int*)tag_data[0];
  //FE
  rval = moab.tag_get_handle("_FEShift",1,MB_TYPE_INTEGER,th_FEShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FEShift,&root_meshset,1,tag_data); CHKERR(rval);
  MoFEMFiniteElement_shift = (int*)tag_data[0];
  //Problem
  rval = moab.tag_get_handle("_ProblemShift",1,MB_TYPE_INTEGER,th_ProblemShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemShift,&root_meshset,1,tag_data); CHKERR(rval);
  p_shift = (int*)tag_data[0];
  //SaftyNets
  int def_bool = 0;
  rval = moab.tag_get_handle("_MoFEMBuild",1,MB_TYPE_INTEGER,th_MoFEMBuild,MB_TAG_CREAT|MB_TAG_MESH,&def_bool); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  rval = moab.tag_get_by_ptr(th_MoFEMBuild,&root_meshset,1,(const void **)&build_MoFEM); CHKERR(rval);
  //Meshsets
  int default_val = -1;
  rval = moab.tag_get_handle(DIRICHLET_SET_TAG_NAME,1, MB_TYPE_INTEGER, 
    nsTag, MB_TAG_SPARSE|MB_TAG_CREAT, &default_val); CHKERR(rval);
  rval = moab.tag_get_handle(NEUMANN_SET_TAG_NAME,1, MB_TYPE_INTEGER, 
    ssTag, MB_TAG_SPARSE|MB_TAG_CREAT, &default_val); CHKERR(rval);
  const int def_bc_data_len = 0;
  std::string tag_name = std::string(DIRICHLET_SET_TAG_NAME)+"__BC_DATA";
  rval = moab.tag_get_handle(tag_name.c_str(),def_bc_data_len,MB_TYPE_OPAQUE,
    nsTag_data,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR(rval);
  tag_name = std::string(NEUMANN_SET_TAG_NAME)+"__BC_DATA";
  rval = moab.tag_get_handle(tag_name.c_str(),def_bc_data_len,MB_TYPE_OPAQUE,
    ssTag_data,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR(rval);
  rval = moab.tag_get_handle(MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER,
    bhTag,MB_TAG_SPARSE|MB_TAG_CREAT,&default_val); CHKERR(rval);
  std::vector<unsigned int> def_uint_zero(12,0);
  rval= moab.tag_get_handle("BLOCK_HEADER",12*sizeof(unsigned int),MB_TYPE_INTEGER,
    bhTag_header,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_uint_zero[0]); CHKERR(rval); 
  // For VTK files
  int def_elem_type = MBMAXTYPE;
  rval = moab.tag_get_handle("ElemType",1,MB_TYPE_INTEGER,th_ElemType,MB_TAG_CREAT|MB_TAG_SPARSE,&def_elem_type); CHKERR(rval); 
  //
  map_from_mesh(3); 
  //
  ShapeDiffMBTET(diffN_TET); 
  // Petsc Logs
  PetscLogEventRegister("FE_preProcess",0,&USER_EVENT_preProcess);
  PetscLogEventRegister("FE_operator",0,&USER_EVENT_operator);
  PetscLogEventRegister("FE_postProcess",0,&USER_EVENT_postProcess);

}
moabField_Core::~moabField_Core() {
}
Interface& moabField_Core::get_moab() {
  return moab;
}
BitFieldId moabField_Core::get_BitFieldId(const string& name) const {
  typedef MoFEMField_multiIndex::index<FieldName_mi_tag>::type field_set_by_name;
  const field_set_by_name &set = moabfields.get<FieldName_mi_tag>();
  field_set_by_name::iterator miit = set.find(name);
  if(miit==set.end()) THROW_AT_LINE("field not in databse (top tip: check spelling)");
  return miit->get_id();
}
string moabField_Core::get_BitFieldId_name(const BitFieldId id) const {
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set = moabfields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set.find(id);
  return miit->get_name();
}
EntityHandle moabField_Core::get_field_meshset(const BitFieldId id) const {
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set = moabfields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set.find(id);
  if(miit==set.end()) THROW_AT_LINE("field not in databse (top tip: check spelling)");
  return miit->meshset;
}
EntityHandle moabField_Core::get_field_meshset(const string& name) const {
  return get_field_meshset(get_BitFieldId(name));
}
BitFieldId moabField_Core::get_field_shift() {
  assert((unsigned int)*f_shift<BitFieldId().set().to_ulong());
  return (BitFieldId)(1<<(((*f_shift)++)-1)); 
}
BitFEId moabField_Core::get_BitFEId() {
  assert((unsigned int)*MoFEMFiniteElement_shift<BitFEId().set().to_ulong());
  return BitFEId(1<<(((*MoFEMFiniteElement_shift)++)-1)); 
}
BitProblemId moabField_Core::get_problem_shift() {
  assert((unsigned int)*p_shift<BitProblemId().set().to_ulong());
  return BitProblemId(1<<(((*p_shift)++)-1)); 
}
PetscErrorCode moabField_Core::clear_map() {
  PetscFunctionBegin;
  cubit_meshsets.clear();
  refined_mofem_entities.clear();
  refined_mofem_elements.clear();
  moabfields.clear();
  ents_moabfield.clear();
  dofs_moabfield.clear();
  finite_elements.clear();
  finite_elements_moabents.clear();
  adjacencies.clear();
  problems.clear();
  PetscFunctionReturn(0);
} 
PetscErrorCode moabField_Core::add_field(const string& name,const BitFieldId id,const FieldSpace space,const ApproximationRank rank,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
  //id
  rval = moab.tag_set_data(th_FieldId,&meshset,1,&id); CHKERR_PETSC(rval);
  //space
  rval = moab.tag_set_data(th_FieldSpace,&meshset,1,&space); CHKERR_PETSC(rval);
  //add meshset to ref_ents // meshset dof on all level sets
  if(space == NoField) {
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,meshset));
    bool success = refined_mofem_entities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(BitRefLevel().set()));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  }
  //id name
  void const* tag_data[] = { name.c_str() };
  int tag_sizes[] = { name.size() };
  rval = moab.tag_set_by_ptr(th_FieldName,&meshset,1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  Tag th_AppOrder,th_FieldData,th_Rank,th_AppDofOrder,th_DofRank;
  //data
  string Tag_data_name = "_App_Data_"+name;
  const int def_len = 0;
  rval = moab.tag_get_handle(Tag_data_name.c_str(),def_len,MB_TYPE_OPAQUE,
    th_FieldData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
  //order
  ApproximationOrder def_ApproximationOrder = -1;
  string Tag_ApproximationOrder_name = "_App_Order_"+name;
  rval = moab.tag_get_handle(Tag_ApproximationOrder_name.c_str(),sizeof(ApproximationOrder),MB_TYPE_OPAQUE,
    th_AppOrder,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_ApproximationOrder); CHKERR_PETSC(rval);
  //dof order
  string Tag_dof_ApproximationOrder_name = "_App_Dof_Order"+name;
  rval = moab.tag_get_handle(Tag_dof_ApproximationOrder_name.c_str(),def_len,MB_TYPE_OPAQUE,
    th_AppDofOrder,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
  //rank
  int def_rank = 1;
  string Tag_rank_name = "_Field_Rank_"+name;
  rval = moab.tag_get_handle(Tag_rank_name.c_str(),sizeof(ApproximationRank),MB_TYPE_OPAQUE,
    th_Rank,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_rank); CHKERR_PETSC(rval);
  rval = moab.tag_set_data(th_Rank,&meshset,1,&rank); CHKERR_PETSC(rval);
  //dof rank
  string Tag_dof_rank_name = "_Field_Dof_Rank_"+name;
  rval = moab.tag_get_handle(Tag_dof_rank_name.c_str(),def_len,MB_TYPE_OPAQUE,
    th_DofRank,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
  //add meshset
  pair<MoFEMField_multiIndex::iterator,bool> p;
  try {
    p = moabfields.insert(MoFEMField(moab,meshset));  
    if(!p.second) SETERRQ(PETSC_COMM_SELF,1,"field not inesrted");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  if(verbose > 0) {
    ostringstream ss;
    ss << "add: " << *p.first << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  //
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_field(const string& name,const FieldSpace space,const ApproximationRank rank,int verb) {
  PetscFunctionBegin;
  *build_MoFEM = 0;
  if(verb==-1) verb = verbose;
  BitFieldId id = get_field_shift();
  ierr = add_field(name,id,space,rank,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::map_from_mesh(int verb) {
  PetscFunctionBegin;
  //ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(verb==-1) verb = verbose;
  ierr = clear_map(); CHKERRQ(ierr);
  Range meshsets;
  rval = moab.get_entities_by_type(0,MBENTITYSET,meshsets,false);  CHKERR_PETSC(rval);
  Range::iterator mit = meshsets.begin();
  for(;mit!=meshsets.end();mit++) {
    CubitMeshSets base_meshset(moab,*mit);
    if((base_meshset.CubitBCType&Cubit_BC_bitset(NodeSet|SideSet|BlockSet)).any()) {
      pair<moabBaseMeshSet_multiIndex::iterator,bool> p = cubit_meshsets.insert(base_meshset);
      if(!p.second) SETERRQ(PETSC_COMM_SELF,1,"meshset not inserted");
      ostringstream ss;
      if(verb > 0) {
	ss << "read cubit" << base_meshset << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      //PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      //PetscSynchronizedFlush(PETSC_COMM_WORLD); 
      ierr = seed_ref_level_MESHSET(*mit,0); CHKERRQ(ierr);
    }
    BitFieldId field_id;
    rval = moab.tag_get_data(th_FieldId,&*mit,1,&field_id); CHKERR_PETSC(rval);
    if(field_id!=0) {
      pair<MoFEMField_multiIndex::iterator,bool> p;
      try {
	p = moabfields.insert(MoFEMField(moab,*mit));
	if(verb > 0) {
	  ostringstream ss;
	  ss << "read field " << *p.first << endl;;
	  PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	}
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } 
      if(p.first->get_space()==NoField) {
	assert(p.first->meshset == *mit);
	//add field to ref ents
	pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*mit));
	NOT_USED(p_ref_ent);
      } else {
	Range ents;
	rval = moab.get_entities_by_handle(*mit,ents,false); CHKERR_PETSC(rval);
	Range::iterator eit = ents.begin();
	for(;eit!=ents.end();eit++) {
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*eit));
	  try {
	    MoFEMEntity moabent(moab,&*p.first,&*p_ref_ent.first);
	    if(moabent.get_order_nb_dofs(moabent.get_max_order())==0) continue; 
	    pair<MoFEMEntity_multiIndex::iterator,bool> p_ent = ents_moabfield.insert(moabent);
	    NOT_USED(p_ent);
	  } catch (const std::exception& ex) {
	    ostringstream ss;
	    ss << ex.what() << endl;
	    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	  }
	}
      }
    }
    BitFieldId fe_id;
    rval = moab.tag_get_data(th_FEId,&*mit,1,&fe_id); CHKERR_PETSC(rval);
    if(fe_id!=0) {
      pair<MoFEMFiniteElement_multiIndex::iterator,bool> p = finite_elements.insert(MoFEMFiniteElement(moab,*mit));
      if(verb > 0) {
     	ostringstream ss;
	ss << "read finite element " << *p.first << endl;;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      NOT_USED(p);
      assert(p.first->meshset == *mit);
      Range ents;
      rval = moab.get_entities_by_type(*mit,MBENTITYSET,ents,false); CHKERR_PETSC(rval);
      rval = moab.get_entities_by_handle(*mit,ents,true); CHKERR_PETSC(rval);
      Range::iterator eit = ents.begin();
      for(;eit!=ents.end();eit++) {
	pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*eit));
	pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
	switch (moab.type_from_handle(*eit)) {
	  case MBTET:
	    p_MoFEMFiniteElement = refined_mofem_elements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_ref_ent.first)));
	    assert(p_MoFEMFiniteElement.first->get_BitRefEdges_ulong()!=-1);
	    break;
	  case MBPRISM:
  	    p_MoFEMFiniteElement = refined_mofem_elements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ref_ent.first)));
	    break;
	  case MBENTITYSET:
  	    p_MoFEMFiniteElement = refined_mofem_elements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_MESHSET(moab,&*p_ref_ent.first)));
	    break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"Only finite elements of type MBTET, MBPRISM and MBENTITYSET are implemented");
	}
      }
    }
    BitProblemId problem_id;
    rval = moab.tag_get_data(th_ProblemId,&*mit,1,&problem_id); CHKERR_PETSC(rval);
    if(problem_id!=0) {
      pair<MoFEMProblem_multiIndex::iterator,bool> p = problems.insert(_MoFEMProblem_(moab,*mit));
      if(verb > 0) {
     	ostringstream ss;
	ss << "read problem " << *p.first << endl;;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    }
  }
  //add prisms, this is needed for refinment
  Range prisms;
  rval = moab.get_entities_by_type(0,MBPRISM,prisms,true);  CHKERR_PETSC(rval);
  Range::iterator pit = prisms.begin();
  for(;pit!=prisms.end();pit++) {
    ierr = add_prism_to_adjacencies_maps_for_prisms(*pit); CHKERRQ(ierr);
  }
  if(verb > 2) {
    list_field();
    list_finite_elements();
    list_problem();
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_ents_to_field_by_TRIs(const EntityHandle meshset,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  Range nodes,tris,edges;
  rval = moab.get_entities_by_type(meshset,MBTRI,tris,false); CHKERR_PETSC(rval);
  switch (space) {
    case L2:
      rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
      if(verb>1) {
	ostringstream ss;
	ss << "add entities to field " << get_BitFieldId_name(id);
	ss << " nb. add tris " << tris.size();
	ss << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      break;
    case H1:
      rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
      rval = moab.get_connectivity(tris,nodes,true); CHKERR_PETSC(rval);
      rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(tris,1,false,edges,Interface::UNION); CHKERR_PETSC(rval);
      rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
      if(verb>1) {
	ostringstream ss;
	ss << "add entities to field " << get_BitFieldId_name(id);
	ss << " nb. add tris " << tris.size();
	ss << " nb. add edges " << edges.size();
	ss << " nb. add nodes " << nodes.size();
	ss << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      break;
    case Hcurl:
      SETERRQ(PETSC_COMM_SELF,1,"sorry, not implemented, Hcurl not implented for TRI");
      break;
    case Hdiv:
      SETERRQ(PETSC_COMM_SELF,1,"sorry, not implemented, Hdiv not implemented for TRI");
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"add_ents_to_field_by_TRIs this field not work for TRIs");
  }
  ierr = seed_ref_level_3D(idm,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_ents_to_field_by_TRIs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TRIs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_ents_to_field_by_TETs(const EntityHandle meshset,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  Range tets,nodes,tris,edges;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,false); CHKERR_PETSC(rval);
  switch (space) {
    case L2:
      rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
      if(verb>1) {
	ostringstream ss;
	ss << "add entities to field " << get_BitFieldId_name(id);
	ss << " nb. add tets " << tets.size();
	ss << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      break;
    case H1:
      rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
      rval = moab.get_connectivity(tets,nodes,true); CHKERR_PETSC(rval);
      rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(tets,2,false,tris,Interface::UNION); CHKERR_PETSC(rval);
      rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(tets,1,false,edges,Interface::UNION); CHKERR_PETSC(rval);
      rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
      if(verb>1) {
	ostringstream ss;
	ss << "add entities to field " << get_BitFieldId_name(id);
	ss << " nb. add tets " << tets.size();
	ss << " nb. add tris " << tris.size();
	ss << " nb. add edges " << edges.size();
	ss << " nb. add nodes " << nodes.size();
	ss << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      break;
    case Hcurl:
      rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(tets,2,false,tris,Interface::UNION); CHKERR_PETSC(rval);
      rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(tets,1,false,edges,Interface::UNION); CHKERR_PETSC(rval);
      rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
      if(verb>1) {
	ostringstream ss;
	ss << "add entities to field " << get_BitFieldId_name(id);
	ss << " nb. add tets " << tets.size();
	ss << " nb. add tris " << tris.size();
	ss << " nb. add edges " << edges.size();
	ss << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      break;
    case Hdiv:
      rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(tets,2,false,tris,Interface::UNION); CHKERR_PETSC(rval);
      rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
      if(verb>1) {
	ostringstream ss;
	ss << "add entities to field " << get_BitFieldId_name(id);
	ss << " nb. add tets " << tets.size();
	ss << " nb. add tris " << tris.size();
	ss << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"add_ents_to_field_by_TETs this field not work for TETs");
  }
  ierr = seed_ref_level_3D(idm,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TETs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::set_field_order(const EntityHandle meshset,const EntityType type,const BitFieldId id,const ApproximationOrder order) {
  PetscFunctionBegin;
  *build_MoFEM = 0;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents); CHKERR_PETSC(rval);
  //check field & meshset
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabfields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit==set_id.end()) SETERRQ(PETSC_COMM_SELF,1,"no id found"); 
  EntityHandle idm = no_handle;
  try {
   idm = get_field_meshset(id);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  //itersection with field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(idm,ents_of_id_meshset,false); CHKERR_PETSC(rval);
  ents = intersect(ents,ents_of_id_meshset);
  vector<const void*> tag_data_order(ents.size());
  rval = moab.tag_get_by_ptr(miit->th_AppOrder,ents,&tag_data_order[0]); CHKERR_PETSC(rval);
  //ent view by field id (in set all MoabEnts has the same FieldId)
  typedef MoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type ent_set_by_id;
  ent_set_by_id& set = ents_moabfield.get<BitFieldId_mi_tag>();
  ent_set_by_id::iterator miit2 = set.lower_bound(id);
  MoFEMEntity_multiIndex_ent_view ents_id_view;
  if(miit2 != set.end()) {
    ent_set_by_id::iterator hi_miit2 = set.upper_bound(id);
    for(;miit2!=hi_miit2;miit2++) {
      ents_id_view.insert(&*miit2);
    }
  }
  //loop over ents
  Range::iterator eit = ents.begin();
  for(unsigned int ee = 0;ee<ents.size();ee++,eit++) {
    MoFEMEntity_multiIndex_ent_view::iterator miit3 = ents_id_view.find(*eit);
    if(miit3!=ents_id_view.end()) {
      const ApproximationOrder old_ApproximationOrder = (*miit3)->get_max_order();
      if(old_ApproximationOrder==order) continue;
      MoFEMEntity_multiIndex::iterator miit4 = ents_moabfield.get<Unique_mi_tag>().find((*miit3)->get_unique_id());
      assert(miit4!=ents_moabfield.end());
      typedef DofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type dof_set_type;
      dof_set_type& set_set = dofs_moabfield.get<Composite_mi_tag2>();
      dof_set_type::iterator miit5 = set_set.lower_bound(boost::make_tuple(miit4->get_name(),miit4->get_ent()));
      dof_set_type::iterator hi_miit6 = set_set.upper_bound(boost::make_tuple(miit4->get_name(),miit4->get_ent()));
      for(;miit5!=hi_miit6;miit5++) {
	if(miit5->get_dof_order()<=order) continue;
	bool success = dofs_moabfield.modify(dofs_moabfield.project<0>(miit5),DofMoFEMEntity_active_change(false));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }
    } else {
      *(ApproximationOrder*)tag_data_order[ee] = order;
      RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refined_mofem_entities.get<MoABEnt_mi_tag>().find(*eit);
      if(miit_ref_ent==refined_mofem_entities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
      try { 
	MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
	if(moabent.get_order_nb_dofs(moabent.get_max_order())==0) continue; 
	pair<MoFEMEntity_multiIndex::iterator,bool> e_miit = ents_moabfield.insert(moabent);
	bool success = ents_moabfield.modify(e_miit.first,MoFEMEntity_change_order(moab,order));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
    }
  }
  if(verbose>2) {
    list_ent_by_id(id);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::set_field_order(const EntityHandle meshset,const EntityType type,const string& name,const ApproximationOrder order) {
  PetscFunctionBegin;
  *build_MoFEM = 0;
  try{
    ierr = set_field_order(meshset,type,get_BitFieldId(name),order); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::dofs_NoField(const BitFieldId id,int &dof_counter) {
  PetscFunctionBegin;
  //field it
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabfields.get<BitFieldId_mi_tag>();
  //find fiels
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) SETERRQ(PETSC_COMM_SELF,1,"field no found");
  //serch if field meshset is in database
  RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refined_mofem_entities.get<MoABEnt_mi_tag>().find(miit->meshset);
  if(miit_ref_ent==refined_mofem_entities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
  pair<MoFEMEntity_multiIndex::iterator,bool> e_miit;
  try {
    //create database entity
    MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
    e_miit = ents_moabfield.insert(moabent);
    //this is nor real field in space (set order to zero)
    bool success = ents_moabfield.modify(e_miit.first,MoFEMEntity_change_order(moab,0));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  assert(e_miit.first->get_ent()==miit->meshset);
  ApproximationRank rank = 0;
  //create dofs on this entity (nb. of dofs is equal to rank)
  for(;rank<e_miit.first->get_max_rank();rank++) {
    pair<DofMoFEMEntity_multiIndex::iterator,bool> d_miit;
    //check if dof is in darabase
    d_miit.first = dofs_moabfield.project<0>(
      dofs_moabfield.get<Unique_mi_tag>().find(DofMoFEMEntity::get_unique_id_calculate(rank,&*(e_miit.first)))
    );
    //if dof is not in databse
    if(d_miit.first==dofs_moabfield.end()) {
      //insert dof
      d_miit = dofs_moabfield.insert(DofMoFEMEntity(&*(e_miit.first),0,rank,rank));
      if(d_miit.second) dof_counter++;
      bool success = dofs_moabfield.modify(d_miit.first,DofMoFEMEntity_active_change(true));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    //check consistency
    assert(d_miit.first->get_ent()==e_miit.first->get_MoFEMField_ptr()->meshset);
    assert(d_miit.first->get_ent_type()==e_miit.first->get_ent_type());
    assert(d_miit.first->get_id()==e_miit.first->get_id());
  }
  if(verbose>2) {
    typedef DofMoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type dof_set_by_id;
    dof_set_by_id &set = dofs_moabfield.get<BitFieldId_mi_tag>();
    dof_set_by_id::iterator miit2 = set.lower_bound(id);
    dof_set_by_id::iterator hi_miit2 = set.upper_bound(id);
    assert(miit2!=hi_miit2);
    for(;miit2!=hi_miit2;miit2++) {
      ostringstream ss;
      ss << *miit2 << endl;;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::dofs_L2H1HcurlHdiv(const BitFieldId id,int &dof_counter,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ents_by_ents;
  //find field
  const field_set_by_id &set_id = moabfields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) SETERRQ(PETSC_COMM_SELF,1,"field no found");
  //ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(miit->meshset,ents_of_id_meshset,false); CHKERR_PETSC(rval);
  //create dofs_moabfield
  Range::iterator eit = ents_of_id_meshset.begin();
  for(;eit!=ents_of_id_meshset.end();eit++) {
    // check if ent is in ref meshset
    RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refined_mofem_entities.get<MoABEnt_mi_tag>().find(*eit);
    if(miit_ref_ent==refined_mofem_entities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
    //pair<MoFEMEntity_multiIndex::iterator,bool> e_miit;
    MoFEMEntity_multiIndex::iterator e_miit;
    try {
      e_miit = ents_moabfield.find(MoFEMEntity(moab,&*miit,&*miit_ref_ent).get_unique_id());
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,1,msg);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss <<  ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    // create mofem entity linked to ref ent
    try {
      e_miit = ents_moabfield.find(MoFEMEntity(moab,&*miit,&*miit_ref_ent).get_unique_id());
    } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    if(e_miit == ents_moabfield.end()) {
      ApproximationOrder order = -1;
      rval = moab.tag_set_data(miit->th_AppOrder,&*eit,1,&order); CHKERR_PETSC(rval);
      pair<MoFEMEntity_multiIndex::iterator,bool> p_e_miit;
      try {
	MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
	p_e_miit = ents_moabfield.insert(moabent);
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      if(!p_e_miit.second) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      bool success = ents_moabfield.modify(p_e_miit.first,MoFEMEntity_change_order(moab,-1));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      e_miit = p_e_miit.first;
    }
    // insert dofmoabent into mofem databse
    int DD = 0;
    int oo = 0;
    for(;oo<=e_miit->get_max_order();oo++) {
      for(int dd = 0;dd<e_miit->get_order_nb_dofs_diff(oo);dd++) {
	for(int rr = 0;rr<e_miit->get_max_rank();rr++,DD++) {
	  pair<DofMoFEMEntity_multiIndex::iterator,bool> d_miit;
	  try {
	    DofMoFEMEntity mdof(&*(e_miit),oo,rr,DD);
	    d_miit = dofs_moabfield.insert(mdof);
	    if(d_miit.second) {
	      dof_counter++;
	      bool success = dofs_moabfield.modify(d_miit.first,DofMoFEMEntity_active_change(true));
	      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
	    } 
	    //check ent
	    assert(d_miit.first->get_ent()==e_miit->get_ent());
	    assert(d_miit.first->get_ent_type()==e_miit->get_ent_type());
	    assert(d_miit.first->get_id()==e_miit->get_id());
	    //check dof
	    if(d_miit.first->get_dof_order()!=oo) {
	      ostringstream ss;
	      ss << "data inconsistency!" << endl;
	      ss << "shuld be " << mdof << endl;
	      ss << "but is " << *d_miit.first << endl;
	      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	    }
	    assert(d_miit.first->get_max_order()==e_miit->get_max_order());
	  } catch (const char* msg) {
	    SETERRQ(PETSC_COMM_SELF,1,msg);
	  }
	}
      }
    }
    assert(DD == e_miit->get_max_rank()*e_miit->get_order_nb_dofs(e_miit->get_max_order()));
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::build_fields(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  field_set_by_id &set_id = moabfields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    int dof_counter = 0;
    if(verbose>0) {
      PetscPrintf(PETSC_COMM_WORLD,"Build Field %s ",miit->get_name().c_str());
    }
    switch (miit->get_space()) {
      case NoField:
	ierr = dofs_NoField(miit->get_id(),dof_counter); CHKERRQ(ierr);
	break;
      case L2:
      case H1:
      case Hcurl:
      case Hdiv:
	ierr = dofs_L2H1HcurlHdiv(miit->get_id(),dof_counter,verb); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    if(verbose>0) {
      PetscPrintf(PETSC_COMM_WORLD,"nb added dofs %d\n",dof_counter);
    }
    if(verb>1) {
      list_ent_by_id(miit->get_id());
      list_dof_by_id(miit->get_id());
    }
  }
  PetscPrintf(PETSC_COMM_WORLD,"Nb. dofs %u\n",dofs_moabfield.size());
  *build_MoFEM = 1<<0;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::list_dof_by_id(const BitFieldId id) const {
  PetscFunctionBegin;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank %d Field dofs %s\n",pcomm->rank(),get_BitFieldId_name(id).c_str()); //!!! Syn ...
  typedef DofMoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type dof_set_by_id;
  const dof_set_by_id &set_id = dofs_moabfield.get<BitFieldId_mi_tag>();
  dof_set_by_id::iterator miit = set_id.lower_bound(id);
  dof_set_by_id::iterator hi_miit = set_id.upper_bound(id);
  DofMoFEMEntity_multiIndex_order_view dofs_order_view;
  for(;miit!=hi_miit;miit++) {
    dofs_order_view.insert(&*miit);
  }
  DofMoFEMEntity_multiIndex_order_view::iterator miit2 =  dofs_order_view.lower_bound(0);
  DofMoFEMEntity_multiIndex_order_view::iterator hi_miit2 =  dofs_order_view.upper_bound(max_ApproximationOrder);
  DofMoFEMEntity_multiIndex_ent_type_view dofs_type_view;
  for(;miit2!=hi_miit2;miit2++) {
    test_moab(moab,(*miit2)->get_ent());
    dofs_type_view.insert(*miit2);
  }
  DofMoFEMEntity_multiIndex_ent_type_view::iterator miit3 =  dofs_order_view.lower_bound(MBVERTEX);
  DofMoFEMEntity_multiIndex_ent_type_view::iterator hi_miit3 =  dofs_order_view.upper_bound(MBTET);
  ostringstream ss;
  for(;miit3!=hi_miit3;miit3++) {
    ss << "rank " << pcomm->rank() << " ";
    ss << **miit3 << endl;
  }
  //PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::list_ent_by_id(const BitFieldId id) const {
  PetscFunctionBegin;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank %d Field ents %s\n",pcomm->rank(),get_BitFieldId_name(id).c_str()); //!!!! Syn ...
  typedef MoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type ent_set_by_id;
  const ent_set_by_id &set = ents_moabfield.get<BitFieldId_mi_tag>();
  ent_set_by_id::iterator miit = set.lower_bound(id);
  ent_set_by_id::iterator hi_miit = set.upper_bound(id);
  ostringstream ss;
  for(;miit!=hi_miit;miit++) {
    ss << "rank " << pcomm->rank() << " ";
    ss << *miit << endl;
  }
  //PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::list_field() const {
  PetscFunctionBegin;
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabfields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_finite_element(const string &MoFEMFiniteElement_name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finite_elements_by_name;
  finite_elements_by_name &MoFEMFiniteElement_name_set = finite_elements.get<MoFEMFiniteElement_name_mi_tag>();
  finite_elements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement!=MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,1,"this MoFEMFiniteElement is there");
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
  //id
  BitFEId id = get_BitFEId();
  rval = moab.tag_set_data(th_FEId,&meshset,1,&id); CHKERR_PETSC(rval);
  //id name
  void const* tag_data[] = { MoFEMFiniteElement_name.c_str() };
  int tag_sizes[] = { MoFEMFiniteElement_name.size() };
  rval = moab.tag_set_by_ptr(th_FEName,&meshset,1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  //tags
  Tag th_FEMatData,th_FEVecData;
  string Tag_mat_name = "_FE_MatData_"+MoFEMFiniteElement_name;
  const int def_val_len = 0;
  rval = moab.tag_get_handle(Tag_mat_name.c_str(),def_val_len,MB_TYPE_OPAQUE,
    th_FEMatData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
  string Tag_vec_name = "_FE_VecData_"+MoFEMFiniteElement_name;
  rval = moab.tag_get_handle(Tag_vec_name.c_str(),def_val_len,MB_TYPE_OPAQUE,
    th_FEVecData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
  //tags uids
  Tag th_DofUidRow,th_DofUidCol,th_DofUidData;
  string Tag_DofUidRow_name = "_DofUidRow_"+MoFEMFiniteElement_name;
  rval = moab.tag_get_handle(Tag_DofUidRow_name.c_str(),def_val_len,MB_TYPE_OPAQUE,
    th_DofUidRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
  string Tag_DofUidCol_name = "_DofUidCol_"+MoFEMFiniteElement_name;
  rval = moab.tag_get_handle(Tag_DofUidCol_name.c_str(),def_val_len,MB_TYPE_OPAQUE,
    th_DofUidCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
  string Tag_DofUidData_name = "_DofUidData_"+MoFEMFiniteElement_name;
  rval = moab.tag_get_handle(Tag_DofUidData_name.c_str(),def_val_len,MB_TYPE_OPAQUE,
    th_DofUidData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
  //add MoFEMFiniteElement
  pair<MoFEMFiniteElement_multiIndex::iterator,bool> p = finite_elements.insert(MoFEMFiniteElement(moab,meshset));
  if(!p.second) SETERRQ(PETSC_COMM_SELF,1,"MoFEMFiniteElement not inserted");
  if(verbose>0) {
    ostringstream ss;
    ss << "add finite element: " << MoFEMFiniteElement_name << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    //list_finite_elements();
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::modify_finite_element_add_field_data(const string &MoFEMFiniteElement_name,const string &name_data) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finite_elements_by_name;
  finite_elements_by_name &MoFEMFiniteElement_name_set = finite_elements.get<MoFEMFiniteElement_name_mi_tag>();
  finite_elements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,1,"this MoFEMFiniteElement is there");
  bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,EntMoFEMFiniteElement_change_bit_add(get_BitFieldId(name_data)));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::modify_finite_element_add_field_row(const string &MoFEMFiniteElement_name,const string &name_row) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finite_elements_by_name;
  finite_elements_by_name &MoFEMFiniteElement_name_set = finite_elements.get<MoFEMFiniteElement_name_mi_tag>();
  finite_elements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,1,"this MoFEMFiniteElement is there");
  bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_row_change_bit_add(get_BitFieldId(name_row)));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::modify_finite_element_add_field_col(const string &MoFEMFiniteElement_name,const string &name_col) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finite_elements_by_name;
  finite_elements_by_name &MoFEMFiniteElement_name_set = finite_elements.get<MoFEMFiniteElement_name_mi_tag>();
  finite_elements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,1,"this MoFEMFiniteElement is there");
  bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_col_change_bit_add(get_BitFieldId(name_col)));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  PetscFunctionReturn(0);
}
BitFEId moabField_Core::get_BitFEId(const string& name) const {
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finite_elements_by_name;
  const finite_elements_by_name& set = finite_elements.get<MoFEMFiniteElement_name_mi_tag>();
  finite_elements_by_name::iterator miit = set.find(name);
  if(miit==set.end()) THROW_AT_LINE(("finite element < "+name+" > not found (top tip: check spelling)").c_str());
  return miit->get_id();
}
string moabField_Core::get_BitFEId_name(const BitFEId id) const {
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finite_elements_by_id;
  const finite_elements_by_id& set = finite_elements.get<BitFEId_mi_tag>();
  finite_elements_by_id::iterator miit = set.find(id);
  assert(miit!=set.end());
  return miit->get_name();
}
EntityHandle moabField_Core::get_meshset_by_BitFEId(const BitFEId id) const {
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finite_elements_by_id;
  const finite_elements_by_id& set = finite_elements.get<BitFEId_mi_tag>();
  finite_elements_by_id::iterator miit = set.find(id);
  if(miit==set.end()) THROW_AT_LINE("finite element not found");
  return miit->meshset;
}
EntityHandle moabField_Core::get_meshset_by_BitFEId(const string& name) const {	
  return get_meshset_by_BitFEId(get_BitFEId(name));
}
PetscErrorCode moabField_Core::list_finite_elements() const {
  PetscFunctionBegin;
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finite_elements_by_id;
  const finite_elements_by_id &BitFEId_set = finite_elements.get<BitFEId_mi_tag>();
  finite_elements_by_id::iterator miit = BitFEId_set.begin();
  for(;miit!=BitFEId_set.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_problem(const BitProblemId id,const string& name) {
  PetscFunctionBegin;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
  rval = moab.tag_set_data(th_ProblemId,&meshset,1,&id); CHKERR_PETSC(rval);
  void const* tag_data[] = { name.c_str() };
  int tag_sizes[] = { name.size() };
  rval = moab.tag_set_by_ptr(th_ProblemName,&meshset,1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  //create entry
  pair<MoFEMProblem_multiIndex::iterator,bool> p = problems.insert(_MoFEMProblem_(moab,meshset));
  NOT_USED(p);
  assert(p.second);
  if(verbose>0) {
    ostringstream ss;
    ss << "add problem: " << name << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_problem(const string& name) {
  PetscFunctionBegin;
  BitProblemId id = get_problem_shift();
  ierr = add_problem(id,name); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
BitProblemId moabField_Core::get_BitProblemId(const string& name) const {
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  const problems_by_name& set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator miit = set.find(name);
  return miit->get_id();
}
PetscErrorCode moabField_Core::list_problem() const {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<BitProblemId_mi_tag>::type problem_set_by_id;
  const problem_set_by_id &set_id = problems.get<BitProblemId_mi_tag>();
  problem_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::erase_inactive_dofs_moabfield() {
  PetscFunctionBegin;
  DofMoFEMEntity_multiIndex_active_view set_active_view;
  for(DofMoFEMEntity_multiIndex::iterator miit = dofs_moabfield.begin();
    miit!=dofs_moabfield.end();miit++) {
      set_active_view.insert(&*miit);
    }
  DofMoFEMEntity_multiIndex_active_view::iterator miit2 = set_active_view.lower_bound(0);
  DofMoFEMEntity_multiIndex_active_view::iterator hi_miit2 = set_active_view.upper_bound(0);
  for(;miit2!=hi_miit2;miit2++) {
    if(verbose>1) {
      ostringstream ss;
      ss << "del: " << **miit2 << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    UId uid = (*miit2)->get_unique_id();
    dofs_moabfield.erase(uid);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  EntityHandle idm = no_handle;
  try {
    idm = get_meshset_by_BitFEId(id);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,false); CHKERR_PETSC(rval);
  rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_TETs(meshset,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit_ref,const string &name,EntityType type) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  typedef RefMoFEMElement_multiIndex::index<EntType_mi_tag>::type refMoabFE_by_type;
  refMoabFE_by_type &ref_MoFEMFiniteElement = refined_mofem_elements.get<EntType_mi_tag>();
  refMoabFE_by_type::iterator miit = ref_MoFEMFiniteElement.lower_bound(type);
  refMoabFE_by_type::iterator hi_miit = ref_MoFEMFiniteElement.upper_bound(type);
  int nb_add_FEs = 0;
  for(;miit!=hi_miit;miit++) {
    if((miit->get_BitRefLevel()&bit_ref)!=bit_ref) continue;
    EntityHandle ent = miit->get_ref_ent();
    rval = moab.add_entities(idm,&ent,1); CHKERR_PETSC(rval);
    nb_add_FEs++;
  }
  {
    ostringstream ss;
    ss << "Add Nb. FEs " << nb_add_FEs << " form BitRef " << bit_ref << endl;
    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const string& name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  rval = moab.add_entities(idm,&meshset,1); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_ents_to_finite_element_by_MESHSETs(const EntityHandle meshset,const string& name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  Range meshsets;
  rval = moab.get_entities_by_type(meshset,MBENTITYSET,meshsets,false); CHKERR_PETSC(rval);
  rval = moab.add_entities(idm,meshsets); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::modify_problem_add_finite_element(const string &name_problem,const string &MoFEMFiniteElement_name) {
  PetscFunctionBegin;
  try {
    typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
    problems_by_name& set = problems.get<MoFEMProblem_mi_tag>();
    problems_by_name::iterator miit = set.find(name_problem);
    if(miit==set.end()) SETERRQ(PETSC_COMM_SELF,1,"this problem is not there");
    BitFEId f_id = get_BitFEId(MoFEMFiniteElement_name);
    bool success = set.modify(miit,problem_MoFEMFiniteElement_change_bit_add(f_id));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::modify_problem_ref_level_add_bit(const string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  problems_by_name& set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator miit = set.find(name_problem);
  if(miit==set.end()) SETERRQ(PETSC_COMM_SELF,1,"this problem is there");
  bool success = set.modify(miit,problem_change_ref_level_bit_add(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::build_finite_element(const EntMoFEMFiniteElement &EntFe,int verb) {
  PetscFunctionBegin;
  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_by_id;
  typedef MoFEMField_multiIndex::index<Meshset_mi_tag>::type field_by_meshset;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ent_by_ent;
  typedef DofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type dof_set_type;
  field_by_id &moabfields_by_id = moabfields.get<BitFieldId_mi_tag>();
  field_by_meshset &moabfields_by_meshset = moabfields.get<Meshset_mi_tag>();
  dof_set_type& dof_set = dofs_moabfield.get<Composite_mi_tag2>();
  EntityHandle fe_ent = EntFe.get_ent();
  pair<EntMoFEMFiniteElement_multiIndex::iterator,bool> p = finite_elements_moabents.insert(EntFe);
  //get id of mofem fields for row, col and data
  enum IntLoop { Row = 0,Col,Data,Last };
  BitFieldId FEAdj_fields[Last] = { 
    EntFe.get_BitFieldId_row(), EntFe.get_BitFieldId_col(), EntFe.get_BitFieldId_data() 
  };
  //get refinment level
  const BitRefLevel& bit_ref_MoFEMFiniteElement = p.first->get_BitRefLevel();
  Range tets,faces,edges,nodes,meshsets,adj_ents,ent_ents;
  Range::iterator eit_eit;
  DofMoFEMEntity_multiIndex_uid_view MoFEMFiniteElement_dof_uid_view[Last];
  //lopp over all fields in database
  for(unsigned int ii = 0;ii<BitFieldId().size();ii++) {
    // common field id for Row, Col and Data
    BitFieldId id_common = 0;
    //check if the field (ii) is added to finite element
    for(int ss = 0;ss<Last;ss++) id_common |= FEAdj_fields[ss]&BitFieldId().set(ii);
    if( id_common.none() ) continue;
    //find in database data associated with the field (ii)
    field_by_id::iterator miit = moabfields_by_id.find(BitFieldId().set(ii));
    if(miit==moabfields_by_id.end()) SETERRQ(PETSC_COMM_SELF,1,"data incosistency");
    //get field (ii) space
    FieldSpace space = miit->get_space();
    //resolve antities on element
    switch (moab.type_from_handle(fe_ent)) {
      case MBTET:
	 switch (space) {
	  case H1: if(nodes.empty()) moab.get_connectivity(&fe_ent,1,nodes,true);
  	   adj_ents.insert(nodes.begin(),nodes.end());
  	  case Hdiv: if(edges.empty()) moab.get_adjacencies(&fe_ent,1,1,false,edges);
  	   adj_ents.insert(edges.begin(),edges.end());
	   for(Range::iterator eeit = edges.begin();eeit!=edges.end();eeit++) p.first->get_side_number_ptr(moab,*eeit);
  	  case Hcurl: if(faces.empty()) moab.get_adjacencies(&fe_ent,1,2,false,faces);
  	   adj_ents.insert(faces.begin(),faces.end());
	   for(Range::iterator fit = faces.begin();fit!=faces.end();fit++) p.first->get_side_number_ptr(moab,*fit);
  	  case L2:
  	   adj_ents.insert(fe_ent);
  	   break;
	  case NoField: {
	      EntityHandle field_meshset = miit->get_meshset();
	      adj_ents.insert(field_meshset);
	   }
	   break;
      	  default:
  	   SETERRQ(PETSC_COMM_SELF,1,"this fild is not implemented for TET finite element");
	 }
	 break;
      case MBPRISM: {
	  try {
	    EntityHandle prism = fe_ent;
	    EntityHandle face_side3,face_side4;
	    rval = moab.side_element(prism,2,3,face_side3); CHKERR_PETSC(rval);
	    rval = moab.side_element(prism,2,4,face_side4); CHKERR_PETSC(rval);
	    EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,face_side3);
	    EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,face_side4);
	    int ee = 0;
	    for(;ee<3;ee++) {
	      EntityHandle edge;
	      rval = moab.side_element(prism,1,ee,edge); CHKERR_PETSC(rval);
	      EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,edge);
	      rval = moab.side_element(prism,1,6+ee,edge); CHKERR_PETSC(rval);
	      EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,edge);
	    }
	  } catch (const char* msg) {
	      SETERRQ(PETSC_COMM_SELF,1,msg);
	  }
	  SideNumber_multiIndex &side_table = EntFe.get_RefMoFEMElement()->get_side_number_table();
	  switch (space) {
	    case H1: if(nodes.empty()) moab.get_connectivity(&fe_ent,1,nodes,true);
	      adj_ents.insert(nodes.begin(),nodes.end());
	    case Hdiv: {
	      SideNumber_multiIndex::nth_index<2>::type::iterator
		siit = side_table.get<2>().lower_bound(MBEDGE), hi_siit = side_table.get<2>().upper_bound(MBEDGE);
	      for(;siit!=hi_siit;siit++) adj_ents.insert(siit->ent);
	    }
	    case Hcurl: {
	      SideNumber_multiIndex::nth_index<2>::type::iterator
		siit = side_table.get<2>().lower_bound(MBTRI), hi_siit = side_table.get<2>().upper_bound(MBTRI);
	      for(;siit!=hi_siit;siit++) adj_ents.insert(siit->ent);
	    }
	    case L2:
	      adj_ents.insert(fe_ent); 
	      break;
	    default:
	      SETERRQ(PETSC_COMM_SELF,1,"this fild is not implemented for PRISM finite element");
	  }
	  }
         break;
      case MBENTITYSET: 
	 //get all meshsets in finite element meshset 
	 rval = moab.get_entities_by_type(fe_ent,MBENTITYSET,ent_ents,false); CHKERR_PETSC(rval);
	 //resolve recusively all ents in the meshset
	 rval = moab.get_entities_by_handle(fe_ent,ent_ents,true); CHKERR_PETSC(rval); 
	 eit_eit = ent_ents.begin();
	 for(;eit_eit!=ent_ents.end();eit_eit++) {
	  switch (space) {
	    case H1: 
	      if(moab.type_from_handle(*eit_eit)!=MBENTITYSET) {
		if(moab.dimension_from_handle(*eit_eit)>0) {
		  rval = moab.get_connectivity(&*eit_eit,1,nodes,true); CHKERR_PETSC(rval);
		  adj_ents.insert(nodes.begin(),nodes.end());
		} else if(moab.type_from_handle(*eit_eit)==MBVERTEX) adj_ents.insert(*eit_eit);
	      }
	    case Hdiv: 
	      if(moab.type_from_handle(*eit_eit)!=MBENTITYSET) {
		if(moab.dimension_from_handle(*eit_eit)>1) {
		  rval = moab.get_adjacencies(&*eit_eit,1,1,false,edges); CHKERR_PETSC(rval);
		  adj_ents.insert(edges.begin(),edges.end());
		} else if(moab.dimension_from_handle(*eit_eit)==1) {
		  adj_ents.insert(*eit_eit);
		}
	      }
	    case Hcurl: 
	      if(moab.type_from_handle(*eit_eit)!=MBENTITYSET) {
		if(moab.dimension_from_handle(*eit_eit)>=2) {
		  rval = moab.get_adjacencies(&*eit_eit,1,2,false,faces); CHKERR_PETSC(rval);
		  adj_ents.insert(faces.begin(),faces.end());
		} else if(moab.dimension_from_handle(*eit_eit)==2) {
		  adj_ents.insert(*eit_eit);
		}
	      }
	    case L2:
	      if(moab.type_from_handle(*eit_eit)!=MBENTITYSET) {
		if(moab.dimension_from_handle(*eit_eit)==3) {
		  adj_ents.insert(*eit_eit);
		}
	      }
	    break;
	    case NoField:
	      if(moab.type_from_handle(*eit_eit)==MBENTITYSET) {
		//if field (ii) has space NoField only add dofs which associated with the meshsets
		if(moabfields_by_meshset.find(*eit_eit)!=moabfields_by_meshset.end()) {
		  adj_ents.insert(*eit_eit);
		}
	      }
	    break;
	    default:
	      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	 }
	 break;
	 default:
	  SETERRQ(PETSC_COMM_SELF,1,"this fild is not implemented for this type finite element");
    }
    Range::iterator eit2 = adj_ents.begin();
    for(;eit2!=adj_ents.end();eit2++) {
      ref_ent_by_ent::iterator ref_ent_miit = refined_mofem_entities.get<MoABEnt_mi_tag>().find(*eit2);
      if(ref_ent_miit==refined_mofem_entities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"ref ent not in database"); 
      const BitRefLevel& bit_ref_ent = ref_ent_miit->get_BitRefLevel();
      if(!(bit_ref_MoFEMFiniteElement&bit_ref_ent).any()) {
	ostringstream ss;
	ss << "inconsitency in database" << " type " << moab.type_from_handle(*eit2) << " bits " << bit_ref_MoFEMFiniteElement << " " << bit_ref_ent;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      dof_set_type::iterator ents_miit2 = dof_set.lower_bound(boost::make_tuple(miit->get_name(),ref_ent_miit->get_ref_ent()));
      dof_set_type::iterator ents_hi_miit2 = dof_set.upper_bound(boost::make_tuple(miit->get_name(),ref_ent_miit->get_ref_ent()));
      for(int ss = 0;ss<Last;ss++) {
	if( !(FEAdj_fields[ss].test(ii)) ) continue;
	dof_set_type::iterator ents_miit3 = ents_miit2;
	for(;ents_miit3!=ents_hi_miit2;ents_miit3++) {
	  MoFEMFiniteElement_dof_uid_view[ss].insert(&*ents_miit3);
	}
      }
    }
  }
  try {
      for(int ss = 0;ss<Last;ss++) {
	bool success;
	switch (ss) {
    	  case Row: success = finite_elements_moabents.modify(p.first,EntMoFEMFiniteElement_row_dofs_change(moab,MoFEMFiniteElement_dof_uid_view[ss])); break;
      	  case Col: success = finite_elements_moabents.modify(p.first,EntMoFEMFiniteElement_col_dofs_change(moab,MoFEMFiniteElement_dof_uid_view[ss])); break;
    	  case Data: success = finite_elements_moabents.modify(p.first,EntMoFEMFiniteElement_data_dofs_change(moab,MoFEMFiniteElement_dof_uid_view[ss])); break;
    	  default: assert(0);
    	}
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  // build data_dofs
  const void* tag_data_uids_data = p.first->tag_data_uids_data;
  const int tag_data_uids_size = p.first->tag_data_uids_size;
  DofMoFEMEntity_multiIndex_active_view data_view;
  ierr = get_MoFEMFiniteElement_dof_uid_view(dofs_moabfield,data_view,Interface::UNION,tag_data_uids_data,tag_data_uids_size); CHKERRQ(ierr);
  DofMoFEMEntity_multiIndex_active_view::iterator viit_data = data_view.begin();
  for(;viit_data!=data_view.end();viit_data++) {
    try {
      switch((*viit_data)->get_space()) {
	case H1:
	case Hdiv:
	case Hcurl:
	case L2: 
	case NoField:
	{
	  SideNumber *side_number_ptr = p.first->get_side_number_ptr(moab,(*viit_data)->get_ent());
	  FEDofMoFEMEntity_multiIndex &data_dofs = const_cast<FEDofMoFEMEntity_multiIndex&>(p.first->data_dofs);
	  FEDofMoFEMEntity FEDof(side_number_ptr,&**viit_data);
	  //add dofs to finite element multi_index database
	  data_dofs.insert(FEDof);
	}
	break;
	default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,1,msg);
    }
  }
  if(verb>2) {
    ostringstream ss;
    ss << "add: FE data"  << endl << *p.first << endl;
    //rows
    DofMoFEMEntity_multiIndex_uid_view MoFEMFiniteElement_row_dof_uid_view;
    ierr = p.first->get_MoFEMFiniteElement_row_dof_uid_view(dofs_moabfield,MoFEMFiniteElement_row_dof_uid_view); CHKERRQ(ierr);
    DofMoFEMEntity_multiIndex_uid_view::iterator miit_row = MoFEMFiniteElement_row_dof_uid_view.begin();
    ss << "rows dofs" << endl;
    for(;miit_row!=MoFEMFiniteElement_row_dof_uid_view.end();miit_row++) ss << **miit_row << endl;
    //cols
    DofMoFEMEntity_multiIndex_uid_view MoFEMFiniteElement_col_dof_uid_view;
    ierr = p.first->get_MoFEMFiniteElement_col_dof_uid_view(dofs_moabfield,MoFEMFiniteElement_col_dof_uid_view); CHKERRQ(ierr);
    DofMoFEMEntity_multiIndex_uid_view::iterator miit_col = MoFEMFiniteElement_col_dof_uid_view.begin();
    ss << "cols dofs" << endl;
    for(;miit_col!=MoFEMFiniteElement_col_dof_uid_view.end();miit_col++) ss << **miit_col << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::build_finite_elements(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMElement_multiIndex::index<MoABEnt_mi_tag>::type ref_MoFEMFiniteElement_by_ent;
  MoFEMFiniteElement_multiIndex::iterator MoFEMFiniteElement_miit = finite_elements.begin();
  // loop Finite Elements
  for(;MoFEMFiniteElement_miit!=finite_elements.end();MoFEMFiniteElement_miit++) {
    if(verbose>0) PetscPrintf(PETSC_COMM_WORLD,"Build Finite Elements %s\n",MoFEMFiniteElement_miit->get_name().c_str());
    //get finite element meshset
    EntityHandle meshset = get_meshset_by_BitFEId(MoFEMFiniteElement_miit->get_id());
    // get entities from finite element meshset // if meshset 
    Range MoFEMFiniteElement_ents;
    rval = moab.get_entities_by_handle(meshset,MoFEMFiniteElement_ents,false); CHKERR_PETSC(rval);
    //loop meshset Ents and add finite elements
    Range::iterator eit = MoFEMFiniteElement_ents.begin();
    for(;eit!=MoFEMFiniteElement_ents.end();eit++) {
      // check if is in refined_mofem_elements database
      ref_MoFEMFiniteElement_by_ent::iterator ref_MoFEMFiniteElement_miit = refined_mofem_elements.get<MoABEnt_mi_tag>().find(*eit); /* iterator is a wrapper*/
      if(ref_MoFEMFiniteElement_miit == refined_mofem_elements.get<MoABEnt_mi_tag>().end()) {
	ostringstream ss;
	ss << "ref MoFEMFiniteElement not in database ent = " << *eit;
	ss << " type " << moab.type_from_handle(*eit);
	ss << " " << *MoFEMFiniteElement_miit;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      ierr = build_finite_element(EntMoFEMFiniteElement(moab,ref_MoFEMFiniteElement_miit->get_RefMoFEMElement(),&*MoFEMFiniteElement_miit),verb); CHKERRQ(ierr);
    }
  }
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Nb. FEs %u\n",finite_elements_moabents.size());
    typedef EntMoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type MoFEMFiniteElement_by_id;
    MoFEMFiniteElement_by_id &MoFEMFiniteElements = finite_elements_moabents.get<BitFEId_mi_tag>();  
    MoFEMFiniteElement_multiIndex::iterator id_MoFEMFiniteElement = finite_elements.begin();
    for(;id_MoFEMFiniteElement!=finite_elements.end();id_MoFEMFiniteElement++) {
      MoFEMFiniteElement_by_id::iterator miit = MoFEMFiniteElements.lower_bound(id_MoFEMFiniteElement->get_id());
      MoFEMFiniteElement_by_id::iterator hi_miit = MoFEMFiniteElements.upper_bound(id_MoFEMFiniteElement->get_id());
      int count = std::distance(miit,hi_miit);
      ostringstream ss;
      ss << *id_MoFEMFiniteElement << " Nb. FEs " << count << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
  }
  *build_MoFEM |= 1<<1;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::build_adjacencies(const BitRefLevel bit) {
  PetscFunctionBegin;
  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,1,"field not build");
  if(!(*build_MoFEM)&(1<<1)) SETERRQ(PETSC_COMM_SELF,1,"fe not build");
  typedef MoFEMEntity_multiIndex::index<Unique_mi_tag>::type ents_by_uid;
  EntMoFEMFiniteElement_multiIndex::iterator fit = finite_elements_moabents.begin();
  for(;fit!=finite_elements_moabents.end();fit++) {
    if(!(fit->get_BitRefLevel()&bit).any()) continue;
    int size_row = fit->tag_row_uids_size/sizeof(UId);
    const UId *uids_row = (UId*)fit->tag_row_uids_data;
    int ii = 0;
    UId uid = 0;
    for(;ii<size_row;ii++) {
      if( uid == (uids_row[ii] >> 8 )) continue;
      uid = uids_row[ii];
      uid = uid >> 8; //look to DofMoFEMEntity::get_unique_id_calculate and MoFEMEntity::get_unique_id_calculate() <- uid is shifted by 8 bits
      ents_by_uid::iterator miit = ents_moabfield.get<Unique_mi_tag>().find(uid);
      assert(dofs_moabfield.get<Unique_mi_tag>().find(uids_row[ii])!=dofs_moabfield.get<Unique_mi_tag>().end());
      assert(dofs_moabfield.get<Unique_mi_tag>().find(uids_row[ii])->get_MoFEMEntity_ptr()->get_unique_id()==uid);
      if(miit ==ents_moabfield.get<Unique_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      pair<MoFEMAdjacencies_multiIndex::iterator,bool> p = adjacencies.insert(MoFEMAdjacencies(&*miit,&*fit));
      bool success = adjacencies.modify(p.first,MoFEMAdjacencies_change_by_what(by_row));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    int size_col = fit->tag_col_uids_size/sizeof(UId);
    const UId *uids_col = (UId*)fit->tag_col_uids_data;
    for(ii = 0,uid = 0;ii<size_col;ii++) {
      if( uid == (uids_col[ii] >> 8 )) continue;
      uid = uids_col[ii];
      uid = uid >> 8; //look to DofMoFEMEntity::get_unique_id_calculate and MoFEMEntity::get_unique_id_calculate() <- uid is shifted by 8 bits
      assert(dofs_moabfield.get<Unique_mi_tag>().find(uids_col[ii])!=dofs_moabfield.get<Unique_mi_tag>().end());
      assert(dofs_moabfield.get<Unique_mi_tag>().find(uids_col[ii])->get_MoFEMEntity_ptr()->get_unique_id()==uid);
      ents_by_uid::iterator miit = ents_moabfield.get<Unique_mi_tag>().find(uid);
      if(miit ==ents_moabfield.get<Unique_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      pair<MoFEMAdjacencies_multiIndex::iterator,bool> p = adjacencies.insert(MoFEMAdjacencies(&*miit,&*fit));
      bool success = adjacencies.modify(p.first,MoFEMAdjacencies_change_by_what(by_col));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    int size_data = fit->tag_data_uids_size/sizeof(UId);
    const UId *uids_data = (UId*)fit->tag_data_uids_data;
    for(ii = 0,uid = 0;ii<size_data;ii++) {
      if( uid == (uids_data[ii] >> 8 )) continue;
      uid = uids_data[ii];
      uid = uid >> 8; //look to DofMoFEMEntity::get_unique_id_calculate and MoFEMEntity::get_unique_id_calculate() <- uid is shifted by 8 bits
      assert(dofs_moabfield.get<Unique_mi_tag>().find(uids_data[ii])!=dofs_moabfield.get<Unique_mi_tag>().end());
      assert(dofs_moabfield.get<Unique_mi_tag>().find(uids_data[ii])->get_MoFEMEntity_ptr()->get_unique_id()==uid);
      ents_by_uid::iterator miit = ents_moabfield.get<Unique_mi_tag>().find(uid);
      if(miit == ents_moabfield.get<Unique_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      pair<MoFEMAdjacencies_multiIndex::iterator,bool> p = adjacencies.insert(MoFEMAdjacencies(&*miit,&*fit));
      bool success = adjacencies.modify(p.first,MoFEMAdjacencies_change_by_what(by_data));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
  }
  if(verbose>1) {
    list_adjacencies();
  }
  if(verbose>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Nb. adjacencies %u\n",adjacencies.size());
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Nb. adjacencies %u\n",adjacencies.size());
    //PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  *build_MoFEM |= 1<<2;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::list_adjacencies() const {
  PetscFunctionBegin;
  MoFEMAdjacencies_multiIndex::iterator miit = adjacencies.begin();
  for(;miit!=adjacencies.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::build_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"adjacencies not build");
  MoFEMProblem_multiIndex::iterator p_miit = problems.begin();
  for(;p_miit!=problems.end();p_miit++) {
    //miit2 iterator for finite elements
    EntMoFEMFiniteElement_multiIndex::iterator miit2 = finite_elements_moabents.begin();
    EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finite_elements_moabents.end();
    DofMoFEMEntity_multiIndex_uid_view dofs_rows;
    DofMoFEMEntity_multiIndex_uid_view dofs_cols;;
    EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
    for(;miit3!=hi_miit2;miit3++) {
      if((miit3->get_id()&p_miit->get_BitFEId()).any()) {
	if((miit3->get_BitRefLevel()&p_miit->get_BitRefLevel())==p_miit->get_BitRefLevel()) {
	  ierr = miit3->get_MoFEMFiniteElement_row_dof_uid_view(dofs_moabfield,dofs_rows); CHKERRQ(ierr);
	  ierr = miit3->get_MoFEMFiniteElement_col_dof_uid_view(dofs_moabfield,dofs_cols); CHKERRQ(ierr);
	}
      }
    }
    //zero rows
    bool success = problems.modify(p_miit,problem_zero_nb_rows_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    //zero cols
    success = problems.modify(p_miit,problem_zero_nb_cols_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    //add dofs for rows
    DofMoFEMEntity_multiIndex_uid_view::iterator miit4 = dofs_rows.begin();
    for(;miit4!=dofs_rows.end();miit4++) {
      if(!(*miit4)->active) continue;
      success = problems.modify(p_miit,problem_row_change(&**miit4));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    //add dofs for cols
    DofMoFEMEntity_multiIndex_uid_view::iterator miit5 = dofs_cols.begin();
    for(;miit5!=dofs_cols.end();miit5++) {
      if(!(*miit5)->active) continue;
      success = problems.modify(p_miit,problem_col_change(&**miit5));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    if(verbose>0) {
      PetscPrintf(PETSC_COMM_WORLD,"Problem %s Nb. rows %u Nb. cols %u\n",p_miit->get_name().c_str(),
	p_miit->numered_dofs_rows.size(),p_miit->numered_dofs_cols.size());
    }
    if(verb>2) {
      EntMoFEMFiniteElement_multiIndex::iterator miit_ss = miit2;
      ostringstream ss;
      ss << "rank " << pcomm->rank() << " ";
      ss << "FEs data for problem " << *p_miit << endl;
      for(;miit_ss!=hi_miit2;miit_ss++) {
	ss << "rank " << pcomm->rank() << " ";
	ss << *miit_ss << endl;
      }
      ss << "rank " << pcomm->rank() << " ";
      ss << "FEs row dofs "<< *p_miit << " Nb. row dof " << p_miit->get_nb_dofs_row() << endl;
      NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_row = p_miit->numered_dofs_rows.begin();
      for(;miit_dd_row!=p_miit->numered_dofs_rows.end();miit_dd_row++) {
	ss << "rank " << pcomm->rank() << " ";
	ss<<*miit_dd_row<<endl;
      }
      ss << "rank " << pcomm->rank() << " ";
      ss << "FEs col dofs "<< *p_miit << " Nb. col dof " << p_miit->get_nb_dofs_col() << endl;
      NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols.begin();
      for(;miit_dd_col!=p_miit->numered_dofs_cols.end();miit_dd_col++) {
	ss << "rank " << pcomm->rank() << " ";
	ss<<*miit_dd_col<<endl;
      }
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      PetscSynchronizedFlush(PETSC_COMM_WORLD); 
    }
  }
  *build_MoFEM |= 1<<3;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::partition_problem(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"adjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,1,"problems not build");
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef NumeredDofMoFEMEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  //find p_miit
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(name);
  if(p_miit==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",name.c_str());
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Partition problem %s\n",p_miit->get_name().c_str());
  }
  DofIdx nb_dofs_row = p_miit->get_nb_dofs_row();
  Mat Adj;
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"\tcreate Adj matrix\n");
  }
  ierr = partition_create_Mat<Idx_mi_tag>(name,&Adj,NULL,true,verb); CHKERRQ(ierr);
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"\t<- done\n");
  }
  //PetscBarrier(PETSC_NULL);
  PetscInt m,n;
  ierr = MatGetSize(Adj,&m,&n); CHKERRQ(ierr);
  if(m!=p_miit->get_nb_dofs_row()) SETERRQ(PETSC_COMM_SELF,1,"row number inconsistency");
  if(n!=p_miit->get_nb_dofs_col()) SETERRQ(PETSC_COMM_SELF,1,"col number inconsistency");
  if(verb>2) {
    MatView(Adj,PETSC_VIEWER_STDOUT_WORLD);
  }
  PetscInt *_values;
  PetscMalloc(nb_dofs_row*sizeof(PetscInt),&_values);
  fill(&_values[0],&_values[nb_dofs_row],1);
  //partitioning
  MatPartitioning part;
  IS is;
  ierr = MatPartitioningCreate(MPI_COMM_WORLD,&part); CHKERRQ(ierr);
  ierr = MatPartitioningSetAdjacency(part,Adj); CHKERRQ(ierr);
  ierr = MatPartitioningSetVertexWeights(part,_values); CHKERRQ(ierr);
  ierr = MatPartitioningSetFromOptions(part); CHKERRQ(ierr);
  ierr = MatPartitioningApply(part,&is); CHKERRQ(ierr);
  if(verb>2) {
    ISView(is,PETSC_VIEWER_STDOUT_WORLD);
  }
  //gather
  IS is_gather,is_num,is_gather_num;
  ierr = ISAllGather(is,&is_gather); CHKERRQ(ierr);
  ierr = ISPartitioningToNumbering(is,&is_num); CHKERRQ(ierr);
  ierr = ISAllGather(is_num,&is_gather_num); CHKERRQ(ierr);
  const int *part_number,*petsc_idx;
  ierr = ISGetIndices(is_gather,&part_number);  CHKERRQ(ierr);
  ierr = ISGetIndices(is_gather_num,&petsc_idx);  CHKERRQ(ierr);
  PetscInt size_is_num,size_is_gather;
  ISGetSize(is_gather,&size_is_gather);
  assert(size_is_gather == (int)nb_dofs_row);
  ISGetSize(is_num,&size_is_num);
  assert(size_is_num == (int)nb_dofs_row);
  //set petsc global indicies
  NumeredDofMoFEMEntitys_by_idx &dofs_row_by_idx_no_const = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_rows.get<Idx_mi_tag>());
  NumeredDofMoFEMEntitys_by_idx &dofs_col_by_idx_no_const = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_cols.get<Idx_mi_tag>());
  DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
  DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
  nb_row_local_dofs = 0;
  nb_col_local_dofs = 0;
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  nb_row_ghost_dofs = 0;
  nb_col_ghost_dofs = 0;
  NumeredDofMoFEMEntitys_by_idx::iterator miit_dofs_row = dofs_row_by_idx_no_const.begin();
  NumeredDofMoFEMEntitys_by_idx::iterator miit_dofs_col = dofs_col_by_idx_no_const.begin();
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"\tloop problem dofs");
  }
  for(;miit_dofs_row!=dofs_row_by_idx_no_const.end();miit_dofs_row++,miit_dofs_col++) {
    if(miit_dofs_col==dofs_col_by_idx_no_const.end()) SETERRQ(PETSC_COMM_SELF,1,"check finite element definition, nb. of rows is not equal to number for columns");
    if(miit_dofs_row->get_unique_id()!=miit_dofs_col->get_unique_id()) SETERRQ(PETSC_COMM_SELF,1,"check finite element definition, rows not equal columns");
    if(miit_dofs_row->dof_idx!=miit_dofs_col->dof_idx) SETERRQ(PETSC_COMM_SELF,1,"check finite element definition, rows not equal columns");
    assert(petsc_idx[miit_dofs_row->dof_idx]>=0);
    assert(petsc_idx[miit_dofs_row->dof_idx]<(int)p_miit->get_nb_dofs_row());
    bool success = dofs_row_by_idx_no_const.modify(miit_dofs_row,NumeredDofMoFEMEntity_part_change(part_number[miit_dofs_row->dof_idx],petsc_idx[miit_dofs_row->dof_idx]));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    success = dofs_col_by_idx_no_const.modify(miit_dofs_col,NumeredDofMoFEMEntity_part_change(part_number[miit_dofs_col->dof_idx],petsc_idx[miit_dofs_col->dof_idx]));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    if(miit_dofs_row->part == pcomm->rank()) {
      assert(miit_dofs_row->part==miit_dofs_col->part);
      assert(miit_dofs_row->petsc_gloabl_dof_idx==miit_dofs_col->petsc_gloabl_dof_idx);
      success = dofs_row_by_idx_no_const.modify(miit_dofs_row,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      success = dofs_col_by_idx_no_const.modify(miit_dofs_col,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
  }
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD," <- done\n");
  }
  if(verbose>0) {
    ostringstream ss;
    ss << "partition_problem: rank = " << pcomm->rank() << " FEs row ghost dofs "<< *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_row() << " nb global row dofs " << p_miit->get_nb_dofs_row() << endl;
    ss << "partition_problem: rank = " << pcomm->rank() << " FEs col ghost dofs " << *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_col() << " nb global col dofs " << p_miit->get_nb_dofs_col() << endl;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  if(verb>2) {
    ostringstream ss;
    ss << "rank = " << pcomm->rank() << " FEs row dofs "<< *p_miit << " Nb. row dof " << p_miit->get_nb_dofs_row() 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_row() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_row = p_miit->numered_dofs_rows.begin();
    for(;miit_dd_row!=p_miit->numered_dofs_rows.end();miit_dd_row++) {
	ss<<*miit_dd_row<<endl;
    }
    ss << "rank = " << pcomm->rank() << " FEs col dofs "<< *p_miit << " Nb. col dof " << p_miit->get_nb_dofs_col() 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_col() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols.begin();
    for(;miit_dd_col!=p_miit->numered_dofs_cols.end();miit_dd_col++) {
	ss<<*miit_dd_col<<endl;
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  //
  ierr = ISRestoreIndices(is_gather,&part_number);  CHKERRQ(ierr);
  ierr = ISRestoreIndices(is_gather_num,&petsc_idx);  CHKERRQ(ierr);
  ierr = ISDestroy(&is_num); CHKERRQ(ierr);
  ierr = ISDestroy(&is_gather_num); CHKERRQ(ierr);
  ierr = ISDestroy(&is_gather); CHKERRQ(ierr);
  ierr = ISDestroy(&is); CHKERRQ(ierr);
  ierr = MatPartitioningDestroy(&part); CHKERRQ(ierr);
  ierr = MatDestroy(&Adj); CHKERRQ(ierr);
  //ierr = PetscFree(_values); CHKERRQ(ierr);
  *build_MoFEM |= 1<<4;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::compose_problem(const string &name,const string &problem_for_rows,const string &problem_for_cols,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"adjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,1,"problems not build");
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofMoFEMEntitys_by_uid;
  //find p_miit
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(name);
  if(p_miit==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",name.c_str());
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Partition problem %s\n",p_miit->get_name().c_str());
  }
  problems_by_name::iterator p_miit_row = problems_set.find(problem_for_rows);
  if(p_miit_row==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",problem_for_rows.c_str());
  problems_by_name::iterator p_miit_col = problems_set.find(problem_for_cols);
  if(p_miit_col==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",problem_for_cols.c_str());
  const NumeredDofMoFEMEntity_multiIndex &dofs_row = p_miit_row->numered_dofs_rows;
  const NumeredDofMoFEMEntity_multiIndex &dofs_col = p_miit_col->numered_dofs_rows;
  MoFEMEntity *MoFEMEntity_ptr = NULL;
  //do rows
  NumeredDofMoFEMEntity_multiIndex::iterator miit_row = dofs_row.begin();
  NumeredDofMoFEMEntity_multiIndex::iterator hi_miit_row = dofs_row.end();
  map<DofIdx,const NumeredDofMoFEMEntity*> rows_problem_map;
  for(;miit_row!=hi_miit_row;miit_row++) {
    if( (MoFEMEntity_ptr == NULL) ? 1 : (MoFEMEntity_ptr->get_unique_id() != miit_row->field_ptr->field_ptr->get_unique_id()) ) {
      MoFEMEntity_ptr = const_cast<MoFEMEntity*>(miit_row->field_ptr->field_ptr);
      typedef MoFEMAdjacencies_multiIndex::index<Composite_mi_tag>::type adj_by_ent;
      adj_by_ent::iterator adj_miit = adjacencies.get<Composite_mi_tag>().lower_bound(boost::make_tuple(MoFEMEntity_ptr->get_meshset(),MoFEMEntity_ptr->get_ent()));
      adj_by_ent::iterator hi_adj_miit = adjacencies.get<Composite_mi_tag>().upper_bound(boost::make_tuple(MoFEMEntity_ptr->get_meshset(),MoFEMEntity_ptr->get_ent()));
      for(;adj_miit!=hi_adj_miit;adj_miit++) {
	if(!(adj_miit->by_other&by_row)) continue; // if it is not row if element
	if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) continue; // if element is not part of prblem
	if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&miit_row->get_BitRefLevel()).none()) continue; // if entity is not problem refinment level
	int size  = adj_miit->EntMoFEMFiniteElement_ptr->tag_row_uids_size/sizeof(UId);
	for(int ii = 0;ii<size;ii++) {
	  UId uid = adj_miit->EntMoFEMFiniteElement_ptr->tag_row_uids_data[ii];
	  NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_row.get<Unique_mi_tag>().find(uid);
	  if(pr_dof == dofs_row.get<Unique_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  DofIdx petsc_global_idx = pr_dof->get_petsc_gloabl_dof_idx();
	  rows_problem_map[petsc_global_idx] = &*pr_dof;
	}
      }
    }
  }
  //do cols
  MoFEMEntity_ptr == NULL;
  NumeredDofMoFEMEntity_multiIndex::iterator miit_col = dofs_col.begin();
  NumeredDofMoFEMEntity_multiIndex::iterator hi_miit_col = dofs_col.end();
  map<DofIdx,const NumeredDofMoFEMEntity*> cols_problem_map;
  for(;miit_col!=hi_miit_col;miit_col++) {
    if( (MoFEMEntity_ptr == NULL) ? 1 : (MoFEMEntity_ptr->get_unique_id() != miit_col->field_ptr->field_ptr->get_unique_id()) ) {
      MoFEMEntity_ptr = const_cast<MoFEMEntity*>(miit_col->field_ptr->field_ptr);
      typedef MoFEMAdjacencies_multiIndex::index<Composite_mi_tag>::type adj_by_ent;
      adj_by_ent::iterator adj_miit = adjacencies.get<Composite_mi_tag>().lower_bound(boost::make_tuple(MoFEMEntity_ptr->get_meshset(),MoFEMEntity_ptr->get_ent()));
      adj_by_ent::iterator hi_adj_miit = adjacencies.get<Composite_mi_tag>().upper_bound(boost::make_tuple(MoFEMEntity_ptr->get_meshset(),MoFEMEntity_ptr->get_ent()));
      for(;adj_miit!=hi_adj_miit;adj_miit++) {
	if(!(adj_miit->by_other&by_col)) continue; // if it is not row if element
	if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) continue; // if element is not part of prblem
	if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&miit_col->get_BitRefLevel()).none()) continue; // if entity is not problem refinment level
	int size  = adj_miit->EntMoFEMFiniteElement_ptr->tag_col_uids_size/sizeof(UId);
	for(int ii = 0;ii<size;ii++) {	  
	  UId uid = adj_miit->EntMoFEMFiniteElement_ptr->tag_col_uids_data[ii];
	  NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_col.get<Unique_mi_tag>().find(uid);
	  if(pr_dof == dofs_col.get<Unique_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  DofIdx petsc_global_idx = pr_dof->get_petsc_gloabl_dof_idx();
	  cols_problem_map[petsc_global_idx] = &*pr_dof;
	}
      }
    }
  }
  // build indices 
  NumeredDofMoFEMEntitys_by_uid &dofs_row_by_uid = const_cast<NumeredDofMoFEMEntitys_by_uid&>(p_miit->numered_dofs_rows.get<Unique_mi_tag>());
  NumeredDofMoFEMEntitys_by_uid &dofs_col_by_uid = const_cast<NumeredDofMoFEMEntitys_by_uid&>(p_miit->numered_dofs_cols.get<Unique_mi_tag>());
  DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
  DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
  nb_row_local_dofs = 0;
  nb_col_local_dofs = 0;
  //rows
  map<DofIdx,const NumeredDofMoFEMEntity*>::iterator miit_map_row = rows_problem_map.begin();
  for(;miit_map_row!=rows_problem_map.end();miit_map_row++) {
    NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_row_by_uid.find(miit_map_row->second->get_unique_id());
    if(pr_dof == dofs_row_by_uid.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    int part_number = miit_map_row->second->get_part();
    int petsc_global_dof = distance(rows_problem_map.begin(),miit_map_row);
    if(verb>1) {
      PetscPrintf(PETSC_COMM_WORLD,"Row Problem Glob Idx %d Problem Glob Idx %d\n",miit_map_row->second->get_petsc_gloabl_dof_idx(),petsc_global_dof);
    }
    bool success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    if(pr_dof->part == pcomm->rank()) {
      success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
  }
  //cols
  map<DofIdx,const NumeredDofMoFEMEntity*>::iterator miit_map_col = cols_problem_map.begin();
  for(;miit_map_col!=cols_problem_map.end();miit_map_col++) {
    NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_col_by_uid.find(miit_map_col->second->get_unique_id());
    if(pr_dof == dofs_col_by_uid.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    int part_number = miit_map_col->second->get_part();
    int petsc_global_dof = distance(cols_problem_map.begin(),miit_map_col);
    if(verb>1) {
      PetscPrintf(PETSC_COMM_WORLD,"Col Problem Glob Idx %d Problem Glob Idx %d\n",miit_map_col->second->get_petsc_gloabl_dof_idx(),petsc_global_dof);
    }
    bool success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    if(pr_dof->part == pcomm->rank()) {
      success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
  }
  if(verbose>0) {
    ostringstream ss;
    ss << "partition_problem: rank = " << pcomm->rank() << " FEs row ghost dofs "<< *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_row() << " nb global row dofs " << p_miit->get_nb_dofs_row() << endl;
    ss << "partition_problem: rank = " << pcomm->rank() << " FEs col ghost dofs " << *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_col() << " nb global col dofs " << p_miit->get_nb_dofs_col() << endl;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  if(verb>2) {
    ostringstream ss;
    ss << "rank = " << pcomm->rank() << " FEs row dofs "<< *p_miit << " Nb. row dof " << p_miit->get_nb_dofs_row() 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_row() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_row = p_miit->numered_dofs_rows.begin();
    for(;miit_dd_row!=p_miit->numered_dofs_rows.end();miit_dd_row++) {
	ss<<*miit_dd_row<<endl;
    }
    ss << "rank = " << pcomm->rank() << " FEs col dofs "<< *p_miit << " Nb. col dof " << p_miit->get_nb_dofs_col() 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_col() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols.begin();
    for(;miit_dd_col!=p_miit->numered_dofs_cols.end();miit_dd_col++) {
	ss<<*miit_dd_col<<endl;
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::partition_ghost_dofs(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"adjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,1,"problems not build");
  if(!(*build_MoFEM&(1<<4))) SETERRQ(PETSC_COMM_SELF,1,"partitions problems not build");
  if(!(*build_MoFEM&(1<<5))) SETERRQ(PETSC_COMM_SELF,1,"partitions finite elements not build");
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef NumeredDofMoFEMEntity_multiIndex::index<Part_mi_tag>::type NumeredDofMoFEMEntitys_by_part;
  typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofMoFEMEntitys_by_unique_id;
  typedef MoFEMAdjacencies_multiIndex::index<Composite_mi_tag>::type adj_by_ent;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  //find p_miit
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(name);
  //
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  //
  nb_row_ghost_dofs = 0;
  nb_col_ghost_dofs = 0;
  if(pcomm->size()>1) {
    NumeredDofMoFEMEntity_multiIndex_uid_view idx_view;
    NumeredMoFEMFiniteElement_multiIndex::index<Composite_unique_mi_tag>::type &numered_finite_elements 
      = const_cast<NumeredMoFEMFiniteElement_multiIndex::index<Composite_unique_mi_tag>::type&>(p_miit->numered_finite_elements.get<Composite_unique_mi_tag>());
    adj_by_ent& adj_by_Composite_mi_tag = adjacencies.get<Composite_mi_tag>();
    NumeredDofMoFEMEntitys_by_part *dof_by_part_no_const[2] = {
      const_cast<NumeredDofMoFEMEntitys_by_part*>(&p_miit->numered_dofs_cols.get<Part_mi_tag>()),
      const_cast<NumeredDofMoFEMEntitys_by_part*>(&p_miit->numered_dofs_rows.get<Part_mi_tag>()) 
    };
    NumeredDofMoFEMEntity_multiIndex *numered_dofs[2]= { 
      const_cast<NumeredDofMoFEMEntity_multiIndex*>(&p_miit->numered_dofs_cols),
      const_cast<NumeredDofMoFEMEntity_multiIndex*>(&p_miit->numered_dofs_rows) 
    };
    NumeredDofMoFEMEntitys_by_unique_id *dof_by_uid_no_const[2] = {
      const_cast<NumeredDofMoFEMEntitys_by_unique_id*>(&p_miit->numered_dofs_cols.get<Unique_mi_tag>()),
      const_cast<NumeredDofMoFEMEntitys_by_unique_id*>(&p_miit->numered_dofs_rows.get<Unique_mi_tag>()) 
    };
    by_what by[2] = { by_row,by_col };
    by_what by_other[2] = { by_col, by_row };
    NumeredDofMoFEMEntity_multiIndex_uid_view ghost_idx_col_view,ghost_idx_row_view;
    NumeredDofMoFEMEntity_multiIndex_uid_view *ghost_idx_view[2] = { &ghost_idx_col_view, &ghost_idx_row_view };
    DofIdx nb_local_dofs[2] = { *((DofIdx*)p_miit->tag_local_nbdof_data_col), *((DofIdx*)p_miit->tag_local_nbdof_data_row) };
    DofIdx *nb_ghost_dofs[2] = { &nb_col_ghost_dofs, &nb_row_ghost_dofs };
    MoFEMEntity *MoFEMEntity_ptr = NULL;
    for(int ss = 0; ss<2; ss++) {
      NumeredDofMoFEMEntitys_by_part::iterator miit_dof_by_part = dof_by_part_no_const[ss]->lower_bound(pcomm->rank());
      NumeredDofMoFEMEntitys_by_part::iterator hi_miit_dof_by_part = dof_by_part_no_const[ss]->upper_bound(pcomm->rank());
      for(;miit_dof_by_part!=hi_miit_dof_by_part;miit_dof_by_part++) {
	//cerr << *miit_dof_by_part << endl;
        if( (MoFEMEntity_ptr == NULL) ? 1 : (MoFEMEntity_ptr->get_unique_id() != miit_dof_by_part->get_MoFEMEntity_ptr()->get_unique_id()) ) {
	  MoFEMEntity_ptr = const_cast<MoFEMEntity *>(miit_dof_by_part->get_MoFEMEntity_ptr());
	  adj_by_ent::iterator adj_miit = adj_by_Composite_mi_tag.lower_bound(boost::make_tuple(MoFEMEntity_ptr->get_meshset(),MoFEMEntity_ptr->get_ent()));
	  adj_by_ent::iterator hi_adj_miit = adj_by_Composite_mi_tag.upper_bound(boost::make_tuple(MoFEMEntity_ptr->get_meshset(),MoFEMEntity_ptr->get_ent()));
	  idx_view.clear();
	  for(;adj_miit!=hi_adj_miit;adj_miit++) {
	    if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) continue;
	    if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&p_miit->get_BitRefLevel()).none()) continue;
	    if(numered_finite_elements.
	      find(boost::make_tuple(adj_miit->get_MoFEMFiniteElement_meshset(),adj_miit->get_MoFEMFiniteElement_entity_handle()))->part != pcomm->rank()) continue;
	    if(adj_miit->by_other&by_row_col) {
	      ierr = adj_miit->get_ent_adj_dofs_bridge(*(numered_dofs[ss]),by_other[ss],idx_view); CHKERRQ(ierr);
	    }
	  }
	  NumeredDofMoFEMEntity_multiIndex_uid_view::iterator adj_idx_miit = idx_view.begin();
	  for(;adj_idx_miit!=idx_view.end();adj_idx_miit++) {
	    if((*adj_idx_miit)->part==pcomm->rank()) continue;
	    ghost_idx_view[ss]->insert(*adj_idx_miit);
	  }
        }
      }
      NumeredDofMoFEMEntity_multiIndex_uid_view::iterator ghost_idx_miit = ghost_idx_view[ss]->begin();
      for(;ghost_idx_miit!=ghost_idx_view[ss]->end();ghost_idx_miit++) {
        NumeredDofMoFEMEntitys_by_unique_id::iterator miit4 = dof_by_uid_no_const[ss]->find((*ghost_idx_miit)->get_unique_id());
        if(miit4->petsc_local_dof_idx!=(DofIdx)-1) SETERRQ(PETSC_COMM_SELF,1,"inconsitent data");
        bool success = dof_by_uid_no_const[ss]->modify(miit4,NumeredDofMoFEMEntity_local_idx_change(nb_local_dofs[ss]++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
        (*nb_ghost_dofs[ss])++;
      }
    }
  }
  if(verb>0) {
    ostringstream ss;
    ss << "partition_ghost_dofs: rank = " << pcomm->rank() 
      << " FEs col ghost dofs "<< *p_miit 
      << " Nb. col ghost dof " << p_miit->get_nb_ghost_dofs_col() 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_col() << endl;
    ss << "partition_ghost_row_dofs: rank = " << pcomm->rank() 
      << " FEs row ghost dofs "<< *p_miit 
      << " Nb. row ghost dof " << p_miit->get_nb_ghost_dofs_row() 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_row() << endl;
    if(verb>1) {
      NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols.begin();
      for(;miit_dd_col!=p_miit->numered_dofs_cols.end();miit_dd_col++) {
	if(miit_dd_col->part==pcomm->rank()) continue;
	if(miit_dd_col->petsc_local_dof_idx==(DofIdx)-1) continue;
	ss<<*miit_dd_col<<endl;
      }
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  *build_MoFEM |= 1<<6;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::partition_finite_elements(const string &name,bool do_skip,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"adjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,1,"partitions not build");
  if(!(*build_MoFEM&(1<<4))) SETERRQ(PETSC_COMM_SELF,1,"partitions problems not build");
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  //find p_miit
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(name);
  if(p_miit == problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str());
  NumeredMoFEMFiniteElement_multiIndex& numered_finite_elements = const_cast<NumeredMoFEMFiniteElement_multiIndex&>(p_miit->numered_finite_elements);
  //MoFEMFiniteElement set
  EntMoFEMFiniteElement_multiIndex::iterator miit2 = finite_elements_moabents.begin();
  EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finite_elements_moabents.end();
  EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
  for(;miit3!=hi_miit2;miit3++) {
    if((miit3->get_BitRefLevel()&p_miit->get_BitRefLevel()).none()) continue;
    if((miit3->get_id()&p_miit->get_BitFEId()).any()) {
      NumeredDofMoFEMEntity_multiIndex_uid_view rows_view,cols_view;
      const void* tag_row_uids_data = miit3->tag_row_uids_data;
      const int tag_row_uids_size = miit3->tag_row_uids_size;
      ierr = get_MoFEMFiniteElement_dof_uid_view(p_miit->numered_dofs_rows,rows_view,Interface::UNION,tag_row_uids_data,tag_row_uids_size); CHKERRQ(ierr);
      pair<NumeredMoFEMFiniteElement_multiIndex::iterator,bool> p = numered_finite_elements.insert(_NumeredMoFEMFiniteElement_(&*miit3));
      _NumeredMoFEMFiniteElement_ &problem_MoFEMFiniteElement = const_cast<_NumeredMoFEMFiniteElement_&>(*p.first);
      if(!p.second) {
	problem_MoFEMFiniteElement.rows_dofs.clear();
	problem_MoFEMFiniteElement.cols_dofs.clear();
      }
      bool skip = true;
      NumeredDofMoFEMEntity_multiIndex_uid_view::iterator viit_rows = rows_view.begin();
      vector<int> parts(pcomm->size(),0);
      for(;viit_rows!=rows_view.end();viit_rows++) {
	try {
	  SideNumber *side_number_ptr = p.first->get_side_number_ptr(moab,(*viit_rows)->get_ent());
	  FENumeredDofMoFEMEntity_multiIndex &rows_dofs = const_cast<FENumeredDofMoFEMEntity_multiIndex&>(p.first->rows_dofs);
	  FENumeredDofMoFEMEntity FEDof(side_number_ptr,&**viit_rows);
	  rows_dofs.insert(FEDof);
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
	}
	if((*viit_rows)->part==pcomm->rank()) skip = false;
	parts[(*viit_rows)->part]++;
      }
      vector<int>::iterator pos = max_element(parts.begin(),parts.end());
      unsigned int max_part = std::distance(parts.begin(),pos);
      bool success = numered_finite_elements.modify(p.first,NumeredMoFEMFiniteElement_change_part(max_part));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      if(do_skip) if(skip) continue;
      //cols
      const void* tag_col_uids_data = miit3->tag_col_uids_data;
      const int tag_col_uids_size = miit3->tag_col_uids_size;
      ierr = get_MoFEMFiniteElement_dof_uid_view(p_miit->numered_dofs_cols,cols_view,Interface::UNION,tag_col_uids_data,tag_col_uids_size); CHKERRQ(ierr);
      NumeredDofMoFEMEntity_multiIndex_uid_view::iterator viit_cols = cols_view.begin();
      for(;viit_cols!=cols_view.end();viit_cols++) {
	try {
	  SideNumber *side_number_ptr = p.first->get_side_number_ptr(moab,(*viit_cols)->get_ent());
	  FENumeredDofMoFEMEntity_multiIndex &cols_dofs = const_cast<FENumeredDofMoFEMEntity_multiIndex&>(p.first->cols_dofs);
	  FENumeredDofMoFEMEntity FEDof(side_number_ptr,&**viit_cols);
	  cols_dofs.insert(FEDof);
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
	}
      }
      if(verb>1) {
	ostringstream ss;
	ss << *p_miit << endl;
	ss << problem_MoFEMFiniteElement << endl;
	typedef FENumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type FENumeredDofMoFEMEntity_multiIndex_by_Unique_mi_tag;
	FENumeredDofMoFEMEntity_multiIndex_by_Unique_mi_tag::iterator miit = problem_MoFEMFiniteElement.rows_dofs.get<Unique_mi_tag>().begin();
	for(;miit!=problem_MoFEMFiniteElement.rows_dofs.get<Unique_mi_tag>().end();miit++) ss << "rows: " << *miit << endl;
	miit = problem_MoFEMFiniteElement.cols_dofs.get<Unique_mi_tag>().begin();
	for(;miit!=problem_MoFEMFiniteElement.cols_dofs.get<Unique_mi_tag>().end();miit++) ss << "cols: " << *miit << endl;
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    }
  }
  if(verb>0) {
    typedef NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_Part_mi_tag>::type NumeredMoFEMFiniteElement_multiIndex_by_part;
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator MoFEMFiniteElement_miit = numered_finite_elements.get<MoFEMFiniteElement_Part_mi_tag>().lower_bound(pcomm->rank());
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator hi_MoMoFEMFiniteElement_miitFEMFE_miit = numered_finite_elements.get<MoFEMFiniteElement_Part_mi_tag>().upper_bound(pcomm->rank());
    int count = std::distance(MoFEMFiniteElement_miit,hi_MoMoFEMFiniteElement_miitFEMFE_miit);
    ostringstream ss;
    ss << *p_miit;
    ss << " Nb. elems " << count << " on proc " << pcomm->rank() << endl;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  *build_MoFEM |= 1<<5;  
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit) {
  PetscFunctionBegin; 
  try {
    Range ents3d;
    rval = moab.get_entities_by_type(meshset,MBTET,ents3d,false); CHKERR_PETSC(rval);
    Range ents;
    rval = moab.get_adjacencies(ents3d,2,true,ents,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(ents3d,1,true,ents,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.get_entities_by_type(meshset,MBPRISM,ents3d,false); CHKERR_PETSC(rval);
    Range::iterator tit = ents3d.begin();
    for(;tit!=ents3d.end();tit++) {
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*tit));
      if(debug > 0) {
	ierr = test_moab(moab,*tit); CHKERRQ(ierr);
      }
      if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
        bool success = refined_mofem_entities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }
      pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      switch (p_ent.first->get_ent_type()) {
        case MBTET: 
	 p_MoFEMFiniteElement = refined_mofem_elements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_ent.first)));	
	  assert(p_MoFEMFiniteElement.first->get_BitRefEdges_ulong()!=-1);
	 break;
	case MBPRISM:
	  p_MoFEMFiniteElement = refined_mofem_elements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
	  assert(p_MoFEMFiniteElement.first->get_BitRefEdges_ulong()!=-1);
	  break;
        case MBENTITYSET:
	  p_MoFEMFiniteElement = refined_mofem_elements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_MESHSET(moab,&*p_ent.first)));
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      if(verbose>2) {
        ostringstream ss;
        ss << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement()) << endl;
        PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    }
    for(int dd = 0;dd<2;dd++) {
      rval = moab.get_entities_by_dimension(meshset,dd,ents); CHKERR_PETSC(rval);
      Range::iterator eit = ents.begin();
      for(;eit!=ents.end();eit++) {
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*eit));
        if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
	  bool success = refined_mofem_entities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
        }
        if(verbose>2) {
  	ostringstream ss;
  	ss << *(p_ent.first) << endl;
  	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
        }
      }
    }
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit) {
  PetscFunctionBegin;
  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,meshset));
  if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
    refined_mofem_entities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
  }
  ptrWrapperRefMoFEMElement pack_fe(new RefMoFEMElement_MESHSET(moab,&*p_ent.first));
  pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement = refined_mofem_elements.insert(pack_fe);
  if(verbose > 0) {
    ostringstream ss;
    ss << "add meshset as ref_ent " << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement()) << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_verices_in_the_middel_of_edges(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMEntity_multiIndex::index<Composite_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refined_mofem_entities.get<Composite_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) ref_parent_ents_view.insert(&*miit);
  Range edges;
  rval = moab.get_entities_by_type(meshset,MBEDGE,edges,recursive);  CHKERR_PETSC(rval);
  if(edges.empty()) {
    Range tets;
    rval = moab.get_entities_by_type(meshset,MBTET,tets,recursive); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(tets,1,true,edges,Interface::UNION); CHKERR_PETSC(rval);
    if(tets.empty()) {
      Range prisms;
      rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,recursive); CHKERR_PETSC(rval);
      for(Range::iterator pit = prisms.begin();pit!=prisms.end();pit++) {
        const EntityHandle* conn; 
        int num_nodes; 
        rval = moab.get_connectivity(*pit,conn,num_nodes,true);  CHKERR_PETSC(rval);
        assert(num_nodes==6);
        //
        Range edge;
        rval = moab.get_adjacencies(&conn[0],2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        edge.clear();
        rval = moab.get_adjacencies(&conn[1],2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        EntityHandle conn_edge2[] = { conn[2], conn[0] };
        edge.clear();
        rval = moab.get_adjacencies(conn_edge2,2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        //
        edge.clear();
        rval = moab.get_adjacencies(&conn[3],2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        edge.clear();
        rval = moab.get_adjacencies(&conn[4],2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        EntityHandle conn_edge8[] = { conn[5], conn[3] };
        edge.clear();
        rval = moab.get_adjacencies(conn_edge8,2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
      }
    }
  }
  // refine edges on the other side of the prism
  typedef AdjacencyMapForBasicMoFEMEntity_multiIndex::index<MoABEnt_mi_tag2>::type AdjacencyMapForBasicMoFEMEntity_by_adj;
  AdjacencyMapForBasicMoFEMEntity_by_adj &adjacencies_maps_for_prisms_by_adj = adjacencies_maps_for_prisms.get<MoABEnt_mi_tag2>();
  Range::iterator eit = edges.begin();
  for(;eit!=edges.end();eit++) {
    AdjacencyMapForBasicMoFEMEntity_by_adj::iterator adj_miit = adjacencies_maps_for_prisms_by_adj.find(*eit);
    if(adj_miit==adjacencies_maps_for_prisms_by_adj.end()) continue;
    EntityHandle prism = adj_miit->ent;
    RefMoFEMElement_multiIndex::iterator miit2 = refined_mofem_elements.get<MoABEnt_mi_tag>().find(prism);
    SideNumber_multiIndex &side_table = miit2->get_side_number_table();
    SideNumber_multiIndex::iterator siit = side_table.find(*eit);
    int side_number = siit->side_number;
    if(side_number==-1) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
    if(prism_adj_edges[side_number]==-1) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
    EntityHandle edge;
    rval = moab.side_element(prism,1,prism_adj_edges[side_number],edge); CHKERR_PETSC(rval);
    if(edge==no_handle) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
    edges.insert(edge);
  }
  if(verb > 0) {
    ostringstream ss;
    ss << "ref level " << bit << " nb. edges to refine " << edges.size() << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  eit = edges.begin();
  for(;eit!=edges.end();eit++) {
    RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator miit_view = ref_parent_ents_view.find(*eit);
    const EntityHandle* conn; 
    int num_nodes; 
    rval = moab.get_connectivity(*eit,conn,num_nodes,true);  CHKERR_PETSC(rval);
    assert(num_nodes==2);
    if(miit_view == ref_parent_ents_view.end()||ref_parent_ents_view.empty()) {
      double coords[num_nodes*3]; 
      moab.get_coords(conn,num_nodes,coords); 
      cblas_daxpy(3,1.,&coords[3],1,coords,1);
      cblas_dscal(3,0.5,coords,1);
      EntityHandle node;
      rval = moab.create_vertex(coords,node); CHKERR_PETSC(rval);
      rval = moab.tag_set_data(th_RefParentHandle,&node,1,&*eit); CHKERR_PETSC(rval);
      rval = moab.tag_set_data(th_RefBitLevel,&node,1,&bit); CHKERR_PETSC(rval);
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,node));
      if(!p_ent.second) SETERRQ(PETSC_COMM_SELF,1,"this entity is there");
      if(verbose>2) {
	ostringstream ss;
	ss << *(p_ent.first) << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    } else {
      const EntityHandle node = (*miit_view)->get_ref_ent();
      bool success = refined_mofem_entities.modify(refined_mofem_entities.get<MoABEnt_mi_tag>().find(node),RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface) {
  //FIXME: refinment is based on entity handlers, should work on global ids of nodes, this will allow parallelize agortihm in the future
  PetscFunctionBegin;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ents_by_ent;
  ref_ents_by_ent &ref_ents_ent = refined_mofem_entities.get<MoABEnt_mi_tag>();
  // find all verices which parent is edge
  typedef RefMoFEMEntity_multiIndex::index<Composite_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refined_mofem_entities.get<Composite_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) ref_parent_ents_view.insert(&*miit);
  typedef RefMoFEMElement_multiIndex::index<MoABEnt_mi_tag>::type ref_MoFEMFiniteElement_by_ent;
  ref_MoFEMFiniteElement_by_ent &ref_MoFEMFiniteElement = refined_mofem_elements.get<MoABEnt_mi_tag>();
  typedef RefMoFEMElement_multiIndex::index<Composite_mi_tag>::type ref_ent_by_composite;
  ref_ent_by_composite &by_composite = refined_mofem_elements.get<Composite_mi_tag>();
  // find oposite intrface nodes
  typedef AdjacencyMapForBasicMoFEMEntity_multiIndex::index<EntType_mi_tag>::type AdjPrism_by_type;
  AdjPrism_by_type::iterator face_prism_miit = adjacencies_maps_for_prisms.get<EntType_mi_tag>().lower_bound(MBTRI);
  AdjPrism_by_type::iterator hi_face_prism_miit = adjacencies_maps_for_prisms.get<EntType_mi_tag>().upper_bound(MBTRI);
  map<EntityHandle,EntityHandle> nodes_face_map_for_faces_adj_to_prism_forward;
  map<EntityHandle,EntityHandle> nodes_face_map_for_faces_adj_to_prism_backward;
  for(;face_prism_miit!=hi_face_prism_miit;face_prism_miit++) {
    EntityHandle prism = face_prism_miit->ent;
    EntityHandle f3;
    rval = moab.side_element(prism,2,3,f3); CHKERR_PETSC(rval);
    EntityHandle f4;
    rval = moab.side_element(prism,2,4,f4); CHKERR_PETSC(rval);
    if((f4 == face_prism_miit->get_adj())&&(f3 != face_prism_miit->get_adj())) continue;
    const EntityHandle* conn_face; 
    int num_nodes_face; 
    moab.get_connectivity(face_prism_miit->get_adj(),conn_face,num_nodes_face,true); 
    assert(num_nodes_face==3);
    const EntityHandle* conn_face_other_side; 
    moab.get_connectivity(f3,conn_face_other_side,num_nodes_face,true); 
    assert(num_nodes_face==3);
    for(int nn = 0;nn<3;nn++) {
      nodes_face_map_for_faces_adj_to_prism_forward[conn_face[nn]] = conn_face_other_side[nn];
      nodes_face_map_for_faces_adj_to_prism_backward[conn_face_other_side[nn]] = conn_face[nn];
    }
  }
  //
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,false); CHKERR_PETSC(rval);
  Range::iterator tit = tets.begin();
  for(;tit!=tets.end();tit++) {
    ref_MoFEMFiniteElement_by_ent::iterator miit2 = ref_MoFEMFiniteElement.find(*tit);
    if(miit2==ref_MoFEMFiniteElement.end()) SETERRQ(PETSC_COMM_SELF,1,"this MoFEMFiniteElement is not there");
    //connectivity
    const EntityHandle* conn; 
    int num_nodes; 
    moab.get_connectivity(*tit,conn,num_nodes,true); 
    assert(num_nodes==4);
    for(int nn = 0;nn<num_nodes;nn++) {
      bool success = refined_mofem_entities.modify(refined_mofem_entities.get<MoABEnt_mi_tag>().find(conn[nn]),RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
    }
    //get edges
    BitRefEdges parent_edges_bit(0);
    EntityHandle edge_new_nodes[6];
    fill(&edge_new_nodes[0],&edge_new_nodes[6],no_handle); 
    int split_edges[6];  
    fill(&split_edges[0],&split_edges[6],-1); 
    map<EntityHandle,const RefMoFEMEntity*> map_ref_nodes_by_edges;
    for(int ee = 0;ee<6;ee++) { 
      EntityHandle edge = no_handle;
      rval = moab.side_element(*tit,1,ee,edge);  CHKERR_PETSC(rval);
      RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator miit_view = ref_parent_ents_view.find(edge);
      if(miit_view != ref_parent_ents_view.end()) {
	if(((*miit_view)->get_BitRefLevel()&bit).any()) {
	  edge_new_nodes[ee] = (*miit_view)->get_ref_ent(); 
	  map_ref_nodes_by_edges[(*miit_view)->get_parent_ent()] = &**miit_view;
	  split_edges[parent_edges_bit.count()] = ee;
	  parent_edges_bit.set(ee,1);
	}
      }
    }
    // swap nodes forward
    EntityHandle _conn_[4];
    copy(&conn[0],&conn[4],&_conn_[0]);
    if(respect_interface) {
      for(int nn = 0;nn<4;nn++) {
	map<EntityHandle,EntityHandle>::iterator mmit = nodes_face_map_for_faces_adj_to_prism_forward.find(_conn_[nn]);
	if(mmit!=nodes_face_map_for_faces_adj_to_prism_forward.end()) {
	  _conn_[nn] = mmit->second;
	}
      } 
    }
    // build connectivity for rf tets
    ref_ents_by_ent::iterator tit_miit;
    EntityHandle new_tets_conns[8*4];
    fill(&new_tets_conns[0],&new_tets_conns[8*4],no_handle);
    int sub_type = -1,nb_new_tets = 0;
    switch (parent_edges_bit.count()) {
      case 0:
	tit_miit = ref_ents_ent.find(*tit);
	if(tit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	refined_mofem_entities.modify(tit_miit,RefMoFEMEntity_change_add_bit(bit));
	break;
      case 1:
	sub_type = 0;
	tet_type_1(_conn_,split_edges[0],edge_new_nodes[split_edges[0]],new_tets_conns);
	nb_new_tets = 2;
	break;
      case 2:
	sub_type = tet_type_2(_conn_,split_edges,edge_new_nodes,new_tets_conns);
	if(sub_type&(4|8|16)) {
	  nb_new_tets = 3;
	  break;
	} else if(sub_type == 1) {
	  nb_new_tets = 4;
	  break; 
	};
	assert(0);
	break;
      case 3:
	sub_type = tet_type_3(_conn_,split_edges,edge_new_nodes,new_tets_conns);
	if(sub_type <= 4 ) {
	  nb_new_tets = 4;
	  break;
	} else if(sub_type <= 7 ) {
	  nb_new_tets = 5;
	  break;
	}
	assert(0);
      case 4:
	sub_type = tet_type_4(_conn_,split_edges,edge_new_nodes,new_tets_conns);
	if(sub_type == 0) {
	  nb_new_tets = 5;
	  break;
	} else if(sub_type <= 7) {
	  nb_new_tets = 6;
	  break;
	}  
	assert(0);
      case 5:
	sub_type = tet_type_5(moab,_conn_,edge_new_nodes,new_tets_conns);
	nb_new_tets = 7;
	break;
      case 6:
	sub_type = 0;
	tet_type_6(moab,_conn_,edge_new_nodes,new_tets_conns);
	nb_new_tets = 8;
	break;
      default:
	assert(0);
    }
    // swap nodes backward 
    if(respect_interface) {
      for(int nn = 0;nn<8*4;nn++) {
	if(new_tets_conns[nn]==no_handle) continue;
	map<EntityHandle,EntityHandle>::iterator mmit = nodes_face_map_for_faces_adj_to_prism_backward.find(new_tets_conns[nn]);
	if(mmit!=nodes_face_map_for_faces_adj_to_prism_backward.end()) {
	  new_tets_conns[nn] = mmit->second;
	}
      }
    }
    // find that tets
    bitset<8> ref_tets_bit(0);
    ref_ent_by_composite::iterator miit_composite = by_composite.lower_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    ref_ent_by_composite::iterator hi_miit_composite = by_composite.upper_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    ref_ent_by_composite::iterator miit_composite2 = miit_composite;
    for(int tt = 0;miit_composite2!=hi_miit_composite;miit_composite2++,tt++) {
      //add this tet to this ref
      refined_mofem_entities.modify(refined_mofem_entities.find(miit_composite2->get_ref_ent()),RefMoFEMEntity_change_add_bit(bit));
      ref_tets_bit.set(tt,1);
      if(verbose>2) {
	ostringstream ss;
	ss << miit_composite2->get_RefMoFEMElement() << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    }
    if(miit_composite!=hi_miit_composite) {
      if(ref_tets_bit.count()!=(unsigned int)nb_new_tets) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    } else {
      //create tets
      EntityHandle ref_tets[8];
      for(int tt = 0;tt<nb_new_tets;tt++) {
	if(!ref_tets_bit.test(tt)) {
	  rval = moab.create_element(MBTET,&new_tets_conns[4*tt],4,ref_tets[tt]); CHKERR_PETSC(rval);
	  double coords[12];
	  ierr = moab.get_coords(&new_tets_conns[4*tt],4,coords); CHKERRQ(ierr);
	  double V = Shape_intVolumeMBTET(diffN_TET,coords); 
	  if(V<=0) {
	    ostringstream ss;
	    ss << "tit " << new_tets_conns[4*tt];
	    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	    assert(V>0); 
	  }
	  int ref_type[] = { parent_edges_bit.count(),sub_type }; 
	  rval = moab.tag_set_data(th_RefType,&ref_tets[tt],1,ref_type); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefParentHandle,&ref_tets[tt],1,&*tit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitLevel,&ref_tets[tt],1,&bit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitEdge,&ref_tets[tt],1,&parent_edges_bit); CHKERR_PETSC(rval);
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,ref_tets[tt]));
	  pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement = refined_mofem_elements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_ent.first)));
	  ref_tets_bit.set(tt);
	  if(verbose>2) {
	    ostringstream ss;
	    ss << "add tet: " << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement()) << endl;
	    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }
	}
      }
      // find parents for new edges and faces
      Range tit_edges,tit_faces;
      rval = moab.get_adjacencies(&*tit,1,1,false,tit_edges); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(&*tit,1,2,false,tit_faces); CHKERR_PETSC(rval);
      Range edges_nodes[6],faces_nodes[4];
      // for edges - add ref nodes
      Range::iterator eit=tit_edges.begin();
      for(int ee = 0;eit!=tit_edges.end();eit++,ee++) {
	rval = moab.get_connectivity(&*eit,1,edges_nodes[ee],true); CHKERR_PETSC(rval);
	map<EntityHandle,const RefMoFEMEntity*>::iterator map_miit = map_ref_nodes_by_edges.find(*eit);
	if(map_miit!=map_ref_nodes_by_edges.end()) {
	  edges_nodes[ee].insert(map_miit->second->get_ref_ent());
	}
      }
      // for faces - add ref nodes 
      Range::iterator fit=tit_faces.begin();
      for(int ff = 0;fit!=tit_faces.end();fit++,ff++) {
	rval = moab.get_connectivity(&*fit,1,faces_nodes[ff],true); CHKERR_PETSC(rval);
	Range fit_edges;
	rval = moab.get_adjacencies(&*fit,1,1,false,fit_edges); CHKERR_PETSC(rval);
	for(Range::iterator eit2 =  fit_edges.begin();eit2 != fit_edges.end();eit2++) {
	  map<EntityHandle,const RefMoFEMEntity*>::iterator map_miit = map_ref_nodes_by_edges.find(*eit2);
	  if(map_miit!=map_ref_nodes_by_edges.end()) {
	    faces_nodes[ff].insert(map_miit->second->get_ref_ent());
	  }
	}
      }
      // add ref nodes to tet
      Range tet_nodes;
      rval = moab.get_connectivity(&*tit,1,tet_nodes,true); CHKERR_PETSC(rval);
      for(map<EntityHandle,const RefMoFEMEntity*>::iterator map_miit = map_ref_nodes_by_edges.begin();
	map_miit != map_ref_nodes_by_edges.end();map_miit++) {
	tet_nodes.insert(map_miit->second->get_ref_ent());
      }
      Range ref_edges;
      rval = moab.get_adjacencies(ref_tets,nb_new_tets,1,true,ref_edges,Interface::UNION); CHKERR_PETSC(rval);
      // check for all ref edge
      for(Range::iterator reit = ref_edges.begin();reit!=ref_edges.end();reit++) {
	Range ref_edges_nodes;
	rval = moab.get_connectivity(&*reit,1,ref_edges_nodes,true); CHKERR_PETSC(rval);
	// check if ref edge is in coarse edge
	int ee = 0;
	for(;ee<6;ee++) {
	  // two nodes are common (node[0],node[1],ref_node (if exist))
	  if(intersect(edges_nodes[ee],ref_edges_nodes).size()==2) {
	    EntityHandle edge = tit_edges[ee];
	    rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&edge); CHKERR_PETSC(rval);
	    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*reit));
	    bool success = refined_mofem_entities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
	    if(p_ent.second) {
	      if(verbose>2) {
		ostringstream ss;
		ss << "edge parent: " << *(p_ent.first) << endl;
		PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }}
	    break;
	  }  
	}
	if(ee<6) continue;
	// check if ref edge is in coarse face
	int ff = 0;
	for(;ff<4;ff++) {
	  // two nodes are common (node[0],node[1],ref_node (if exist))
	  if(intersect(faces_nodes[ff],ref_edges_nodes).size()==2) {
	    EntityHandle face = tit_faces[ff];
	    rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&face); CHKERR_PETSC(rval);
	    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*reit));
	    bool success = refined_mofem_entities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
	    if(p_ent.second) {
	      if(verbose>2) {
		ostringstream ss;
		ss << "face parent: " << *(p_ent.first) << endl;
		PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }}
	    break;
	  }
	}
	if(ff<4) continue;
	// check if ref edge is in coarse tetrahedral
	if(intersect(tet_nodes,ref_edges_nodes).size()==2) {
	  rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&*tit); CHKERR_PETSC(rval);
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*reit));
	  bool success = refined_mofem_entities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
	  if(p_ent.second) {
	    if(verbose>2) {
	      ostringstream ss;
	      ss << "tet parent: " << *(p_ent.first) << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }}
	  continue;
	}
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      Range ref_faces;
      rval = moab.get_adjacencies(ref_tets,nb_new_tets,2,true,ref_faces,Interface::UNION); CHKERR_PETSC(rval);
      // check for all ref edges
      for(Range::iterator rfit = ref_faces.begin();rfit!=ref_faces.end();rfit++) {
	Range ref_faces_nodes;
	rval = moab.get_connectivity(&*rfit,1,ref_faces_nodes,true); CHKERR_PETSC(rval);
	// check if ref face is in coarse face
	int ff = 0;
	for(;ff<4;ff++) {
	  if(intersect(faces_nodes[ff],ref_faces_nodes).size()==3) {
	    EntityHandle face = tit_faces[ff];
	    rval = moab.tag_set_data(th_RefParentHandle,&*rfit,1,&face); CHKERR_PETSC(rval);
	    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*rfit));
	    bool success = refined_mofem_entities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
	    if(p_ent.second) {
	      if(verbose>2) {
		ostringstream ss;
		ss << "face parent: " << *(p_ent.first) << endl;
		PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }}
	    break;
	  }
	}
	if(ff<4) continue;
	// check if ref face is in coarse tetrahedral
	if(intersect(tet_nodes,ref_faces_nodes).size()==3) {
	  rval = moab.tag_set_data(th_RefParentHandle,&*rfit,1,&*tit); CHKERR_PETSC(rval);
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,*rfit));
	  bool success = refined_mofem_entities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
	  if(p_ent.second) {
	    if(verbose>2) {
	      ostringstream ss;
	      ss << "tet parent: " << *(p_ent.first) << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }}
	  continue;
	}
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  //FIXME: refinment is based on entity handlers, should work on global ids of nodes, this will allow parallelize agortihm in the future
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ENTs_by_ent;
  typedef RefMoFEMElement_multiIndex::index<Composite_mi_tag>::type ref_fe_by_composite;
  ref_fe_by_composite &ref_fe_by_comp = refined_mofem_elements.get<Composite_mi_tag>();
  // find all verices which parent is edge
  typedef RefMoFEMEntity_multiIndex::index<Composite_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents_by_comp = refined_mofem_entities.get<Composite_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents_by_comp.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents_by_comp.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) ref_parent_ents_view.insert(&*miit);
  //
  Range prisms;
  rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,false); CHKERR_PETSC(rval);
  Range::iterator pit = prisms.begin();
  for(;pit!=prisms.end();pit++) {
    ref_ENTs_by_ent::iterator miit_prism = refined_mofem_entities.get<MoABEnt_mi_tag>().find(*pit);   
    if(miit_prism==refined_mofem_entities.end()) SETERRQ(PETSC_COMM_SELF,1,"this prism is not in ref database");
    if(verb>3) {
      ostringstream ss;
      ss << "ref prism " << *miit << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    //prism connectivity
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*pit,conn,num_nodes,true); CHKERR_PETSC(rval);
    assert(num_nodes==6);
    // edges connectivity
    EntityHandle edges[6];
    for(int ee = 0;ee<3; ee++) {
      rval = moab.side_element(*pit,1,ee,edges[ee]); CHKERR_PETSC(rval);
    }
    for(int ee = 6;ee<9; ee++) {
      rval = moab.side_element(*pit,1,ee,edges[ee-3]); CHKERR_PETSC(rval);
    }
    // detetct split edges
    BitRefEdges split_edges(0);
    EntityHandle edge_nodes[6];
    fill(&edge_nodes[0],&edge_nodes[6],no_handle);
    for(int ee = 0;ee<6;ee++) {
      RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator miit_view = ref_parent_ents_view.find(edges[ee]);
      if(miit_view != ref_parent_ents_view.end()) {
	if(((*miit_view)->get_BitRefLevel()&bit).any()) {
	  edge_nodes[ee] = (*miit_view)->get_ref_ent(); 
	  split_edges.set(ee);
	}
      }
    }
    if(split_edges.count()==0) {
      refined_mofem_entities.modify(miit_prism,RefMoFEMEntity_change_add_bit(bit));
      if(verb>6) PetscPrintf(PETSC_COMM_WORLD,"no refinment");
      continue;
    } 
    //check consitency
    if(verb>3) {
      ostringstream ss;
      ss << "prism split edges " << split_edges << " count " << split_edges.count() << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    // prism ref
    EntityHandle new_prism_conn[4*6];
    fill(&new_prism_conn[0],&new_prism_conn[4*6],no_handle);
    int nb_new_prisms = 0;
    switch (split_edges.count()) {
      case 0:
	break;
      case 2:
	ierr = prism_type_1(conn,split_edges,edge_nodes,new_prism_conn); CHKERRQ(ierr);
	nb_new_prisms = 2;
	break;
      case 4:
	ierr = prism_type_2(conn,split_edges,edge_nodes,new_prism_conn); CHKERRQ(ierr);
	nb_new_prisms = 3;
	break;
      case 6:
	ierr = prism_type_3(conn,split_edges,edge_nodes,new_prism_conn); CHKERRQ(ierr);
	nb_new_prisms = 4;
	break;
      default:
	ostringstream ss;
	ss << split_edges << " : [ " 
	  << conn[0] << " "
	  << conn[1] << " "
	  << conn[2] << " "
	  << conn[3] << " "
	  << conn[4] << " "
	  << conn[5] << " ]";
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    // find that prism
    bitset<4> ref_prism_bit(0);
    ref_fe_by_composite::iterator miit_composite = ref_fe_by_comp.lower_bound(boost::make_tuple(*pit,split_edges.to_ulong()));
    ref_fe_by_composite::iterator hi_miit_composite = ref_fe_by_comp.upper_bound(boost::make_tuple(*pit,split_edges.to_ulong()));
    ref_fe_by_composite::iterator miit_composite2 = miit_composite;
    for(int pp = 0;miit_composite2!=hi_miit_composite;miit_composite2++,pp++) {
      //add this tet to this ref
      refined_mofem_entities.modify(refined_mofem_entities.find(miit_composite2->get_ref_ent()),RefMoFEMEntity_change_add_bit(bit));
      ref_prism_bit.set(pp,1);
      if(verb>2) {
	ostringstream ss;
	ss << "is refined " << *(miit_composite2->get_RefMoFEMElement()) << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    }
    if(miit_composite!=hi_miit_composite) {
      if(ref_prism_bit.count()!=(unsigned int)nb_new_prisms) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    } else {
      EntityHandle ref_prisms[4];
      // create prism
      for(int pp = 0;pp<nb_new_prisms;pp++) {
	if(verb>3) {
	  ostringstream ss;
	  ss << "ref prism " << ref_prism_bit << endl;
	  PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	}
	if(!ref_prism_bit.test(pp)) {
	  rval = moab.create_element(MBPRISM,&new_prism_conn[6*pp],6,ref_prisms[pp]); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefParentHandle,&ref_prisms[pp],1,&*pit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitLevel,&ref_prisms[pp],1,&bit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitEdge,&ref_prisms[pp],1,&split_edges); CHKERR_PETSC(rval);
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,ref_prisms[pp]));
	  pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
	  try {
	    p_MoFEMFiniteElement = refined_mofem_elements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
	  } catch (const char* msg) {
	    SETERRQ(PETSC_COMM_SELF,1,msg);
	  }
	  ref_prism_bit.set(pp);
	  ierr = add_prism_to_adjacencies_maps_for_prisms(ref_prisms[pp]); CHKERRQ(ierr);
	  if(verb>2) {
	    ostringstream ss;
	    ss << "add prism: " << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement()) << endl;
	    if(verb>7) {
	      for(int nn = 0;nn<6;nn++) {
		ss << new_prism_conn[nn] << " ";
	      }
	      ss << endl;
	    }
	    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ENTs_by_ent;
  ref_ENTs_by_ent::iterator miit = refined_mofem_entities.find(meshset);
  if(miit==refined_mofem_entities.end()) SETERRQ(PETSC_COMM_SELF,1,"this meshset is not in ref database");
  ierr = refine_get_childern(meshset,bit,meshset,MBEDGE,recursive,verb); CHKERRQ(ierr);
  ierr = refine_get_childern(meshset,bit,meshset,MBTRI,recursive,verb); CHKERRQ(ierr);
  ierr = refine_get_childern(meshset,bit,meshset,MBTET,recursive,verb); CHKERRQ(ierr);
  refined_mofem_entities.modify(miit,RefMoFEMEntity_change_add_bit(bit));
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::refine_get_finite_elements(const BitRefLevel &bit,const EntityHandle meshset) {
  PetscFunctionBegin;
  RefMoFEMElement_multiIndex::iterator miit = refined_mofem_elements.begin();
  for(;miit!=refined_mofem_elements.end();miit++) {
    BitRefLevel bit2 = miit->get_BitRefLevel(); 
    if((bit2&bit).any()) {
      switch (miit->get_ent_type()) {
	case MBTET:
	case MBPRISM:
	break;
	case MBENTITYSET:
	continue;
	default:
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      EntityHandle ent = miit->get_ref_ent();
      int type = miit->get_ent_type();
      rval = moab.tag_set_data(th_ElemType,&ent,1,&type); CHKERR_PETSC(rval);
      rval = moab.add_entities(meshset,&ent,1); CHKERR_PETSC(rval);
    }
  }	
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::refine_get_ents(const BitRefLevel &bit,const EntityHandle meshset) {
  PetscFunctionBegin;
  RefMoFEMEntity_multiIndex::iterator miit = refined_mofem_entities.begin();
  for(;miit!=refined_mofem_entities.end();miit++) {
    BitRefLevel bit2 = miit->get_BitRefLevel(); 
    if((bit2&bit)==bit) {
      EntityHandle ent = miit->get_ref_ent();
      rval = moab.add_entities(meshset,&ent,1); CHKERR_PETSC(rval);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::refine_get_childern(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refined_mofem_entities.get<Composite_mi_tag2>();
  Range ents;
  rval = moab.get_entities_by_handle(parent,ents,recursive);  CHKERR_PETSC(rval);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    if(verb>2) {
      ostringstream ss;
      ss << "ent " << *eit << endl;;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(*eit,child_type));
    ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(*eit,child_type));
    for(;miit!=hi_miit;miit++) {
      if(verb>2) {
	ostringstream ss;
	ss << "any bit " << *miit << endl;;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      if((miit->get_BitRefLevel()&child_bit).any()) {
	EntityHandle ref_ent = miit->get_ref_ent();
	rval = moab.add_entities(child,&ref_ent,1); CHKERR_PETSC(rval);
	if(verb>1) {
	  ostringstream ss;
	  ss << "good bit " << *miit << endl;;
	  PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::problem_get_FE(const string &problem_name,const string &fe_name,const EntityHandle meshset) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(problem_name);
  NumeredMoFEMFiniteElement_multiIndex &numered_finite_elements = const_cast<NumeredMoFEMFiniteElement_multiIndex&>(p_miit->numered_finite_elements);
  NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator miit = numered_finite_elements.get<MoFEMFiniteElement_name_mi_tag>().lower_bound(fe_name);
  for(;miit!=numered_finite_elements.get<MoFEMFiniteElement_name_mi_tag>().upper_bound(fe_name);miit++) {
    EntityHandle ent = miit->get_ent();
    rval = moab.add_entities(meshset,&ent,1); CHKERR_PETSC(rval);
    int part = miit->get_part();
    rval = moab.tag_set_data(th_Part,&ent,1,&part); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  moabBaseMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    ierr = miit->get_Cubit_msId_entities_by_dimension(moab,dimension,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  moabBaseMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    ierr = miit->get_Cubit_msId_entities_by_dimension(moab,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_Cubit_msId_entities_by_dimension(msId,Cubit_BC_bitset(CubitBCType),dimension,entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType,
  Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_Cubit_msId_entities_by_dimension(msId,Cubit_BC_bitset(CubitBCType),entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset) {
  PetscFunctionBegin;
  moabBaseMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    meshset = miit->meshset;
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_CubitBCType_meshsets(const unsigned int CubitBCType,Range &meshsets) {
  PetscFunctionBegin;
  moabBaseMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<CubitMeshSets_mi_tag>().lower_bound(CubitBCType);
  moabBaseMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator 
    hi_miit = cubit_meshsets.get<CubitMeshSets_mi_tag>().upper_bound(CubitBCType);
  for(;miit!=hi_miit;miit++) {
    meshsets.insert(miit->meshset);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::VecCreateGhost(const string &name,RowColData rc,Vec *V) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_local_idx;
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(name);
  DofIdx nb_dofs,nb_local_dofs,nb_ghost_dofs;
  dofs_by_local_idx *dofs;
  switch (rc) {
    case Row:
      nb_dofs = p_miit->get_nb_dofs_row();
      nb_local_dofs = p_miit->get_nb_local_dofs_row();
      nb_ghost_dofs = p_miit->get_nb_ghost_dofs_row();
      dofs = const_cast<dofs_by_local_idx*>(&p_miit->numered_dofs_rows.get<PetscLocalIdx_mi_tag>());
      break;
    case Col:
      nb_dofs = p_miit->get_nb_dofs_col();
      nb_local_dofs = p_miit->get_nb_local_dofs_col();
      nb_ghost_dofs = p_miit->get_nb_ghost_dofs_col();
      dofs = const_cast<dofs_by_local_idx*>(&p_miit->numered_dofs_cols.get<PetscLocalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  dofs_by_local_idx::iterator miit = dofs->lower_bound(nb_local_dofs);
  dofs_by_local_idx::iterator hi_miit = dofs->upper_bound(nb_local_dofs+nb_ghost_dofs);
  int count = std::distance(miit,hi_miit);
  if(count != nb_ghost_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
  vector<DofIdx> ghost_idx(count);
  vector<DofIdx>::iterator vit = ghost_idx.begin();
  for(;miit!=hi_miit;miit++,vit++) *vit = miit->petsc_gloabl_dof_idx;
  ierr = ::VecCreateGhost(PETSC_COMM_WORLD,nb_local_dofs,nb_dofs,nb_ghost_dofs,&ghost_idx[0],V); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::MatCreateMPIAIJWithArrays(const string &name,Mat *Aij,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = partition_create_Mat<Part_mi_tag>(name,NULL,Aij,false,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::set_local_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_local_idx;
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(name);
  dofs_by_local_idx *dofs;
  DofIdx nb_local_dofs,nb_ghost_dofs;
  switch (rc) {
    case Row:
      nb_local_dofs = p_miit->get_nb_local_dofs_row();
      nb_ghost_dofs = p_miit->get_nb_ghost_dofs_row();
      dofs = const_cast<dofs_by_local_idx*>(&p_miit->numered_dofs_rows.get<PetscLocalIdx_mi_tag>());
      break;
    case Col:
      nb_local_dofs = p_miit->get_nb_local_dofs_col();
      nb_ghost_dofs = p_miit->get_nb_ghost_dofs_col();
      dofs = const_cast<dofs_by_local_idx*>(&p_miit->numered_dofs_cols.get<PetscLocalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  Vec Vlocal;
  ierr = VecGhostGetLocalForm(V,&Vlocal); CHKERRQ(ierr);
  PetscInt size;
  ierr = VecGetLocalSize(V,&size); CHKERRQ(ierr);
  if(size!=nb_local_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency: check ghost vector, problem with nb. of local nodes"); 
  ierr = VecGetLocalSize(Vlocal,&size); CHKERRQ(ierr);
  if(size!=nb_local_dofs+nb_ghost_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency: check ghost vector, problem with nb. of ghost nodes"); 
  dofs_by_local_idx::iterator miit = dofs->lower_bound(0);
  dofs_by_local_idx::iterator hi_miit = dofs->upper_bound(nb_local_dofs+nb_ghost_dofs);
  PetscScalar *array;
  VecGetArray(Vlocal,&array);
  DofIdx ii = 0;
  switch (scatter_mode) {
    case SCATTER_FORWARD:
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) array[ii] = miit->get_FieldData();
	  break;
	case ADD_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) array[ii] += miit->get_FieldData();	  
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
    break;
    case SCATTER_REVERSE:
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) miit->get_FieldData() = array[ii];
	  break;
	case ADD_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) miit->get_FieldData() += array[ii];	  
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
    break;
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  VecRestoreArray(Vlocal,&array);
  VecDestroy(&Vlocal);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::set_global_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type dofs_by_global_idx;
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(name);
  dofs_by_global_idx *dofs;
  DofIdx nb_dofs;
  switch (rc) {
    case Row:
      nb_dofs = p_miit->get_nb_dofs_row();
      dofs = const_cast<dofs_by_global_idx*>(&p_miit->numered_dofs_rows.get<PetscGlobalIdx_mi_tag>());
      break;
    case Col:
      nb_dofs = p_miit->get_nb_dofs_col();
      dofs = const_cast<dofs_by_global_idx*>(&p_miit->numered_dofs_cols.get<PetscGlobalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  dofs_by_global_idx::iterator miit = dofs->lower_bound(0);
  dofs_by_global_idx::iterator hi_miit = dofs->upper_bound(nb_dofs);
  DofIdx ii = 0;
  switch (scatter_mode) {
    case SCATTER_REVERSE: {
      VecScatter ctx;
      Vec V_glob;
      ierr = VecScatterCreateToAll(V,&ctx,&V_glob); CHKERRQ(ierr);
      ierr = VecScatterBegin(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      PetscInt size;
      ierr = VecGetSize(V_glob,&size); CHKERRQ(ierr);
      if(size!=nb_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency: nb. of dofs and decalared nb. dofs in database");
      if(size!=distance(miit,hi_miit)) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency: nb. of dofs and decalared nb. dofs in database");
      PetscScalar *array;
      ierr = VecGetArray(V_glob,&array); CHKERRQ(ierr);
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) miit->get_FieldData() = array[ii];
	  break;
	case ADD_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) miit->get_FieldData() += array[ii];	  
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      ierr = VecRestoreArray(V_glob,&array); CHKERRQ(ierr);
      ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
      ierr = VecDestroy(&V_glob); CHKERRQ(ierr);
    break;
    }
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::set_other_global_VecCreateGhost(
  const string &name,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,
  int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type dofs_by_name;
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(name);
  if(p_miit==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found",name.c_str());
  dofs_by_name *dofs;
  DofIdx nb_dofs;
  switch (rc) {
    case Row:
      nb_dofs = p_miit->get_nb_dofs_row();
      dofs = const_cast<dofs_by_name*>(&p_miit->numered_dofs_rows.get<FieldName_mi_tag>());
      break;
    case Col:
      nb_dofs = p_miit->get_nb_dofs_col();
      dofs = const_cast<dofs_by_name*>(&p_miit->numered_dofs_cols.get<FieldName_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit = moabfields.get<FieldName_mi_tag>().find(cpy_field_name);
  if(cpy_fit==moabfields.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"cpy field < %s > not found, (top tip: check spelling)",cpy_field_name.c_str());
  }
  dofs_by_name::iterator miit = dofs->lower_bound(field_name);
  if(miit==dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"cpy field < %s > not found, (top tip: check spelling)",field_name.c_str());
  }
  dofs_by_name::iterator hi_miit = dofs->upper_bound(field_name);
  if(miit->get_space() != cpy_fit->get_space()) {
    SETERRQ(PETSC_COMM_SELF,1,"fiedls has to have same space");
  }
  if(miit->get_max_rank() != cpy_fit->get_max_rank()) {
    SETERRQ(PETSC_COMM_SELF,1,"fiedls has to have same rank");
  }
  switch (scatter_mode) {
    case SCATTER_REVERSE: {
      Vec V_glob;
      VecScatter ctx;
      ierr = VecScatterCreateToAll(V,&ctx,&V_glob); CHKERRQ(ierr);
      ierr = VecScatterBegin(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      PetscInt size;
      ierr = VecGetSize(V_glob,&size); CHKERRQ(ierr);
      if(size!=nb_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency: nb. of dofs and decalared nb. dofs in database");
      PetscScalar *array;
      VecGetArray(V_glob,&array);
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++) {
	    if(miit->get_petsc_gloabl_dof_idx()>=size) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsitency: nb. of dofs and decalared nb. dofs in database");
	    }
	    DofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator diiiit;
	    diiiit = dofs_moabfield.get<Composite_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	    if(diiiit==dofs_moabfield.get<Composite_mi_tag>().end()) {
	      EntityHandle ent = miit->get_ent();
	      rval = moab.add_entities(cpy_fit->get_meshset(),&ent,1); CHKERR_PETSC(rval);
	      //create field moabent
	      ApproximationOrder order = miit->get_max_order();
	      pair<MoFEMEntity_multiIndex::iterator,bool> p_e_miit;
	      try {
		MoFEMEntity moabent(moab,cpy_fit->get_MoFEMField_ptr(),miit->get_RefMoFEMEntity_ptr());
		p_e_miit = ents_moabfield.insert(moabent);
	      } catch (const std::exception& ex) {
		ostringstream ss;
		ss << ex.what() << endl;
		SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	      }
	      if(p_e_miit.first->get_max_order()<order) {
		bool success = ents_moabfield.modify(p_e_miit.first,MoFEMEntity_change_order(moab,order));
		if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
	      }
	      //create field moabdof
	      DofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator hi_diit,diit;
	      diit = dofs_moabfield.get<Composite_mi_tag2>().lower_bound(boost::make_tuple(field_name,miit->get_ent()));
	      hi_diit = dofs_moabfield.get<Composite_mi_tag2>().upper_bound(boost::make_tuple(field_name,miit->get_ent()));
	      for(;diit!=hi_diit;diit++) {
		DofMoFEMEntity mdof(&*(p_e_miit.first),diit->get_dof_order(),diit->get_dof_rank(),diit->get_EntDofIdx());
		pair<DofMoFEMEntity_multiIndex::iterator,bool> cpy_p_diit;
		cpy_p_diit = dofs_moabfield.insert(mdof);
		if(cpy_p_diit.second) {
		  bool success = dofs_moabfield.modify(cpy_p_diit.first,DofMoFEMEntity_active_change(true));
		  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
		}
	      }
	      diiiit = dofs_moabfield.get<Composite_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	      if(diiiit==dofs_moabfield.get<Composite_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
	    }
	    diiiit->get_FieldData() = array[miit->get_petsc_gloabl_dof_idx()];
	    if(verb > 1) {
	      ostringstream ss;
	      ss << *diiiit << "set " << array[miit->get_petsc_gloabl_dof_idx()] << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }
	  }
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      ierr = VecRestoreArray(V_glob,&array); CHKERRQ(ierr);
      ierr = VecDestroy(&V_glob); CHKERRQ(ierr);
      ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    }
    break;
    case SCATTER_FORWARD:
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++) {
	    DofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator diiiit;
	    diiiit = dofs_moabfield.get<Composite_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	    if(diiiit==dofs_moabfield.get<Composite_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no data to fill the vector");
	    ierr = VecSetValue(V,miit->get_petsc_gloabl_dof_idx(),diiiit->get_FieldData(),INSERT_VALUES); CHKERRQ(ierr);
	  }
	  ierr = VecAssemblyBegin(V); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(V); CHKERRQ(ierr);
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
    break;  
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_msId_3dENTS_sides(const int msId,const Cubit_BC_bitset CubitBCType,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  moabBaseMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    ierr = moabField_Core::get_msId_3dENTS_sides(miit->meshset,recursive,verb); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_msId_3dENTS_sides(const EntityHandle SideSet,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Skinner skin(&moab);
  Range triangles;
  rval = moab.get_entities_by_type(SideSet,MBTRI,triangles,recursive);  CHKERR_PETSC(rval);
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"Nb. of triangles in set %u\n",triangles.size());
  }
  Range nodes; // nodes from triangles
  rval = moab.get_connectivity(triangles,nodes,true); CHKERR_PETSC(rval);
  Range edges; // edges from truangles
  rval = moab.get_adjacencies(triangles,1,true,edges,Interface::UNION); CHKERR_PETSC(rval);
  Range ents3d; // 3d ents form nodes
  rval = moab.get_adjacencies(nodes,3,true,ents3d,Interface::UNION); CHKERR_PETSC(rval);
  //
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(ents3d,false,skin_faces); CHKERR(rval);
  Range skin_edges; // skin edges from triangles
  rval = skin.find_skin(triangles,false,skin_edges); CHKERR(rval);
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"skin_edges %u\n",skin_edges.size());
  //
  Range skin_faces_edges; // edges from skin faces of 3d ents
  rval = moab.get_adjacencies(skin_faces,1,true,skin_faces_edges,Interface::UNION); CHKERR_PETSC(rval);
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"skin_faces_edges %u\n",skin_faces_edges.size());
  // 
  // skin edges are internal edge <- skin_faces_edges contains edges which are on the body boundary <- that is the trick
  skin_edges = subtract(skin_edges,skin_faces_edges); // from skin edges subtract edges from skin faces of 3d ents
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"subtract skin_edges %u\n",skin_edges.size());
  Range skin_nodes;
  rval = moab.get_connectivity(skin_edges,skin_nodes,true); CHKERR_PETSC(rval);
  nodes = subtract(nodes,skin_nodes); // nodes adjacent to all splitted face edges except those on internal edge
  ents3d.clear();
  // ents3 that are adjacent to nodes on splitted faces but not those which are on the nodes on internal edge
  rval = moab.get_adjacencies(nodes,3,true,ents3d,Interface::UNION); CHKERR_PETSC(rval);
  Range side_ents3d;
  unsigned int nb_side_ents3d = side_ents3d.size();
  side_ents3d.insert(*ents3d.begin());
  Range adj_tris,adj_ents3d;
  do {
    nb_side_ents3d = side_ents3d.size();
    if(verb>2) PetscPrintf(PETSC_COMM_WORLD,"nb_side_ents3d %u\n",nb_side_ents3d);
    rval = moab.get_adjacencies(side_ents3d,2,true,adj_tris,Interface::UNION); CHKERR_PETSC(rval);
    if(verb>2) PetscPrintf(PETSC_COMM_WORLD,"adj_tris %u\n",adj_tris.size());
    adj_tris = subtract(adj_tris,triangles);
    if(verb>2) PetscPrintf(PETSC_COMM_WORLD,"adj_tris %u\n",adj_tris.size());
    rval = moab.get_adjacencies(adj_tris,3,true,adj_ents3d,Interface::UNION); CHKERR_PETSC(rval);
    if(verb>2) PetscPrintf(PETSC_COMM_WORLD,"adj_ents3d %u\n",adj_ents3d.size());
    adj_ents3d = intersect(adj_ents3d,ents3d);
    if(verb>2) PetscPrintf(PETSC_COMM_WORLD,"adj_ents3d %u\n",adj_ents3d.size());
    side_ents3d.insert(adj_ents3d.begin(),adj_ents3d.end());
  } while (nb_side_ents3d != side_ents3d.size());
  Range other_side = subtract(ents3d,side_ents3d);
  vector<EntityHandle> children;
  rval = moab.get_child_meshsets(SideSet,children);  CHKERR_PETSC(rval);
  if(children.empty()) {
    children.resize(3);
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,children[0]); CHKERR_PETSC(rval);
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,children[1]); CHKERR_PETSC(rval);
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,children[2]); CHKERR_PETSC(rval);
    rval = moab.add_child_meshset(SideSet,children[0]); CHKERR_PETSC(rval);
    rval = moab.add_child_meshset(SideSet,children[1]); CHKERR_PETSC(rval);
    rval = moab.add_child_meshset(children[0],children[2]); CHKERR_PETSC(rval);
    rval = moab.add_child_meshset(children[1],children[2]); CHKERR_PETSC(rval);

  } else { assert(children.size()==2); }
  EntityHandle &child_side = children[0];
  EntityHandle &child_other_side = children[1];
  EntityHandle &child_nodes_and_skin_edges = children[2];
  rval = moab.add_entities(child_side,side_ents3d); CHKERR_PETSC(rval);
  rval = moab.add_entities(child_other_side,other_side); CHKERR_PETSC(rval);
  rval = moab.add_entities(child_nodes_and_skin_edges,nodes); CHKERR_PETSC(rval);
  rval = moab.add_entities(child_nodes_and_skin_edges,skin_edges); CHKERR_PETSC(rval);
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"Nb. of side ents3d in set %u\n",side_ents3d.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. of other side ents3d in set %u\n",side_ents3d.size());
  }
  if(verb>3) {
    ierr = moab.write_file("side.vtk","VTK","",&children[0],1); CHKERRQ(ierr);
    ierr = moab.write_file("other_side.vtk","VTK","",&children[1],1); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_msId_3dENTS_split_sides(
  const EntityHandle meshset,const BitRefLevel &bit,
  const int msId,const Cubit_BC_bitset CubitBCType,const bool add_iterfece_entities,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  moabBaseMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    ierr = moabField_Core::get_msId_3dENTS_split_sides(
      meshset,bit,miit->meshset,add_iterfece_entities,recursive,verb); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_msId_3dENTS_split_sides(
  const EntityHandle meshset,const BitRefLevel &bit,
  const EntityHandle SideSet,const bool add_iterfece_entities,const bool recursive,
  int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  vector<EntityHandle> children;
  rval = moab.get_child_meshsets(SideSet,children);  CHKERR_PETSC(rval);
  if(children.size()!=2) SETERRQ(PETSC_COMM_SELF,1,"No children in set");
  vector<EntityHandle> children_nodes_and_skin_edges;
  rval = moab.get_child_meshsets(children[0],children_nodes_and_skin_edges);  CHKERR_PETSC(rval);
  if(children_nodes_and_skin_edges.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"No children in set");
  Range triangles;
  rval = moab.get_entities_by_type(SideSet,MBTRI,triangles,recursive);  CHKERR_PETSC(rval);
  Range side_ents3d;
  rval = moab.get_entities_by_type(children[0],MBTET,side_ents3d,false);  CHKERR_PETSC(rval);
  Range other_ents3d;
  rval = moab.get_entities_by_type(children[1],MBTET,other_ents3d,false);  CHKERR_PETSC(rval);
  Range nodes;
  rval = moab.get_entities_by_type(children_nodes_and_skin_edges[0],MBVERTEX,nodes,false);  CHKERR_PETSC(rval);
  /*Range edges;
  rval = moab.get_entities_by_type(children_nodes_and_skin_edges[0],MBEDGE,edges,false);  CHKERR_PETSC(rval);
  Range edges_nodes;
  rval = moab.get_connectivity(edges,edges_nodes,true); CHKERR_PETSC(rval);
  nodes = subtract(nodes,edges_nodes);*/
  if(verb>3) {
    PetscPrintf(PETSC_COMM_WORLD,"triangles %u\n",triangles.size());
    PetscPrintf(PETSC_COMM_WORLD,"side_ents3d %u\n",side_ents3d.size());
    PetscPrintf(PETSC_COMM_WORLD,"nodes %u\n",nodes.size());
  }
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ents_by_ent_type;
  ref_ents_by_ent_type &ref_ents_by_ent = refined_mofem_entities.get<MoABEnt_mi_tag>();
  map<EntityHandle,EntityHandle> map_nodes;
  Range::iterator nit = nodes.begin();
  double coord[3];
  for(;nit!=nodes.end();nit++) {
    rval = moab.get_coords(&*nit,1,coord); CHKERR_PETSC(rval);	
    EntityHandle new_node;
    rval = moab.create_vertex(coord,new_node); CHKERR(rval);
    map_nodes[*nit] = new_node;
    ref_ents_by_ent_type::iterator miit_ref_ent = ref_ents_by_ent.find(*nit);
    if(miit_ref_ent == ref_ents_by_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    rval = moab.tag_set_data(th_RefParentHandle,&new_node,1,&*nit); CHKERR_PETSC(rval);
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,new_node));
    refined_mofem_entities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(bit));
    refined_mofem_entities.modify(miit_ref_ent,RefMoFEMEntity_change_add_bit(bit));
  }
  EntityHandle meshset_for_bit_level;
  rval = moab.create_meshset(MESHSET_SET,meshset_for_bit_level); CHKERR_PETSC(rval);
  Range meshset_ents;
  rval = moab.get_entities_by_handle(meshset,meshset_ents,false); CHKERR_PETSC(rval);
  if(intersect(meshset_ents,side_ents3d).size() != side_ents3d.size()) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  rval = moab.add_entities(meshset_for_bit_level,subtract(meshset_ents,side_ents3d)); CHKERR_PETSC(rval);
  Range new_tets;
  Range::iterator tit = side_ents3d.begin();
  for(;tit!=side_ents3d.end();tit++) {
    ref_ents_by_ent_type::iterator miit_ref_ent = ref_ents_by_ent.find(*tit);
    if(miit_ref_ent==ref_ents_by_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"tet not in database");
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*tit,conn,num_nodes,true); CHKERR_PETSC(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    int ii = 0;
    for(; ii<num_nodes; ii++) {
      map<EntityHandle,EntityHandle>::iterator mit = map_nodes.find(conn[ii]);
      if(mit != map_nodes.end()) {
	new_conn[ii] = mit->second;
	if(verb>6) {
	  PetscPrintf(PETSC_COMM_WORLD,"nodes %u -> %d\n",conn[ii],new_conn[ii]);
	}
      } else {
	new_conn[ii] = conn[ii];
	nb_new_conn++;
      }
    }
    if(nb_new_conn==0) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
    //rval = moab.set_connectivity(*tit,new_conn,num_nodes); CHKERR_PETSC(rval);
    EntityHandle tet;
    rval = moab.create_element(MBTET,new_conn,4,tet); CHKERR_PETSC(rval);
    rval = moab.tag_set_data(th_RefParentHandle,&tet,1,&*tit); CHKERR_PETSC(rval);
    rval = moab.add_entities(meshset_for_bit_level,&tet,1); CHKERR_PETSC(rval);
    new_tets.insert(tet);
  }
  Range new_ents; 
  // create new entities by adjecies form new tets
  rval = moab.get_adjacencies(new_tets,1,true,new_ents,Interface::UNION); CHKERR_PETSC(rval);
  rval = moab.get_adjacencies(new_tets,2,true,new_ents,Interface::UNION); CHKERR_PETSC(rval);
  Range ents; 
  // edges and triangles
  rval = moab.get_adjacencies(triangles,1,false,ents,Interface::UNION); CHKERR_PETSC(rval);
  Range new_ents_in_database;
  ents.insert(triangles.begin(),triangles.end());
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    int ii = 0;
    for(;ii<num_nodes; ii++) {
      map<EntityHandle,EntityHandle>::iterator mit = map_nodes.find(conn[ii]);
      if(mit != map_nodes.end()) {
	new_conn[ii] = mit->second;
	nb_new_conn++;
	if(verb>6) {
	  PetscPrintf(PETSC_COMM_WORLD,"nodes %u -> %d\n",conn[ii],new_conn[ii]);
	}
      } else {
	new_conn[ii] = conn[ii];
      }
    }
    if(nb_new_conn==0) continue;
    ref_ents_by_ent_type::iterator miit_ref_ent = ref_ents_by_ent.find(*eit);
    if(miit_ref_ent == ref_ents_by_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
    Range new_ent;
    switch (moab.type_from_handle(*eit)) {
      case MBTRI: {
	  rval = moab.get_adjacencies(new_conn,3,2,true,new_ent); CHKERR_PETSC(rval);
	  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"new_ent %u\n",new_ent.size());
	  if(add_iterfece_entities) {
	    EntityHandle prism_conn[6] = { conn[0],conn[1],conn[2], new_conn[0],new_conn[1],new_conn[2] };
	    //cerr << conn[0] << " " << conn[1] << " " << conn[2] << " ::: " << new_conn[0] << " " << new_conn[1] << " " << new_conn[2] << endl;
	    EntityHandle prism = no_handle;
	    rval = moab.create_element(MBPRISM,prism_conn,6,prism); CHKERR_PETSC(rval);
	    ierr = add_prism_to_adjacencies_maps_for_prisms(prism,verb/*nb_new_conn < 3 ? 1 : 0*/); CHKERRQ(ierr);
	    rval = moab.add_entities(meshset_for_bit_level,&prism,1); CHKERR_PETSC(rval);
	  }
	} break;
      case MBEDGE: {
	  rval = moab.get_adjacencies(new_conn,2,1,true,new_ent); CHKERR_PETSC(rval);
	  if(new_ent.size()!=1) {
	    SETERRQ1(PETSC_COMM_SELF,1,"database insonistency, new_ent.size() = %u",new_ent.size());
	  }
	} break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
    }
    if(new_ent.size()!=1) {
      SETERRQ1(PETSC_COMM_SELF,1,"database insonistency, new_ent.size() = %u",new_ent.size());
    }
    if(new_ents.find(*new_ent.begin())==new_ents.end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
    rval = moab.tag_set_data(th_RefParentHandle,&*new_ent.begin(),1,&*eit); CHKERR_PETSC(rval);
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,new_ent[0]));
    refined_mofem_entities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(bit));
    new_ents_in_database.insert(new_ent.begin(),new_ent.end());
  }
  //all other entities, some ents like triangles and faces on the side of tets
  Range side_adj_faces_and_edges;
  rval = moab.get_adjacencies(side_ents3d,1,true,side_adj_faces_and_edges,Interface::UNION); CHKERR_PETSC(rval);
  rval = moab.get_adjacencies(side_ents3d,2,true,side_adj_faces_and_edges,Interface::UNION); CHKERR_PETSC(rval);
  side_adj_faces_and_edges = subtract(side_adj_faces_and_edges,new_ents_in_database);
  eit = side_adj_faces_and_edges.begin();
  for(;eit!=side_adj_faces_and_edges.end();eit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    int ii = 0;
    for(;ii<num_nodes; ii++) {
      map<EntityHandle,EntityHandle>::iterator mit = map_nodes.find(conn[ii]);
      if(mit != map_nodes.end()) {
	new_conn[ii] = mit->second;
	nb_new_conn++;
	if(verb>6) {
	  PetscPrintf(PETSC_COMM_WORLD,"nodes %u -> %d\n",conn[ii],new_conn[ii]);
	}
      } else {
	new_conn[ii] = conn[ii];
      }
    }
    if(nb_new_conn==0) continue;
    ref_ents_by_ent_type::iterator miit_ref_ent = ref_ents_by_ent.find(*eit);
    if(miit_ref_ent == ref_ents_by_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
    Range new_ent;
    switch (moab.type_from_handle(*eit)) {
      case MBTRI: {
	  rval = moab.get_adjacencies(new_conn,3,2,true,new_ent); CHKERR_PETSC(rval);
	}
	break;
      case MBEDGE: {
	  rval = moab.get_adjacencies(new_conn,2,1,true,new_ent); CHKERR_PETSC(rval);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"database insonistency");   
    }
    if(new_ent.size()!=1) {
      SETERRQ1(PETSC_COMM_SELF,1,"database insonistency, new_ent.size() = %u",new_ent.size());
    }
    if(new_ents.find(*new_ent.begin())==side_adj_faces_and_edges.end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
    rval = moab.tag_set_data(th_RefParentHandle,&*new_ent.begin(),1,&*eit); CHKERR_PETSC(rval);
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,new_ent[0]));
    refined_mofem_entities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(bit));
    if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"new_ent %u\n",new_ent.size());
    new_ents_in_database.insert(new_ent.begin(),new_ent.end());
  }
  //finalise by adding new tets and prism ti bitlelvel
  ierr = seed_ref_level_3D(meshset_for_bit_level,bit); CHKERRQ(ierr);
  rval = moab.delete_entities(&meshset_for_bit_level,1); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::add_prism_to_adjacencies_maps_for_prisms(const EntityHandle prism,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  vector<EntityHandle> Ents(8,no_handle);
  try {
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_mofem_entities.insert(RefMoFEMEntity(moab,prism));
    pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
    if(p_ent.second) {
      p_MoFEMFiniteElement = refined_mofem_elements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
      int num_nodes;
      const EntityHandle* conn;
      rval = moab.get_connectivity(prism,conn,num_nodes,true); CHKERR_THROW(rval);
      Range face_side3,face_side4;
      rval = moab.get_adjacencies(conn,3,2,false,face_side3); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(&conn[3],3,2,false,face_side4); CHKERR_PETSC(rval);
      if(face_side3.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prism don't have side face 3");
      if(face_side4.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prims don't have side face 4");
      p_MoFEMFiniteElement.first->get_side_number_ptr(moab,*face_side3.begin());
      p_MoFEMFiniteElement.first->get_side_number_ptr(moab,*face_side4.begin());
      // set bit common for faces with side number 3 and 4
      RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refined_mofem_entities.get<MoABEnt_mi_tag>().find(*face_side3.begin());
      if(miit_ref_ent!=refined_mofem_entities.get<MoABEnt_mi_tag>().end()) {
	BitRefLevel bit = miit_ref_ent->get_BitRefLevel();
	if(face_side4.empty()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
	miit_ref_ent = refined_mofem_entities.get<MoABEnt_mi_tag>().find(*face_side4.begin());
	if(miit_ref_ent!=refined_mofem_entities.get<MoABEnt_mi_tag>().end()) {
	  bit &= miit_ref_ent->get_BitRefLevel();
	  refined_mofem_entities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	}
      }
    } 
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::loop_finite_elements(const string &problem_name,const string &fe_name,FEMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  //
  ierr = loop_finite_elements(problem_name,fe_name,method,pcomm->rank(),pcomm->rank(),verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::loop_finite_elements(
  const string &problem_name,const string &fe_name,FEMethod &method,int lower_rank,int upper_rank,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  // find p_miit
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(problem_name);
  if(p_miit == problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem not in databse %s",problem_name.c_str());
  // finite element
  typedef NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type FEs_by_composite;
  FEs_by_composite &numered_finite_elements = 
    (const_cast<NumeredMoFEMFiniteElement_multiIndex&>(p_miit->numered_finite_elements)).get<Composite_mi_tag>();
  FEs_by_composite::iterator miit = numered_finite_elements.lower_bound(boost::make_tuple(fe_name,lower_rank));
  FEs_by_composite::iterator hi_miit = numered_finite_elements.upper_bound(boost::make_tuple(fe_name,upper_rank));
  method.fe_name = fe_name;
  ierr = method.set_problem(&*p_miit); CHKERRQ(ierr);
  ierr = method.set_moabfields(&moabfields); CHKERRQ(ierr);
  ierr = method.set_ents_multiIndex(&ents_moabfield); CHKERRQ(ierr);
  ierr = method.set_dofs_multiIndex(&dofs_moabfield); CHKERRQ(ierr);
  ierr = method.set_fes_multiIndex(&finite_elements); CHKERRQ(ierr);
  ierr = method.set_fes_data_multiIndex(&finite_elements_moabents); CHKERRQ(ierr);
  ierr = method.set_adjacencies(&adjacencies); CHKERRQ(ierr);
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);
  for(;miit!=hi_miit;miit++) {
    ierr = method.set_fe(&*miit); CHKERRQ(ierr);
    ierr = method.set_data_multIndex( const_cast<FEDofMoFEMEntity_multiIndex*>(&(miit->fe_ptr->data_dofs)) ); CHKERRQ(ierr);
    ierr = method.set_row_multIndex( const_cast<FENumeredDofMoFEMEntity_multiIndex*>(&(miit->rows_dofs)) ); CHKERRQ(ierr);
    ierr = method.set_col_multIndex( const_cast<FENumeredDofMoFEMEntity_multiIndex*>(&(miit->cols_dofs)) ); CHKERRQ(ierr);
    try {
      PetscLogEventBegin(USER_EVENT_operator,0,0,0,0);
      ierr = method(); CHKERRQ(ierr);
      PetscLogEventEnd(USER_EVENT_operator,0,0,0,0);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  }
  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type dofs_by_name;
  // find p_miit
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(problem_name);
  if(p_miit == problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem not in databse %s",problem_name.c_str());
  ierr = method.set_moabfields(&moabfields); CHKERRQ(ierr);
  ierr = method.set_ents_multiIndex(&ents_moabfield); CHKERRQ(ierr);
  ierr = method.set_dofs_multiIndex(&dofs_moabfield); CHKERRQ(ierr);
  ierr = method.set_fes_multiIndex(&finite_elements); CHKERRQ(ierr);
  ierr = method.set_fes_data_multiIndex(&finite_elements_moabents); CHKERRQ(ierr);
  ierr = method.set_adjacencies(&adjacencies); CHKERRQ(ierr);
  dofs_by_name *dofs;
  switch (rc) {
    case Row:
      dofs = const_cast<dofs_by_name*>(&p_miit->numered_dofs_rows.get<FieldName_mi_tag>());
      break;
    case Col:
      dofs = const_cast<dofs_by_name*>(&p_miit->numered_dofs_cols.get<FieldName_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  dofs_by_name::iterator miit = dofs->lower_bound(field_name);
  dofs_by_name::iterator hi_miit = dofs->upper_bound(field_name);
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(;miit!=hi_miit;miit++) {
    ierr = method.set_dof(&*miit); CHKERRQ(ierr);
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_problems_database(const string &problem_name,const MoFEMProblem **problem_ptr) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type problems_by_name;
  problems_by_name &problems_set = problems.get<MoFEMProblem_mi_tag>();
  problems_by_name::iterator p_miit = problems_set.find(problem_name);
  if(p_miit == problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found, (top tip: check spelling)",problem_name.c_str());
  *problem_ptr = &*p_miit;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_dofs_moabfield(const DofMoFEMEntity_multiIndex **dofs_moabfield_ptr) {
  PetscFunctionBegin;
  *dofs_moabfield_ptr = &dofs_moabfield;
  PetscFunctionReturn(0);
}
PetscErrorCode moabField_Core::get_finite_elements(const MoFEMFiniteElement_multiIndex **finite_elements_ptr) {
  PetscFunctionBegin;
  *finite_elements_ptr = &finite_elements;
  PetscFunctionReturn(0);
}


}
