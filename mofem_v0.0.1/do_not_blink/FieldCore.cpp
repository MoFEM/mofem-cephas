/** \file FieldCore.cpp
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

#include<FieldCore.hpp>
#include<FEM.h>

namespace MoFEM {

const static int debug = 1;

FieldCore::FieldCore(Interface& _moab,int _verbose): 
  moab(_moab),verbose(_verbose) {
  const EntityHandle root_meshset = moab.get_root_set();
  // Version
  stringstream strs_version;
  strs_version << "MoFEM_version_" << MoFEM_VERSION_MAJOR << "." << MoFEM_VERSION_MINOR;
  Tag th_version;
  string version = strs_version.str();
  rval = moab.tag_get_handle("_MoFEM_VERSION",version.size()*sizeof(char),MB_TYPE_OPAQUE,
    th_version,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,NULL); CHKERR_THROW(rval);
  const char *ptr_version = version.c_str();
  rval = moab.tag_set_data(th_version,&root_meshset,1,ptr_version); CHKERR_THROW(rval);
  //tags saved in vtk-files
  const int def_part = -1;
  rval = moab.tag_get_handle("PARTITION",1,MB_TYPE_INTEGER,th_Part,MB_TAG_CREAT|MB_TAG_SPARSE,&def_part); CHKERR_THROW(rval);
  //Tags Ref
  rval = moab.tag_get_handle("_RefParentHandle",1,MB_TYPE_HANDLE,th_RefParentHandle,MB_TAG_CREAT|MB_TAG_SPARSE,&root_meshset); CHKERR_THROW(rval); CHKERR_THROW(rval);
  const int def_type[] = {0,0};
  rval = moab.tag_get_handle("_RefType",2,MB_TYPE_INTEGER,th_RefType,MB_TAG_CREAT|MB_TAG_SPARSE,def_type); CHKERR_THROW(rval);
  BitRefLevel def_bit_level = 0;
  rval = moab.tag_get_handle("_RefBitLevel",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_RefBitLevel,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); CHKERR_THROW(rval);
  BitRefEdges def_bit_egde = 0;
  rval = moab.tag_get_handle("_RefBitEdge",sizeof(BitRefEdges),MB_TYPE_OPAQUE,
    th_RefBitEdge,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_egde); CHKERR_THROW(rval);
  //Tags Field
  const unsigned long int def_id = 0;
  rval = moab.tag_get_handle("_FieldId",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FieldId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  FieldSpace def_space = LastSpace;
  rval = moab.tag_get_handle("_FieldSpace",sizeof(FieldSpace),MB_TYPE_OPAQUE,
    th_FieldSpace,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_space); CHKERR_THROW(rval);
  const int def_val_len = 0;
  rval = moab.tag_get_handle("_FieldName",def_val_len,MB_TYPE_OPAQUE,
    th_FieldName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FieldName_DataNamePrefix",def_val_len,MB_TYPE_OPAQUE,
    th_FieldName_DataNamePrefix,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  //Tags FE
  rval = moab.tag_get_handle("_FEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_FEId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FEName",def_val_len,MB_TYPE_OPAQUE,
    th_FEName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FEIdCol",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FEIdRow",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FEIdData",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  //Tags Problem
  rval = moab.tag_get_handle("_ProblemId",sizeof(BitProblemId),MB_TYPE_OPAQUE,
    th_ProblemId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemFEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_ProblemFEId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemName",def_val_len,MB_TYPE_OPAQUE,
    th_ProblemName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  DofIdx def_nbdofs = 0;
  rval = moab.tag_get_handle("_ProblemNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemLocalNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemGhostNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemLocalNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemGhostNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  //Global Variables
  //Fields
  int def_shift = 1;
  rval = moab.tag_get_handle("_FieldShift",1,MB_TYPE_INTEGER,th_FieldShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR_THROW(rval);
  const void* tag_data[1];
  rval = moab.tag_get_by_ptr(th_FieldShift,&root_meshset,1,tag_data); CHKERR_THROW(rval);
  f_shift = (int*)tag_data[0];
  //FE
  rval = moab.tag_get_handle("_FEShift",1,MB_TYPE_INTEGER,th_FEShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_FEShift,&root_meshset,1,tag_data); CHKERR_THROW(rval);
  MoFEMFiniteElement_shift = (int*)tag_data[0];
  //Problem
  rval = moab.tag_get_handle("_ProblemShift",1,MB_TYPE_INTEGER,th_ProblemShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemShift,&root_meshset,1,tag_data); CHKERR_THROW(rval);
  p_shift = (int*)tag_data[0];
  //SaftyNets
  int def_bool = 0;
  rval = moab.tag_get_handle("_MoFEMBuild",1,MB_TYPE_INTEGER,th_MoFEMBuild,MB_TAG_CREAT|MB_TAG_MESH,&def_bool); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  rval = moab.tag_get_by_ptr(th_MoFEMBuild,&root_meshset,1,(const void **)&build_MoFEM); CHKERR_THROW(rval);
  //Meshsets
  int default_val = -1;
  rval = moab.tag_get_handle(DIRICHLET_SET_TAG_NAME,1, MB_TYPE_INTEGER, 
    nsTag, MB_TAG_SPARSE|MB_TAG_CREAT, &default_val); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(NEUMANN_SET_TAG_NAME,1, MB_TYPE_INTEGER, 
    ssTag, MB_TAG_SPARSE|MB_TAG_CREAT, &default_val); CHKERR_THROW(rval);
  const int def_bc_data_len = 0;
  std::string tag_name = std::string(DIRICHLET_SET_TAG_NAME)+"__BC_DATA";
  rval = moab.tag_get_handle(tag_name.c_str(),def_bc_data_len,MB_TYPE_OPAQUE,
    nsTag_data,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  tag_name = std::string(NEUMANN_SET_TAG_NAME)+"__BC_DATA";
  rval = moab.tag_get_handle(tag_name.c_str(),def_bc_data_len,MB_TYPE_OPAQUE,
    ssTag_data,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER,
    bhTag,MB_TAG_SPARSE|MB_TAG_CREAT,&default_val); CHKERR_THROW(rval);
  std::vector<unsigned int> def_uint_zero(12,0);
  rval= moab.tag_get_handle("BLOCK_HEADER",12*sizeof(unsigned int),MB_TYPE_INTEGER,
    bhTag_header,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_uint_zero[0]); CHKERR_THROW(rval); 
  Tag block_attribs;
  int def_Block_Attributes_length = 0;
  rval = moab.tag_get_handle("Block_Attributes",def_Block_Attributes_length,MB_TYPE_DOUBLE,
    block_attribs,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval); 
  Tag entity_name_tag;
  rval = moab.tag_get_handle(
    NAME_TAG_NAME,NAME_TAG_SIZE,MB_TYPE_OPAQUE,entity_name_tag,MB_TAG_SPARSE|MB_TAG_CREAT); CHKERR_THROW(rval);
  // For VTK files
  int def_elem_type = MBMAXTYPE;
  rval = moab.tag_get_handle("ElemType",1,MB_TYPE_INTEGER,th_ElemType,MB_TAG_CREAT|MB_TAG_SPARSE,&def_elem_type); CHKERR_THROW(rval); 
  //
  initialiseDatabseInformationFromMesh(verbose); 
  //
  ShapeDiffMBTET(diffN_TET); 
  // Petsc Logs
  PetscLogEventRegister("FE_preProcess",0,&USER_EVENT_preProcess);
  PetscLogEventRegister("FE_operator",0,&USER_EVENT_operator);
  PetscLogEventRegister("FE_postProcess",0,&USER_EVENT_postProcess);
}
FieldCore::~FieldCore() {
}
Interface& FieldCore::get_moab() {
  return moab;
}
BitFieldId FieldCore::get_BitFieldId(const string& name) const {
  typedef MoFEMField_multiIndex::index<FieldName_mi_tag>::type field_set_by_name;
  const field_set_by_name &set = moabFields.get<FieldName_mi_tag>();
  field_set_by_name::iterator miit = set.find(name);
  if(miit==set.end()) THROW_AT_LINE("field < "+name+" > not in databse (top tip: check spelling)");
  return miit->get_id();
}
string FieldCore::get_BitFieldId_name(const BitFieldId id) const {
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set.find(id);
  return miit->get_name();
}
EntityHandle FieldCore::get_field_meshset(const BitFieldId id) const {
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set.find(id);
  if(miit==set.end()) THROW_AT_LINE("field not in databse (top tip: check spelling)");
  return miit->meshset;
}
EntityHandle FieldCore::get_field_meshset(const string& name) const {
  return get_field_meshset(get_BitFieldId(name));
}
BitFieldId FieldCore::get_field_shift() {
  assert((unsigned int)*f_shift<BitFieldId().set().to_ulong());
  return (BitFieldId)(1<<(((*f_shift)++)-1)); 
}
BitFEId FieldCore::get_BitFEId() {
  assert((unsigned int)*MoFEMFiniteElement_shift<BitFEId().set().to_ulong());
  return BitFEId(1<<(((*MoFEMFiniteElement_shift)++)-1)); 
}
BitProblemId FieldCore::get_problem_shift() {
  assert((unsigned int)*p_shift<BitProblemId().set().to_ulong());
  return BitProblemId(1<<(((*p_shift)++)-1)); 
}
PetscErrorCode FieldCore::clear_map() {
  PetscFunctionBegin;
  cubit_meshsets.clear();
  refinedMoFemEntities.clear();
  refinedMoFemElements.clear();
  moabFields.clear();
  entsMoabField.clear();
  dofsMoabField.clear();
  finiteElements.clear();
  finiteElementsMoFEMEnts.clear();
  entFEAdjacencies.clear();
  moFEMProblems.clear();
  PetscFunctionReturn(0);
} 
PetscErrorCode FieldCore::add_field(const string& name,const BitFieldId id,const FieldSpace space,const ApproximationRank rank,enum MoFEMTypes bh,int verb) {
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
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,meshset));
    bool success = refinedMoFemEntities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(BitRefLevel().set()));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  }
  //name
  void const* tag_data[] = { name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = name.size();
  rval = moab.tag_set_by_ptr(th_FieldName,&meshset,1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  //name data prefix
  string name_data_prefix("_App_Data");
  void const* tag_prefix_data[] = { name_data_prefix.c_str() };
  int tag_prefix_sizes[1]; tag_prefix_sizes[0] = name_data_prefix.size();
  rval = moab.tag_set_by_ptr(th_FieldName_DataNamePrefix,&meshset,1,tag_prefix_data,tag_prefix_sizes); CHKERR_PETSC(rval);
  Tag th_AppOrder,th_FieldData,th_Rank,th_AppDofOrder,th_DofRank;
  //data
  string Tag_data_name = name_data_prefix+name;
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
    p = moabFields.insert(MoFEMField(moab,meshset));  
    if(bh == MF_EXCL) {
      if(!p.second) SETERRQ1(PETSC_COMM_SELF,1,
	"field not inesrted %s (top tip, it could be already there)",
	MoFEMField(moab,meshset).get_name().c_str());
    }
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
PetscErrorCode FieldCore::add_field(const string& name,const FieldSpace space,const ApproximationRank rank,enum MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  *build_MoFEM = 0;
  if(verb==-1) verb = verbose;
  BitFieldId id = get_field_shift();
  ierr = add_field(name,id,space,rank,bh,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::initialiseDatabseInformationFromMesh(int verb) {
  PetscFunctionBegin;
  //ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(verb==-1) verb = verbose;
  ierr = clear_map(); CHKERRQ(ierr);
  Range meshsets;
  rval = moab.get_entities_by_type(0,MBENTITYSET,meshsets,false);  CHKERR_PETSC(rval);
  //loop all meshsets in moab database
  Range::iterator mit = meshsets.begin();
  for(;mit!=meshsets.end();mit++) {
    try {
      //check if meshset is cubit meshset
      CubitMeshSets base_meshset(moab,*mit);
      if((base_meshset.CubitBCType&Cubit_BC_bitset(NodeSet|SideSet|BlockSet)).any()) {
	pair<moabCubitMeshSet_multiIndex::iterator,bool> p = cubit_meshsets.insert(base_meshset);
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
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,1,msg);
    }
    BitFieldId field_id;
    //get bit id form field tag
    rval = moab.tag_get_data(th_FieldId,&*mit,1,&field_id); CHKERR_PETSC(rval);
    //check if meshset if field meshset
    if(field_id!=0) {
      pair<MoFEMField_multiIndex::iterator,bool> p;
      try {
	p = moabFields.insert(MoFEMField(moab,*mit));
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
	pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*mit));
	NOT_USED(p_ref_ent);
      } else {
	Range ents;
	rval = moab.get_entities_by_handle(*mit,ents,false); CHKERR_PETSC(rval);
	if(verb > 1) {
	  ostringstream ss;
	  ss << "read field ents " << ents.size() << endl;;
	  PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	}
	Range::iterator eit = ents.begin();
	for(;eit!=ents.end();eit++) {
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*eit));
	  try {
	    MoFEMEntity moabent(moab,&*p.first,&*p_ref_ent.first);
	    if(moabent.get_order_nb_dofs(moabent.get_max_order())==0) continue; 
	    pair<MoFEMEntity_multiIndex::iterator,bool> p_ent = entsMoabField.insert(moabent);
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
    //get bit id from fe tag
    rval = moab.tag_get_data(th_FEId,&*mit,1,&fe_id); CHKERR_PETSC(rval);
    //check if meshset is finite element meshset
    if(fe_id!=0) {
      pair<MoFEMFiniteElement_multiIndex::iterator,bool> p = finiteElements.insert(MoFEMFiniteElement(moab,*mit));
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
	pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*eit));
	pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
	switch (moab.type_from_handle(*eit)) {
	  case MBVERTEX:
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_VERTEX(moab,&*p_ref_ent.first)));
	  case MBEDGE:
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_EDGE(moab,&*p_ref_ent.first)));
	  case MBTRI:
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TRI(moab,&*p_ref_ent.first)));
	  case MBTET:
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_ref_ent.first)));
	    break;
	  case MBPRISM:
  	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ref_ent.first)));
	    break;
	  case MBENTITYSET:
  	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_MESHSET(moab,&*p_ref_ent.first)));
	    break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"Only finite elements of type MBTET, MBPRISM and MBENTITYSET are implemented");
	}
      }
    }
    BitProblemId problem_id;
    //get bit id form problem tag
    rval = moab.tag_get_data(th_ProblemId,&*mit,1,&problem_id); CHKERR_PETSC(rval);
    //check if meshset if problem meshset
    if(problem_id!=0) {
      pair<MoFEMProblem_multiIndex::iterator,bool> p = moFEMProblems.insert(MoFEMProblem(moab,*mit));
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
    ierr = add_prism_to_basicEntAdjacencies(*pit); CHKERRQ(ierr);
  }
  if(verb > 2) {
    list_field();
    list_finite_elements();
    list_problem();
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_field_by_EDGEs(const EntityHandle meshset,const BitFieldId id,int verb) {
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
  Range nodes,edges;
  rval = moab.get_entities_by_type(meshset,MBEDGE,edges,true); CHKERR_PETSC(rval);
  switch (space) {
    case L2:
      rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
      if(verb>1) {
	ostringstream ss;
	ss << "add entities to field " << get_BitFieldId_name(id);
	ss << " nb. add edges " << edges.size();
	ss << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      break;
    case H1:
      rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
      rval = moab.get_connectivity(edges,nodes,true); CHKERR_PETSC(rval);
      rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
      if(verb>1) {
	ostringstream ss;
	ss << "add entities to field " << get_BitFieldId_name(id);
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
PetscErrorCode FieldCore::add_ents_to_field_by_EDGEs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_EDGEs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_field_by_TRIs(const EntityHandle meshset,const BitFieldId id,int verb) {
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
  rval = moab.get_entities_by_type(meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
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
PetscErrorCode FieldCore::add_ents_to_field_by_TRIs(const EntityHandle meshset,const string& name,int verb) {
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
PetscErrorCode FieldCore::add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const BitFieldId id,int verb) {
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
  Range nodes;
  rval = moab.get_entities_by_type(meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
  switch (space) {
    case H1:
      rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
      if(verb>1) {
	ostringstream ss;
	ss << "add entities to field " << get_BitFieldId_name(id);
	ss << " nb. add nodes " << nodes.size();
	ss << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"add_ents_to_field_by_TRIs this field not work for TRIs");
  }
  ierr = seed_ref_level_3D(idm,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_VERTICEs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_field_by_TETs(const EntityHandle meshset,const BitFieldId id,int verb) {
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
  rval = moab.get_entities_by_type(meshset,MBTET,tets,true); CHKERR_PETSC(rval);
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
PetscErrorCode FieldCore::add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name,int verb) {
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
PetscErrorCode FieldCore::set_field_order(const EntityHandle meshset,const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents); CHKERR_PETSC(rval);
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents for order change %d\n",ents.size());
  }
  //check field & meshset
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
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
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents for order change in the field %d\n",ents.size());
  }
  vector<const void*> tag_data_order(ents.size());
  rval = moab.tag_get_by_ptr(miit->th_AppOrder,ents,&tag_data_order[0]); CHKERR_PETSC(rval);
  //ent view by field id (in set all MoabEnts has the same FieldId)
  typedef MoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type ent_set_by_id;
  ent_set_by_id& set = entsMoabField.get<BitFieldId_mi_tag>();
  ent_set_by_id::iterator miit2 = set.lower_bound(id);
  MoFEMEntity_multiIndex_ent_view ents_id_view;
  if(miit2 != set.end()) {
    ent_set_by_id::iterator hi_miit2 = set.upper_bound(id);
    for(;miit2!=hi_miit2;miit2++) {
      ents_id_view.insert(&*miit2);
    }
  }
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents in the multi index field %d\n",ents_id_view.size());
  }
  //loop over ents
  int nb_ents_set_order_up = 0;
  int nb_ents_set_order_down = 0;
  int nb_ents_set_order_new = 0;
  Range::iterator eit = ents.begin();
  for(unsigned int ee = 0;ee<ents.size();ee++,eit++) {
    MoFEMEntity_multiIndex_ent_view::iterator miit3 = ents_id_view.find(*eit);
    if(miit3!=ents_id_view.end()) {
      const ApproximationOrder old_ApproximationOrder = (*miit3)->get_max_order();
      if(old_ApproximationOrder==order) continue;
      MoFEMEntity_multiIndex::iterator miit4 = entsMoabField.get<Unique_mi_tag>().find((*miit3)->get_unique_id());
      assert(miit4!=entsMoabField.end());
      if(miit4->get_max_order()<order) nb_ents_set_order_up++;
      if(miit4->get_max_order()>order) nb_ents_set_order_down++;
      typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent>::type dof_set_type;
      dof_set_type& set_set = dofsMoabField.get<Composite_Name_And_Ent>();
      dof_set_type::iterator miit5 = set_set.lower_bound(boost::make_tuple(miit4->get_name_ref(),miit4->get_ent()));
      dof_set_type::iterator hi_miit6 = set_set.upper_bound(boost::make_tuple(miit4->get_name_ref(),miit4->get_ent()));
      for(;miit5!=hi_miit6;miit5++) {
	if(miit5->get_dof_order()<=order) continue;
	bool success = dofsMoabField.modify(dofsMoabField.project<0>(miit5),DofMoFEMEntity_active_change(false));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }
      bool success = entsMoabField.modify(entsMoabField.project<0>(miit4),MoFEMEntity_change_order(moab,order));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    } else {
      *(ApproximationOrder*)tag_data_order[ee] = order;
      RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*eit);
      if(miit_ref_ent==refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
      try { 
	MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
	//if(moabent.get_order_nb_dofs(moabent.get_max_order())==0) continue; 
	pair<MoFEMEntity_multiIndex::iterator,bool> e_miit = entsMoabField.insert(moabent);
	bool success = entsMoabField.modify(e_miit.first,MoFEMEntity_change_order(moab,order));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
	nb_ents_set_order_new++;
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
    }
  }
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents for which order was incresed %d (order %d)\n",nb_ents_set_order_up,order);
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents for which order was reduced %d (order %d)\n",nb_ents_set_order_down,order);
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents for which order set %d (order %d)\n",nb_ents_set_order_new,order);
  }
  if(verb>4) {
    list_ent_by_id(id);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::set_field_order(const EntityHandle meshset,const EntityType type,const string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try{
    ierr = set_field_order(meshset,type,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::dofs_NoField(const BitFieldId id,int &dof_counter) {
  PetscFunctionBegin;
  //field it
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  //find fiels
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) SETERRQ(PETSC_COMM_SELF,1,"field no found");
  //serch if field meshset is in database
  RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(miit->meshset);
  if(miit_ref_ent==refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
  pair<MoFEMEntity_multiIndex::iterator,bool> e_miit;
  try {
    //create database entity
    MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
    e_miit = entsMoabField.insert(moabent);
    //this is nor real field in space (set order to zero)
    bool success = entsMoabField.modify(e_miit.first,MoFEMEntity_change_order(moab,0));
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
    d_miit.first = dofsMoabField.project<0>(
      dofsMoabField.get<Unique_mi_tag>().find(DofMoFEMEntity::get_unique_id_calculate(rank,&*(e_miit.first)))
    );
    //if dof is not in databse
    if(d_miit.first==dofsMoabField.end()) {
      //insert dof
      d_miit = dofsMoabField.insert(DofMoFEMEntity(&*(e_miit.first),0,rank,rank));
      if(d_miit.second) dof_counter++;
      bool success = dofsMoabField.modify(d_miit.first,DofMoFEMEntity_active_change(true));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    //check consistency
    assert(d_miit.first->get_ent()==e_miit.first->get_MoFEMField_ptr()->meshset);
    assert(d_miit.first->get_ent_type()==e_miit.first->get_ent_type());
    assert(d_miit.first->get_id()==e_miit.first->get_id());
  }
  if(verbose>2) {
    typedef DofMoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type dof_set_by_id;
    dof_set_by_id &set = dofsMoabField.get<BitFieldId_mi_tag>();
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
PetscErrorCode FieldCore::dofs_L2H1HcurlHdiv(const BitFieldId id,int &dof_counter,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ents_by_ents;
  //find field
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) SETERRQ(PETSC_COMM_SELF,1,"field no found");
  //ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(miit->meshset,ents_of_id_meshset,false); CHKERR_PETSC(rval);
  //create dofsMoabField
  Range::iterator eit = ents_of_id_meshset.begin();
  for(;eit!=ents_of_id_meshset.end();eit++) {
    // check if ent is in ref meshset
    RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*eit);
    if(miit_ref_ent==refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
    //pair<MoFEMEntity_multiIndex::iterator,bool> e_miit;
    MoFEMEntity_multiIndex::iterator e_miit;
    try {
      e_miit = entsMoabField.find(MoFEMEntity(moab,&*miit,&*miit_ref_ent).get_unique_id());
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,1,msg);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss <<  ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    // create mofem entity linked to ref ent
    try {
      e_miit = entsMoabField.find(MoFEMEntity(moab,&*miit,&*miit_ref_ent).get_unique_id());
    } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    if(e_miit == entsMoabField.end()) {
      ApproximationOrder order = -1;
      rval = moab.tag_set_data(miit->th_AppOrder,&*eit,1,&order); CHKERR_PETSC(rval);
      pair<MoFEMEntity_multiIndex::iterator,bool> p_e_miit;
      try {
	MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
	p_e_miit = entsMoabField.insert(moabent);
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      if(!p_e_miit.second) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      bool success = entsMoabField.modify(p_e_miit.first,MoFEMEntity_change_order(moab,-1));
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
	    d_miit = dofsMoabField.insert(mdof);
	    if(d_miit.second) {
	      dof_counter++;
	      bool success = dofsMoabField.modify(d_miit.first,DofMoFEMEntity_active_change(true));
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
PetscErrorCode FieldCore::build_fields(int verb) {
  //PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
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
  PetscPrintf(PETSC_COMM_WORLD,"Nb. dofs %u\n",dofsMoabField.size());
  *build_MoFEM = 1<<0;
  //PetscFunctionReturn(0);
  return 0;
}
PetscErrorCode FieldCore::list_dof_by_id(const BitFieldId id) const {
  PetscFunctionBegin;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank %d Field dofs %s\n",pcomm->rank(),get_BitFieldId_name(id).c_str()); //!!! Syn ...
  typedef DofMoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type dof_set_by_id;
  const dof_set_by_id &set_id = dofsMoabField.get<BitFieldId_mi_tag>();
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
PetscErrorCode FieldCore::list_ent_by_id(const BitFieldId id) const {
  PetscFunctionBegin;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank %d Field ents %s\n",pcomm->rank(),get_BitFieldId_name(id).c_str()); //!!!! Syn ...
  typedef MoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type ent_set_by_id;
  const ent_set_by_id &set = entsMoabField.get<BitFieldId_mi_tag>();
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
PetscErrorCode FieldCore::list_field() const {
  PetscFunctionBegin;
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_finite_element(const string &MoFEMFiniteElement_name,enum MoFEMTypes bh) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<MoFEMFiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(bh == MF_EXCL) {
    if(it_MoFEMFiniteElement!=MoFEMFiniteElement_name_set.end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"this < %s > is there",MoFEMFiniteElement_name.c_str());
    }
  } else {
    if(it_MoFEMFiniteElement!=MoFEMFiniteElement_name_set.end()) PetscFunctionReturn(0);
  }
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
  //id
  BitFEId id = get_BitFEId();
  rval = moab.tag_set_data(th_FEId,&meshset,1,&id); CHKERR_PETSC(rval);
  //id name
  void const* tag_data[] = { MoFEMFiniteElement_name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = MoFEMFiniteElement_name.size();
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
  pair<MoFEMFiniteElement_multiIndex::iterator,bool> p = finiteElements.insert(MoFEMFiniteElement(moab,meshset));
  if(!p.second) SETERRQ(PETSC_COMM_SELF,1,"MoFEMFiniteElement not inserted");
  if(verbose>0) {
    ostringstream ss;
    ss << "add finite element: " << MoFEMFiniteElement_name << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    //list_finiteElements();
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::modify_finite_element_add_field_data(const string &MoFEMFiniteElement_name,const string &name_data) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<MoFEMFiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,1,"this MoFEMFiniteElement is there");
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,EntMoFEMFiniteElement_change_bit_add(get_BitFieldId(name_data)));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::modify_finite_element_add_field_row(const string &MoFEMFiniteElement_name,const string &name_row) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<MoFEMFiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"this < %s > is nor there",MoFEMFiniteElement_name.c_str());
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_row_change_bit_add(get_BitFieldId(name_row)));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::modify_finite_element_add_field_col(const string &MoFEMFiniteElement_name,const string &name_col) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<MoFEMFiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,1,"this MoFEMFiniteElement is there");
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_col_change_bit_add(get_BitFieldId(name_col)));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::modify_finite_element_off_field_data(const string &MoFEMFiniteElement_name,const string &name_data) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<MoFEMFiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,1,"this MoFEMFiniteElement is there");
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,EntMoFEMFiniteElement_change_bit_off(get_BitFieldId(name_data)));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::modify_finite_element_off_field_row(const string &MoFEMFiniteElement_name,const string &name_row) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<MoFEMFiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"this < %s > is nor there",MoFEMFiniteElement_name.c_str());
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_row_change_bit_off(get_BitFieldId(name_row)));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::modify_finite_element_off_field_col(const string &MoFEMFiniteElement_name,const string &name_col) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<MoFEMFiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,1,"this MoFEMFiniteElement is there");
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_col_change_bit_off(get_BitFieldId(name_col)));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
BitFEId FieldCore::get_BitFEId(const string& name) const {
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finiteElements_by_name;
  const finiteElements_by_name& set = finiteElements.get<MoFEMFiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator miit = set.find(name);
  if(miit==set.end()) THROW_AT_LINE(("finite element < "+name+" > not found (top tip: check spelling)").c_str());
  return miit->get_id();
}
string FieldCore::get_BitFEId_name(const BitFEId id) const {
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
  const finiteElements_by_id& set = finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = set.find(id);
  assert(miit!=set.end());
  return miit->get_name();
}
EntityHandle FieldCore::get_meshset_by_BitFEId(const BitFEId id) const {
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
  const finiteElements_by_id& set = finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = set.find(id);
  if(miit==set.end()) THROW_AT_LINE("finite element not found");
  return miit->meshset;
}
EntityHandle FieldCore::get_meshset_by_BitFEId(const string& name) const {	
  return get_meshset_by_BitFEId(get_BitFEId(name));
}
PetscErrorCode FieldCore::list_finite_elements() const {
  PetscFunctionBegin;
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
  const finiteElements_by_id &BitFEId_set = finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = BitFEId_set.begin();
  for(;miit!=BitFEId_set.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_problem(const BitProblemId id,const string& name) {
  PetscFunctionBegin;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
  rval = moab.tag_set_data(th_ProblemId,&meshset,1,&id); CHKERR_PETSC(rval);
  void const* tag_data[] = { name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = name.size();
  rval = moab.tag_set_by_ptr(th_ProblemName,&meshset,1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  //create entry
  pair<MoFEMProblem_multiIndex::iterator,bool> p = moFEMProblems.insert(MoFEMProblem(moab,meshset));
  NOT_USED(p);
  assert(p.second);
  if(verbose>0) {
    ostringstream ss;
    ss << "add problem: " << name << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_problem(const string& name) {
  PetscFunctionBegin;
  BitProblemId id = get_problem_shift();
  ierr = add_problem(id,name); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
BitProblemId FieldCore::get_BitProblemId(const string& name) const {
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  const moFEMProblems_by_name& set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator miit = set.find(name);
  return miit->get_id();
}
PetscErrorCode FieldCore::list_problem() const {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<BitProblemId_mi_tag>::type problem_set_by_id;
  const problem_set_by_id &set_id = moFEMProblems.get<BitProblemId_mi_tag>();
  problem_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_EDGEs(const Range& edges,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  rval = moab.add_entities(idm,edges.subset_by_type(MBEDGE)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_EDGEs(const Range& edges,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_EDGEs(edges,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_VERTICEs(const Range& vert,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  rval = moab.add_entities(idm,vert.subset_by_type(MBVERTEX)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_VERTICEs(const Range& vert,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_VERTICEs(vert,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_TRIs(const Range& tris,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  rval = moab.add_entities(idm,tris.subset_by_type(MBTRI)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_TRIs(const Range& tris,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_TRIs(tris,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const BitFEId id,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  EntityHandle idm = no_handle;
  try {
    idm = get_meshset_by_BitFEId(id);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,recursive); CHKERR_PETSC(rval);
  rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_TETs(const Range& tets,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  rval = moab.add_entities(idm,tets.subset_by_type(MBTET)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_TETs(const Range& tets,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_TETs(tets,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const string &name,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_TETs(meshset,get_BitFEId(name),recursive);  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit_ref,const string &name,EntityType type,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  typedef RefMoFEMElement_multiIndex::index<EntType_mi_tag>::type refMoabFE_by_type;
  refMoabFE_by_type &ref_MoFEMFiniteElement = refinedMoFemElements.get<EntType_mi_tag>();
  refMoabFE_by_type::iterator miit = ref_MoFEMFiniteElement.lower_bound(type);
  refMoabFE_by_type::iterator hi_miit = ref_MoFEMFiniteElement.upper_bound(type);
  if(verb > 1) {
    PetscPrintf(PETSC_COMM_WORLD,"nb. ref elements in database %d\n",distance(miit,hi_miit));
  }
  int nb_add_FEs = 0;
  for(;miit!=hi_miit;miit++) {
    if((miit->get_BitRefLevel()&bit_ref)!=bit_ref) continue;
    EntityHandle ent = miit->get_ref_ent();
    rval = moab.add_entities(idm,&ent,1); CHKERR_PETSC(rval);
    nb_add_FEs++;
  }
  if(verb > 0) {
    ostringstream ss;
    ss << "Add Nb. FEs " << nb_add_FEs << " form BitRef " << bit_ref << endl;
    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const string& name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  rval = moab.add_entities(idm,&meshset,1); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_MESHSETs(const EntityHandle meshset,const string& name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_meshset_by_BitFEId(id);
  Range meshsets;
  rval = moab.get_entities_by_type(meshset,MBENTITYSET,meshsets,false); CHKERR_PETSC(rval);
  rval = moab.add_entities(idm,meshsets); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::modify_problem_add_finite_element(const string &name_problem,const string &MoFEMFiniteElement_name) {
  PetscFunctionBegin;
  try {
    typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
    moFEMProblems_by_name& set = moFEMProblems.get<MoFEMProblem_mi_tag>();
    moFEMProblems_by_name::iterator miit = set.find(name_problem);
    ostringstream ss;
    ss << name_problem;
    if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,1,"this problem <%s> is not there",ss.str().c_str());
    BitFEId f_id = get_BitFEId(MoFEMFiniteElement_name);
    bool success = set.modify(miit,problem_MoFEMFiniteElement_change_bit_add(f_id));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::modify_problem_ref_level_add_bit(const string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name& set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator miit = set.find(name_problem);
  ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,1,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,problem_change_ref_level_bit_add(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::build_finite_elements(const EntMoFEMFiniteElement &EntFe,int verb) {
  PetscFunctionBegin;
  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_by_id;
  typedef MoFEMField_multiIndex::index<Meshset_mi_tag>::type field_by_meshset;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ent_by_ent;
  typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent>::type dof_set_type;
  field_by_id &moabFields_by_id = moabFields.get<BitFieldId_mi_tag>();
  field_by_meshset &moabFields_by_meshset = moabFields.get<Meshset_mi_tag>();
  dof_set_type& dof_set = dofsMoabField.get<Composite_Name_And_Ent>();
  EntityHandle fe_ent = EntFe.get_ent();
  pair<EntMoFEMFiniteElement_multiIndex::iterator,bool> p = finiteElementsMoFEMEnts.insert(EntFe);
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
    field_by_id::iterator miit = moabFields_by_id.find(BitFieldId().set(ii));
    if(miit==moabFields_by_id.end()) SETERRQ(PETSC_COMM_SELF,1,"data incosistency");
    //get field (ii) space
    FieldSpace space = miit->get_space();
    //resolve antities on element
    switch (moab.type_from_handle(fe_ent)) {
      case MBVERTEX:
	switch (space) {
	  case H1: 
	    adj_ents.insert(fe_ent);
	    break;
      	  default:
  	   SETERRQ(PETSC_COMM_SELF,1,"this fild is not implemented for TET finite element");
	}
	break;
      case MBEDGE:
	switch (space) {
	  case H1: if(nodes.empty()) moab.get_connectivity(&fe_ent,1,nodes,true);
	    adj_ents.insert(nodes.begin(),nodes.end());
	    adj_ents.insert(fe_ent);
	    break;
      	  default:
  	   SETERRQ(PETSC_COMM_SELF,1,"this fild is not implemented for TET finite element");
	}
	break;
      case MBTRI: 
	switch (space) {
	  case H1: if(nodes.empty()) moab.get_connectivity(&fe_ent,1,nodes,true);
	    adj_ents.insert(nodes.begin(),nodes.end());
	    for(Range::iterator eeit = edges.begin();eeit!=edges.end();eeit++) p.first->get_side_number_ptr(moab,*eeit);
	    adj_ents.insert(faces.begin(),faces.end());
	    for(Range::iterator fit = faces.begin();fit!=faces.end();fit++) p.first->get_side_number_ptr(moab,*fit);
	    adj_ents.insert(fe_ent);
	    break;
      	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"this fild is not implemented for TET finite element");
	}
	break;
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
	      SideNumber *side_ptr = EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,edge);
	      if(side_ptr->side_number!=ee) SETERRQ1(PETSC_COMM_SELF,1,"data insonsitenct for edge %d",ee);
	      rval = moab.side_element(prism,1,6+ee,edge); CHKERR_PETSC(rval);
	      side_ptr = EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,edge);
	      if(side_ptr->side_number!=ee+6) {
		if(side_ptr->side_number!=ee) {
		  SETERRQ1(PETSC_COMM_SELF,1,"data insonsitenct for edge %d",ee);
		} else {
		  side_ptr->brother_side_number = ee+6;
		}
	      }
	    }
	    int nn = 0;
	    for(;nn<3;nn++) {
	      EntityHandle node;
	      rval = moab.side_element(prism,0,nn,node); CHKERR_PETSC(rval);
	      SideNumber *side_ptr = EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,node);
	      if(side_ptr->side_number!=nn) SETERRQ1(PETSC_COMM_SELF,1,"data insonsitenct for node %d",nn);
	      rval = moab.side_element(prism,0,nn+3,node); CHKERR_PETSC(rval);
	      side_ptr = EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,node);
	      if(side_ptr->side_number!=nn+3) {
		if(side_ptr->side_number!=nn) {
		  SETERRQ1(PETSC_COMM_SELF,1,"data insonsitenct for node %d",nn);
		} else {
		  side_ptr->brother_side_number = nn+3; 
		}
	      }
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
	    case NoField:
	      if(moab.type_from_handle(*eit_eit)==MBENTITYSET) {
		//if field (ii) has space NoField only add dofs which associated with the meshsets
		if(moabFields_by_meshset.find(*eit_eit)!=moabFields_by_meshset.end()) {
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
	  SETERRQ(PETSC_COMM_SELF,1,"this finite element type is not implemented");
    }
    Range::iterator eit2 = adj_ents.begin();
    for(;eit2!=adj_ents.end();eit2++) {
      ref_ent_by_ent::iterator ref_ent_miit = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*eit2);
      if(ref_ent_miit==refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"ref ent not in database"); 
      const BitRefLevel& bit_ref_ent = ref_ent_miit->get_BitRefLevel();
      if(!(bit_ref_MoFEMFiniteElement&bit_ref_ent).any()) {
	ostringstream ss;
	ss << "top tip: check if you seed mesh with the elements for bit ref level1" << endl;
	ss << "inconsitency in database entity" << " type " << moab.type_from_handle(*eit2) << " bits ENT " << bit_ref_ent << endl;
	ss << "inconsitency in database entity" << " type " << moab.type_from_handle(p.first->get_ent()) << " bits FE  " << bit_ref_MoFEMFiniteElement << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      dof_set_type::iterator ents_miit2 = dof_set.lower_bound(boost::make_tuple(miit->get_name_ref(),ref_ent_miit->get_ref_ent()));
      dof_set_type::iterator ents_hi_miit2 = dof_set.upper_bound(boost::make_tuple(miit->get_name_ref(),ref_ent_miit->get_ref_ent()));
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
    	  case Row: success = finiteElementsMoFEMEnts.modify(p.first,EntMoFEMFiniteElement_row_dofs_change(moab,MoFEMFiniteElement_dof_uid_view[ss])); break;
      	  case Col: success = finiteElementsMoFEMEnts.modify(p.first,EntMoFEMFiniteElement_col_dofs_change(moab,MoFEMFiniteElement_dof_uid_view[ss])); break;
    	  case Data: success = finiteElementsMoFEMEnts.modify(p.first,EntMoFEMFiniteElement_data_dofs_change(moab,MoFEMFiniteElement_dof_uid_view[ss])); break;
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
  ierr = get_MoFEMFiniteElement_dof_uid_view(dofsMoabField,data_view,Interface::UNION,tag_data_uids_data,tag_data_uids_size); CHKERRQ(ierr);
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
    ierr = p.first->get_MoFEMFiniteElement_row_dof_uid_view(dofsMoabField,MoFEMFiniteElement_row_dof_uid_view); CHKERRQ(ierr);
    DofMoFEMEntity_multiIndex_uid_view::iterator miit_row = MoFEMFiniteElement_row_dof_uid_view.begin();
    ss << "rows dofs" << endl;
    for(;miit_row!=MoFEMFiniteElement_row_dof_uid_view.end();miit_row++) ss << **miit_row << endl;
    //cols
    DofMoFEMEntity_multiIndex_uid_view MoFEMFiniteElement_col_dof_uid_view;
    ierr = p.first->get_MoFEMFiniteElement_col_dof_uid_view(dofsMoabField,MoFEMFiniteElement_col_dof_uid_view); CHKERRQ(ierr);
    DofMoFEMEntity_multiIndex_uid_view::iterator miit_col = MoFEMFiniteElement_col_dof_uid_view.begin();
    ss << "cols dofs" << endl;
    for(;miit_col!=MoFEMFiniteElement_col_dof_uid_view.end();miit_col++) ss << **miit_col << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::build_finite_elements(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMElement_multiIndex::index<MoABEnt_mi_tag>::type ref_MoFEMFiniteElement_by_ent;
  MoFEMFiniteElement_multiIndex::iterator MoFEMFiniteElement_miit = finiteElements.begin();
  // loop Finite Elements
  for(;MoFEMFiniteElement_miit!=finiteElements.end();MoFEMFiniteElement_miit++) {
    if(verbose>0) PetscPrintf(PETSC_COMM_WORLD,"Build Finite Elements %s\n",MoFEMFiniteElement_miit->get_name().c_str());
    //get finite element meshset
    EntityHandle meshset = get_meshset_by_BitFEId(MoFEMFiniteElement_miit->get_id());
    // get entities from finite element meshset // if meshset 
    Range MoFEMFiniteElement_ents;
    rval = moab.get_entities_by_handle(meshset,MoFEMFiniteElement_ents,false); CHKERR_PETSC(rval);
    //loop meshset Ents and add finite elements
    Range::iterator eit = MoFEMFiniteElement_ents.begin();
    for(;eit!=MoFEMFiniteElement_ents.end();eit++) {
      // check if is in refinedMoFemElements database
      ref_MoFEMFiniteElement_by_ent::iterator ref_MoFEMFiniteElement_miit = refinedMoFemElements.get<MoABEnt_mi_tag>().find(*eit); /* iterator is a wrapper*/
      if(ref_MoFEMFiniteElement_miit == refinedMoFemElements.get<MoABEnt_mi_tag>().end()) {
	ostringstream ss;
	ss << "ref MoFEMFiniteElement not in database ent = " << *eit;
	ss << " type " << moab.type_from_handle(*eit);
	ss << " " << *MoFEMFiniteElement_miit;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      ierr = build_finite_elements(EntMoFEMFiniteElement(moab,ref_MoFEMFiniteElement_miit->get_RefMoFEMElement(),&*MoFEMFiniteElement_miit),verb); CHKERRQ(ierr);
    }
  }
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Nb. FEs %u\n",finiteElementsMoFEMEnts.size());
    typedef EntMoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type MoFEMFiniteElement_by_id;
    MoFEMFiniteElement_by_id &MoFEMFiniteElements = finiteElementsMoFEMEnts.get<BitFEId_mi_tag>();  
    MoFEMFiniteElement_multiIndex::iterator id_MoFEMFiniteElement = finiteElements.begin();
    for(;id_MoFEMFiniteElement!=finiteElements.end();id_MoFEMFiniteElement++) {
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
PetscErrorCode FieldCore::build_adjacencies(const BitRefLevel bit) {
  PetscFunctionBegin;
  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,1,"field not build");
  if(!(*build_MoFEM)&(1<<1)) SETERRQ(PETSC_COMM_SELF,1,"fe not build");
  typedef MoFEMEntity_multiIndex::index<Unique_mi_tag>::type ents_by_uid;
  EntMoFEMFiniteElement_multiIndex::iterator fit = finiteElementsMoFEMEnts.begin();
  for(;fit!=finiteElementsMoFEMEnts.end();fit++) {
    if(!(fit->get_BitRefLevel()&bit).any()) continue;
    int size_row = fit->tag_row_uids_size/sizeof(UId);
    const UId *uids_row = (UId*)fit->tag_row_uids_data;
    int ii = 0;
    UId uid = 0;
    for(;ii<size_row;ii++) {
      if( uid == (uids_row[ii] >> 8 )) continue;
      uid = uids_row[ii];
      uid = uid >> 8; //look to DofMoFEMEntity::get_unique_id_calculate and MoFEMEntity::get_unique_id_calculate() <- uid is shifted by 8 bits
      ents_by_uid::iterator miit = entsMoabField.get<Unique_mi_tag>().find(uid);
      assert(dofsMoabField.get<Unique_mi_tag>().find(uids_row[ii])!=dofsMoabField.get<Unique_mi_tag>().end());
      assert(dofsMoabField.get<Unique_mi_tag>().find(uids_row[ii])->get_MoFEMEntity_ptr()->get_unique_id()==uid);
      if(miit ==entsMoabField.get<Unique_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      pair<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator,bool> p = entFEAdjacencies.insert(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap(&*miit,&*fit));
      bool success = entFEAdjacencies.modify(p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_by_what(by_row));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    int size_col = fit->tag_col_uids_size/sizeof(UId);
    const UId *uids_col = (UId*)fit->tag_col_uids_data;
    for(ii = 0,uid = 0;ii<size_col;ii++) {
      if( uid == (uids_col[ii] >> 8 )) continue;
      uid = uids_col[ii];
      uid = uid >> 8; //look to DofMoFEMEntity::get_unique_id_calculate and MoFEMEntity::get_unique_id_calculate() <- uid is shifted by 8 bits
      assert(dofsMoabField.get<Unique_mi_tag>().find(uids_col[ii])!=dofsMoabField.get<Unique_mi_tag>().end());
      assert(dofsMoabField.get<Unique_mi_tag>().find(uids_col[ii])->get_MoFEMEntity_ptr()->get_unique_id()==uid);
      ents_by_uid::iterator miit = entsMoabField.get<Unique_mi_tag>().find(uid);
      if(miit ==entsMoabField.get<Unique_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      pair<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator,bool> p = entFEAdjacencies.insert(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap(&*miit,&*fit));
      bool success = entFEAdjacencies.modify(p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_by_what(by_col));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    int size_data = fit->tag_data_uids_size/sizeof(UId);
    const UId *uids_data = (UId*)fit->tag_data_uids_data;
    for(ii = 0,uid = 0;ii<size_data;ii++) {
      if( uid == (uids_data[ii] >> 8 )) continue;
      uid = uids_data[ii];
      uid = uid >> 8; //look to DofMoFEMEntity::get_unique_id_calculate and MoFEMEntity::get_unique_id_calculate() <- uid is shifted by 8 bits
      assert(dofsMoabField.get<Unique_mi_tag>().find(uids_data[ii])!=dofsMoabField.get<Unique_mi_tag>().end());
      assert(dofsMoabField.get<Unique_mi_tag>().find(uids_data[ii])->get_MoFEMEntity_ptr()->get_unique_id()==uid);
      ents_by_uid::iterator miit = entsMoabField.get<Unique_mi_tag>().find(uid);
      if(miit == entsMoabField.get<Unique_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      pair<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator,bool> p = entFEAdjacencies.insert(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap(&*miit,&*fit));
      bool success = entFEAdjacencies.modify(p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_by_what(by_data));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
  }
  if(verbose>1) {
    list_adjacencies();
  }
  if(verbose>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Nb. entFEAdjacencies %u\n",entFEAdjacencies.size());
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Nb. entFEAdjacencies %u\n",entFEAdjacencies.size());
    //PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  *build_MoFEM |= 1<<2;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::list_adjacencies() const {
  PetscFunctionBegin;
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator miit = entFEAdjacencies.begin();
  for(;miit!=entFEAdjacencies.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::build_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"entFEAdjacencies not build");
  MoFEMProblem_multiIndex::iterator p_miit = moFEMProblems.begin();
  for(;p_miit!=moFEMProblems.end();p_miit++) {
    if(p_miit->get_BitRefLevel().none()) SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinment level not set",p_miit->get_name().c_str());
    //miit2 iterator for finite elements
    EntMoFEMFiniteElement_multiIndex::iterator miit2 = finiteElementsMoFEMEnts.begin();
    EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finiteElementsMoFEMEnts.end();
    DofMoFEMEntity_multiIndex_uid_view dofs_rows;
    DofMoFEMEntity_multiIndex_uid_view dofs_cols;;
    EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
    for(;miit3!=hi_miit2;miit3++) {
      if((miit3->get_id()&p_miit->get_BitFEId()).any()) {
	if((miit3->get_BitRefLevel()&p_miit->get_BitRefLevel())==p_miit->get_BitRefLevel()) {
	  ierr = miit3->get_MoFEMFiniteElement_row_dof_uid_view(dofsMoabField,dofs_rows); CHKERRQ(ierr);
	  ierr = miit3->get_MoFEMFiniteElement_col_dof_uid_view(dofsMoabField,dofs_cols); CHKERRQ(ierr);
	}
      }
    }
    //zero rows
    bool success = moFEMProblems.modify(p_miit,problem_zero_nb_rows_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    //zero cols
    success = moFEMProblems.modify(p_miit,problem_zero_nb_cols_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    //add dofs for rows
    DofMoFEMEntity_multiIndex_uid_view::iterator miit4 = dofs_rows.begin();
    for(;miit4!=dofs_rows.end();miit4++) {
      if(!(*miit4)->active) continue;
      success = moFEMProblems.modify(p_miit,problem_row_change(&**miit4));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    //add dofs for cols
    DofMoFEMEntity_multiIndex_uid_view::iterator miit5 = dofs_cols.begin();
    for(;miit5!=dofs_cols.end();miit5++) {
      if(!(*miit5)->active) continue;
      success = moFEMProblems.modify(p_miit,problem_col_change(&**miit5));
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
PetscErrorCode FieldCore::simple_partition_problem(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,1,"moFEMProblems not build");
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  // find p_miit
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > is not found (top tip: check spelling)",name.c_str());
  typedef boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
  NumeredDofMoFEMEntitys_by_idx &dofs_row_by_idx = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_rows.get<Idx_mi_tag>());
  NumeredDofMoFEMEntitys_by_idx &dofs_col_by_idx = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_cols.get<Idx_mi_tag>());
  boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Idx_mi_tag>::type::iterator miit_row,hi_miit_row;
  boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Idx_mi_tag>::type::iterator miit_col,hi_miit_col;
  DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
  nb_row_local_dofs = 0;
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  nb_row_ghost_dofs = 0;
  DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
  nb_col_local_dofs = 0;
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  nb_col_ghost_dofs = 0;
  for(unsigned int part = 0;part<pcomm->size();part++) {
    DofIdx nb_dofs_row = dofs_row_by_idx.size();
    assert(p_miit->get_nb_dofs_row()==nb_dofs_row);
    DofIdx nb_dofs_row_on_proc = (DofIdx)ceil(nb_dofs_row/pcomm->size());
    DofIdx lower_dof_row = nb_dofs_row_on_proc*part;
    miit_row = dofs_row_by_idx.lower_bound(lower_dof_row);
    DofIdx upper_dof_row = part==pcomm->size()-1 ? nb_dofs_row-1 : nb_dofs_row_on_proc*(part+1)-1;
    hi_miit_row = dofs_row_by_idx.upper_bound(upper_dof_row);
    // loop rows
    for(;miit_row!=hi_miit_row;miit_row++) {
      bool success = dofs_row_by_idx.modify(miit_row,NumeredDofMoFEMEntity_part_change(part,miit_row->dof_idx));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      if(miit_row->part == pcomm->rank()) {
	success = dofs_row_by_idx.modify(miit_row,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }
    }
    DofIdx nb_dofs_col = dofs_col_by_idx.size();
    assert(p_miit->get_nb_dofs_col()==nb_dofs_col);
    DofIdx nb_dofs_col_on_proc = (DofIdx)ceil(nb_dofs_col/pcomm->size());
    DofIdx lower_dof_col = nb_dofs_col_on_proc*part;
    miit_col = dofs_col_by_idx.lower_bound(lower_dof_col);
    DofIdx upper_dof_col = part==pcomm->size()-1 ? nb_dofs_col-1 : nb_dofs_col_on_proc*(part+1)-1;
    hi_miit_col = dofs_col_by_idx.upper_bound(upper_dof_col);
    // loop cols
    for(;miit_col!=hi_miit_col;miit_col++) {
      bool success = dofs_col_by_idx.modify(miit_col,NumeredDofMoFEMEntity_part_change(part,miit_col->dof_idx));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      if(miit_col->part == pcomm->rank()) {
	success = dofs_col_by_idx.modify(miit_col,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }
    }
  }
  if(verbose>0) {
    ostringstream ss;
    ss << "simple_partition_problem: rank = " << pcomm->rank() << " FEs row ghost dofs "<< *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_row() << " nb global row dofs " << p_miit->get_nb_dofs_row() << endl;
    ss << "simple_partition_problem: rank = " << pcomm->rank() << " FEs col ghost dofs " << *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_col() << " nb global col dofs " << p_miit->get_nb_dofs_col() << endl;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  *build_MoFEM |= 1<<4;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::partition_problem(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,1,"moFEMProblems not build");
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef NumeredDofMoFEMEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",name.c_str());
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Partition problem %s\n",p_miit->get_name().c_str());
  }
  DofIdx nb_dofs_row = p_miit->get_nb_dofs_row();
  Mat Adj;
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"\tcreate Adj matrix\n");
  }
  ierr = partition_create_Mat<Idx_mi_tag>(name,&Adj,MATMPIADJ,true,verb); CHKERRQ(ierr);
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
  if(debug>0) {
    NumeredDofMoFEMEntitys_by_idx::iterator dit,hi_dit;
    dit = p_miit->numered_dofs_rows.get<Idx_mi_tag>().begin();
    hi_dit = p_miit->numered_dofs_rows.get<Idx_mi_tag>().end();
    for(;dit!=hi_dit;dit++) {
      if(dit->get_part()==pcomm->rank()) {
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "rank " << pcomm->rank() << " " << *dit;
	  SETERRQ1(PETSC_COMM_SELF,1,"local dof index for row not set\n %s",ss.str().c_str());
	}
      }
    }
    dit = p_miit->numered_dofs_cols.get<Idx_mi_tag>().begin();
    hi_dit = p_miit->numered_dofs_cols.get<Idx_mi_tag>().end();
    for(;dit!=hi_dit;dit++) {
      if(dit->get_part()==pcomm->rank()) {
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "rank " << pcomm->rank() << " " << *dit;
	  SETERRQ1(PETSC_COMM_SELF,1,"local dof index for col not set\n %s",ss.str().c_str());
	}
      }
    }
  }
  *build_MoFEM |= 1<<4;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::compose_problem(const string &name,const string &problem_for_rows,const string &problem_for_cols,int verb) {
  PetscFunctionBegin;
  ierr = compose_problem(name,problem_for_rows,false,problem_for_cols,false,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::compose_problem(const string &name,const string &problem_for_rows,bool copy_rows,const string &problem_for_cols,bool copy_cols,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,1,"moFEMProblems not build");
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofMoFEMEntitys_by_uid;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",name.c_str());
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Partition problem %s\n",p_miit->get_name().c_str());
  }
  //find p_miit_row
  moFEMProblems_by_name::iterator p_miit_row = moFEMProblems_set.find(problem_for_rows);
  if(p_miit_row==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",problem_for_rows.c_str());
  moFEMProblems_by_name::iterator p_miit_col = moFEMProblems_set.find(problem_for_cols);
  //find p_mit_col
  if(p_miit_col==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",problem_for_cols.c_str());
  const NumeredDofMoFEMEntity_multiIndex &dofs_col = p_miit_col->numered_dofs_cols;
  //do rows
  map<DofIdx,const NumeredDofMoFEMEntity*> rows_problem_map;
  const NumeredDofMoFEMEntity_multiIndex &dofs_row = p_miit_row->numered_dofs_rows;
  if(!copy_rows) {
  MoFEMEntity *MoFEMEntity_ptr = NULL;
  NumeredDofMoFEMEntity_multiIndex::iterator miit_row = dofs_row.begin();
  NumeredDofMoFEMEntity_multiIndex::iterator hi_miit_row = dofs_row.end();
  for(;miit_row!=hi_miit_row;miit_row++) {
    if( (MoFEMEntity_ptr == NULL) ? 1 : (MoFEMEntity_ptr->get_unique_id() != miit_row->field_ptr->field_ptr->get_unique_id()) ) {
      MoFEMEntity_ptr = const_cast<MoFEMEntity*>(miit_row->field_ptr->field_ptr);
      typedef MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type adj_by_ent;
      adj_by_ent::iterator adj_miit = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(MoFEMEntity_ptr->get_unique_id());
      adj_by_ent::iterator hi_adj_miit = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(MoFEMEntity_ptr->get_unique_id());
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
  if(rows_problem_map.size() != (unsigned int)p_miit->get_nb_dofs_row()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  } 
  //do cols
  map<DofIdx,const NumeredDofMoFEMEntity*> cols_problem_map;
  if(!copy_cols) {
  MoFEMEntity *MoFEMEntity_ptr = NULL;
  NumeredDofMoFEMEntity_multiIndex::iterator miit_col = dofs_col.begin();
  NumeredDofMoFEMEntity_multiIndex::iterator hi_miit_col = dofs_col.end();
  for(;miit_col!=hi_miit_col;miit_col++) {
    if( (MoFEMEntity_ptr == NULL) ? 1 : (MoFEMEntity_ptr->get_unique_id() != miit_col->field_ptr->field_ptr->get_unique_id()) ) {
      MoFEMEntity_ptr = const_cast<MoFEMEntity*>(miit_col->field_ptr->field_ptr);
      typedef MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type adj_by_ent;
      adj_by_ent::iterator adj_miit = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(MoFEMEntity_ptr->get_unique_id());
      adj_by_ent::iterator hi_adj_miit = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(MoFEMEntity_ptr->get_unique_id());
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
  if(cols_problem_map.size() != (unsigned int)p_miit->get_nb_dofs_col()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  // build indices 
  NumeredDofMoFEMEntitys_by_uid &dofs_row_by_uid = const_cast<NumeredDofMoFEMEntitys_by_uid&>(p_miit->numered_dofs_rows.get<Unique_mi_tag>());
  DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
  nb_row_local_dofs = 0;
  if(!copy_rows) {
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
    if(pr_dof->get_part() == pcomm->rank()) {
      success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
  }
  } else {
    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(p_miit_row,diit)) {
      bool success = moFEMProblems.modify(moFEMProblems.project<0>(p_miit),problem_row_change(diit->get_DofMoFEMEntity_ptr()));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      int part_number = diit->get_part();
      int petsc_global_dof = diit->get_petsc_gloabl_dof_idx();
      NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_row_by_uid.find(diit->get_unique_id());
      if(pr_dof == dofs_row_by_uid.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      if(pr_dof->get_part() == pcomm->rank()) {
	success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }
    }
  }
  //cols
  NumeredDofMoFEMEntitys_by_uid &dofs_col_by_uid = const_cast<NumeredDofMoFEMEntitys_by_uid&>(p_miit->numered_dofs_cols.get<Unique_mi_tag>());
  DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
  nb_col_local_dofs = 0;
  if(!copy_cols) {
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
    if(pr_dof->get_part() == pcomm->rank()) {
      success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
  }
  } else {
    for(_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(p_miit_col,diit)) {
      bool success = moFEMProblems.modify(moFEMProblems.project<0>(p_miit),problem_col_change(diit->get_DofMoFEMEntity_ptr()));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      int part_number = diit->get_part();
      int petsc_global_dof = diit->get_petsc_gloabl_dof_idx();
      NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_col_by_uid.find(diit->get_unique_id());
      success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      if(pr_dof->get_part() == pcomm->rank()) {
	success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }
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
  if(debug>0) {
    typedef NumeredDofMoFEMEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
    NumeredDofMoFEMEntitys_by_idx::iterator dit,hi_dit;
    dit = p_miit->numered_dofs_rows.get<Idx_mi_tag>().begin();
    hi_dit = p_miit->numered_dofs_rows.get<Idx_mi_tag>().end();
    for(;dit!=hi_dit;dit++) {
      if(dit->get_part()==pcomm->rank()) {
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "rank " << pcomm->rank() << " " << *dit;
	  SETERRQ1(PETSC_COMM_SELF,1,"local dof index for row not set\n %s",ss.str().c_str());
	}
      }
    }
    dit = p_miit->numered_dofs_cols.get<Idx_mi_tag>().begin();
    hi_dit = p_miit->numered_dofs_cols.get<Idx_mi_tag>().end();
    for(;dit!=hi_dit;dit++) {
      if(dit->get_part()==pcomm->rank()) {
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "rank " << pcomm->rank() << " " << *dit;
	  SETERRQ1(PETSC_COMM_SELF,1,"local dof index for col not set\n %s",ss.str().c_str());
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::partition_ghost_dofs(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,1,"moFEMProblems not build");
  if(!(*build_MoFEM&(1<<4))) SETERRQ(PETSC_COMM_SELF,1,"partitions moFEMProblems not build");
  if(!(*build_MoFEM&(1<<5))) SETERRQ(PETSC_COMM_SELF,1,"partitions finite elements not build");
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  //
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  //
  nb_row_ghost_dofs = 0;
  nb_col_ghost_dofs = 0;
  if(pcomm->size()>1) {
    NumeredDofMoFEMEntity_multiIndex_uid_view ghost_idx_col_view,ghost_idx_row_view;
    NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_Part_mi_tag>::type::iterator fe_it,hi_fe_it;
    fe_it = p_miit->numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().lower_bound(pcomm->rank());
    hi_fe_it = p_miit->numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().upper_bound(pcomm->rank());
    for(;fe_it!=hi_fe_it;fe_it++) {
      if(fe_it->rows_dofs.empty()) SETERRQ(PETSC_COMM_SELF,1,"no row dofs on this element why it is here?");
      if(fe_it->cols_dofs.empty()) SETERRQ(PETSC_COMM_SELF,1,"no col dofs on this element why it is here?");
      typedef FENumeredDofMoFEMEntity_multiIndex::iterator dof_it;
      dof_it rowdofit,hi_rowdofit;
      rowdofit = fe_it->rows_dofs.begin();
      hi_rowdofit = fe_it->rows_dofs.end();
      for(;rowdofit!=hi_rowdofit;rowdofit++) {
	if(rowdofit->get_part()==pcomm->rank()) continue; 
	ghost_idx_row_view.insert(rowdofit->get_NumeredDofMoFEMEntity_ptr());
      }
      dof_it coldofit,hi_coldofit;
      coldofit = fe_it->cols_dofs.begin();
      hi_coldofit = fe_it->cols_dofs.end();
      for(;coldofit!=hi_coldofit;coldofit++) {
	if(coldofit->get_part()==pcomm->rank()) continue;
	ghost_idx_col_view.insert(coldofit->get_NumeredDofMoFEMEntity_ptr());
      }
    }
    DofIdx *nb_ghost_dofs[2] = { &nb_col_ghost_dofs, &nb_row_ghost_dofs };
    DofIdx nb_local_dofs[2] = { *((DofIdx*)p_miit->tag_local_nbdof_data_col), *((DofIdx*)p_miit->tag_local_nbdof_data_row) };
    NumeredDofMoFEMEntity_multiIndex_uid_view *ghost_idx_view[2] = { &ghost_idx_col_view, &ghost_idx_row_view };
    typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofMoFEMEntitys_by_unique_id;
    NumeredDofMoFEMEntitys_by_unique_id *dof_by_uid_no_const[2] = {
      const_cast<NumeredDofMoFEMEntitys_by_unique_id*>(&p_miit->numered_dofs_cols.get<Unique_mi_tag>()),
      const_cast<NumeredDofMoFEMEntitys_by_unique_id*>(&p_miit->numered_dofs_rows.get<Unique_mi_tag>())    
    };
    for(int ss = 0;ss<2;ss++) {
      NumeredDofMoFEMEntity_multiIndex_uid_view::iterator ghost_idx_miit = ghost_idx_view[ss]->begin();
      for(;ghost_idx_miit!=ghost_idx_view[ss]->end();ghost_idx_miit++) {
        NumeredDofMoFEMEntitys_by_unique_id::iterator diit = dof_by_uid_no_const[ss]->find((*ghost_idx_miit)->get_unique_id());
        if(diit->petsc_local_dof_idx!=(DofIdx)-1) SETERRQ(PETSC_COMM_SELF,1,"inconsitent data, ghost dof already set");
        bool success = dof_by_uid_no_const[ss]->modify(diit,NumeredDofMoFEMEntity_local_idx_change(nb_local_dofs[ss]++));
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
  if(debug>0) {
    NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_Part_mi_tag>::type::iterator fit,hi_fit;
    fit = p_miit->numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().lower_bound(pcomm->rank());
    hi_fit = p_miit->numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().upper_bound(pcomm->rank());
    for(;fit!=hi_fit;fit++) {
      const UId* row_uids = fit->fe_ptr->tag_row_uids_data;
      int row_size = (fit->fe_ptr->tag_row_uids_size)/sizeof(UId);
      const UId* col_uids = fit->fe_ptr->tag_col_uids_data;
      int col_size = (fit->fe_ptr->tag_col_uids_size)/sizeof(UId);
      NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator dit;
      for(int dd = 0;dd<row_size;dd++) {
	dit = p_miit->numered_dofs_rows.get<Unique_mi_tag>().find(row_uids[dd]); 
	if(dit == p_miit->numered_dofs_rows.get<Unique_mi_tag>().end()) {
	  SETERRQ(PETSC_COMM_SELF,1,"uid not in row databse");
	} 
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "proc " << pcomm->rank() << " : " << *dit;
	  if(dit->get_part()==pcomm->rank()) {
	    SETERRQ1(PETSC_COMM_SELF,1,"local dof for row not set\n%s",ss.str().c_str());
	  } else {
	    SETERRQ1(PETSC_COMM_SELF,1,"ghost dof for row not set\n%s",ss.str().c_str());
	  }
	}
      }
      for(int dd = 0;dd<col_size;dd++) {
	dit = p_miit->numered_dofs_cols.get<Unique_mi_tag>().find(col_uids[dd]); 
	if(dit == p_miit->numered_dofs_cols.get<Unique_mi_tag>().end()) {
	  SETERRQ(PETSC_COMM_SELF,1,"uid not in col databse");
	} 
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "proc " << pcomm->rank() << " : " << *dit;
	  if(dit->get_part()==pcomm->rank()) {
	    SETERRQ1(PETSC_COMM_SELF,1,"local dof for col not set\n%s",ss.str().c_str());
	  } else {
	    SETERRQ1(PETSC_COMM_SELF,1,"ghost dof for col not set\n%s",ss.str().c_str());
	  }
	}
      }
    }
  }
  *build_MoFEM |= 1<<6;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::partition_finite_elements(const string &name,bool do_skip,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,1,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,1,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,1,"partitions not build");
  if(!(*build_MoFEM&(1<<4))) SETERRQ(PETSC_COMM_SELF,1,"partitions moFEMProblems not build");
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str());
  NumeredMoFEMFiniteElement_multiIndex& numeredFiniteElements = const_cast<NumeredMoFEMFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
  //MoFEMFiniteElement set
  EntMoFEMFiniteElement_multiIndex::iterator miit2 = finiteElementsMoFEMEnts.begin();
  EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finiteElementsMoFEMEnts.end();
  EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
  for(;miit3!=hi_miit2;miit3++) {
    if((miit3->get_BitRefLevel()&p_miit->get_BitRefLevel()).none()) continue;
    if((miit3->get_id()&p_miit->get_BitFEId()).any()) {
      NumeredDofMoFEMEntity_multiIndex_uid_view rows_view,cols_view;
      //rows_view
      const void* tag_row_uids_data = miit3->tag_row_uids_data;
      const int tag_row_uids_size = miit3->tag_row_uids_size;
      ierr = get_MoFEMFiniteElement_dof_uid_view(p_miit->numered_dofs_rows,rows_view,Interface::UNION,tag_row_uids_data,tag_row_uids_size); CHKERRQ(ierr);
      if(rows_view.empty()) continue;
      //cols_vies
      const void* tag_col_uids_data = miit3->tag_col_uids_data;
      const int tag_col_uids_size = miit3->tag_col_uids_size;
      ierr = get_MoFEMFiniteElement_dof_uid_view(p_miit->numered_dofs_cols,cols_view,Interface::UNION,tag_col_uids_data,tag_col_uids_size); CHKERRQ(ierr);
      if(cols_view.empty()) continue;
      pair<NumeredMoFEMFiniteElement_multiIndex::iterator,bool> p = numeredFiniteElements.insert(NumeredMoFEMFiniteElement(&*miit3));
      NumeredMoFEMFiniteElement &problem_MoFEMFiniteElement = const_cast<NumeredMoFEMFiniteElement&>(*p.first);
      if(!p.second) {
	problem_MoFEMFiniteElement.rows_dofs.clear();
	problem_MoFEMFiniteElement.cols_dofs.clear();
      }
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
	parts[(*viit_rows)->part]++;
      }
      vector<int>::iterator pos = max_element(parts.begin(),parts.end());
      unsigned int max_part = std::distance(parts.begin(),pos);
      bool success = numeredFiniteElements.modify(p.first,NumeredMoFEMFiniteElement_change_part(max_part));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      if(do_skip) if(max_part!=pcomm->rank()) continue; 
      //cols
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
	PetscSynchronizedFlush(PETSC_COMM_WORLD); 
      }
    }
  }
  if(verb>0) {
    typedef NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_Part_mi_tag>::type NumeredMoFEMFiniteElement_multiIndex_by_part;
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator MoFEMFiniteElement_miit = numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().lower_bound(pcomm->rank());
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator hi_MoMoFEMFiniteElement_miitFEMFE_miit = numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().upper_bound(pcomm->rank());
    int count = std::distance(MoFEMFiniteElement_miit,hi_MoMoFEMFiniteElement_miitFEMFE_miit);
    ostringstream ss;
    ss << *p_miit;
    ss << " Nb. elems " << count << " on proc " << pcomm->rank() << endl;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  if(debug>0) {
    typedef NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_Part_mi_tag>::type NumeredMoFEMFiniteElement_multiIndex_by_part;
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator fit = numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().lower_bound(pcomm->rank());
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator hi_fit = numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().upper_bound(pcomm->rank());
    for(;fit!=hi_fit;fit++) {
      const UId* row_uids = fit->fe_ptr->tag_row_uids_data;
      int row_size = (fit->fe_ptr->tag_row_uids_size)/sizeof(UId);
      const UId* col_uids = fit->fe_ptr->tag_col_uids_data;
      int col_size = (fit->fe_ptr->tag_col_uids_size)/sizeof(UId);
      NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator dit;
      for(int dd = 0;dd<row_size;dd++) {
	dit = p_miit->numered_dofs_rows.get<Unique_mi_tag>().find(row_uids[dd]); 
	if(dit == p_miit->numered_dofs_rows.get<Unique_mi_tag>().end()) {
	  SETERRQ(PETSC_COMM_SELF,1,"uid not in row databse");
	} 
	if(dit->get_petsc_local_dof_idx()<0) {
	  if(dit->get_part()==pcomm->rank()) {
	    ostringstream ss;
	    ss << "proc " << pcomm->rank() << " : " << *dit;
	    SETERRQ1(PETSC_COMM_SELF,1,"local dof for row not set\n%s",ss.str().c_str());
	  } 
	}
      }
      for(int dd = 0;dd<col_size;dd++) {
	dit = p_miit->numered_dofs_cols.get<Unique_mi_tag>().find(col_uids[dd]); 
	if(dit == p_miit->numered_dofs_cols.get<Unique_mi_tag>().end()) {
	  SETERRQ(PETSC_COMM_SELF,1,"uid not in col databse");
	} 
	if(dit->get_petsc_local_dof_idx()<0) {
	  if(dit->get_part()==pcomm->rank()) {
	    ostringstream ss;
	    ss << "proc " << pcomm->rank() << " : " << *dit;
	    SETERRQ1(PETSC_COMM_SELF,1,"local dof for col not set\n%s",ss.str().c_str());
	  } 
	}
      }
    }
  }
  *build_MoFEM |= 1<<5;  
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::seed_finite_elements(const EntityHandle meshset,int verb) {
  PetscFunctionBegin;
  Range entities;
  ierr = moab.get_entities_by_handle(meshset,entities,true); CHKERRQ(ierr);
  ierr = seed_finite_elements(entities,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::seed_finite_elements(const Range &entities,int verb) {
  PetscFunctionBegin;
  for(Range::iterator eit = entities.begin();eit!=entities.end();eit++) {
    RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator 
      eiit = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*eit);
    if(eiit == refinedMoFemEntities.get<MoABEnt_mi_tag>().end())  SETERRQ(PETSC_COMM_SELF,1,"entity is not in database");
    if(eiit->get_BitRefLevel().none()) continue;
    pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
    switch (eiit->get_ent_type()) {
      case MBVERTEX: 
	p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_VERTEX(moab,&*eiit)));	
	break;
      case MBEDGE: 
	p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_EDGE(moab,&*eiit)));	
	break;
      case MBTRI: 
	p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TRI(moab,&*eiit)));	
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::seed_ref_level_2D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin; 
  if(verb==-1) verb = verbose;
  try {
    Range ents2d;
    rval = moab.get_entities_by_type(meshset,MBTRI,ents2d,true); CHKERR_PETSC(rval);
    Range ents;
    rval = moab.get_adjacencies(ents2d,1,true,ents,Interface::UNION); CHKERR_PETSC(rval);
    if(verb > 1) {
      PetscPrintf(PETSC_COMM_WORLD,"nb. 2d entities for seed %d\n",ents2d.size());
    }
    Range::iterator tit = ents2d.begin();
    for(;tit!=ents2d.end();tit++) {
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*tit));
      if(debug > 0) {
	ierr = test_moab(moab,*tit); CHKERRQ(ierr);
      }
      if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
        bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }
      pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      switch (p_ent.first->get_ent_type()) {
        case MBTRI: 
	 p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TRI(moab,&*p_ent.first)));	
	  assert(p_MoFEMFiniteElement.first->get_BitRefEdges_ulong()!=-1);
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
    for(int dd = 0;dd<1;dd++) {
      rval = moab.get_entities_by_dimension(meshset,dd,ents); CHKERR_PETSC(rval);
      Range::iterator eit = ents.begin();
      for(;eit!=ents.end();eit++) {
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*eit));
        if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
	  bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
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
PetscErrorCode FieldCore::seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin; 
  if(verb==-1) verb = verbose;
  try {
    Range ents3d;
    rval = moab.get_entities_by_type(meshset,MBTET,ents3d,false); CHKERR_PETSC(rval);
    Range ents;
    rval = moab.get_adjacencies(ents3d,2,true,ents,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(ents3d,1,true,ents,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.get_entities_by_type(meshset,MBPRISM,ents3d,false); CHKERR_PETSC(rval);
    if(verb > 1) {
      PetscPrintf(PETSC_COMM_WORLD,"nb. 3d entities for seed %d\n",ents3d.size());
    }
    Range::iterator tit = ents3d.begin();
    for(;tit!=ents3d.end();tit++) {
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*tit));
      if(debug > 0) {
	ierr = test_moab(moab,*tit); CHKERRQ(ierr);
      }
      if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
        bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }
      pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      switch (p_ent.first->get_ent_type()) {
        case MBTET: 
	 p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_ent.first)));	
	 break;
	case MBPRISM:
	  p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
	  break;
        case MBENTITYSET:
	  p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_MESHSET(moab,&*p_ent.first)));
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
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*eit));
        if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
	  bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
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
PetscErrorCode FieldCore::seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit) {
  PetscFunctionBegin;
  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,meshset));
  if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
    refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
  }
  ptrWrapperRefMoFEMElement pack_fe(new RefMoFEMElement_MESHSET(moab,&*p_ent.first));
  pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement = refinedMoFemElements.insert(pack_fe);
  if(verbose > 0) {
    ostringstream ss;
    ss << "add meshset as ref_ent " << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement()) << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_verices_in_the_middel_of_edges(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
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
  ierr = add_verices_in_the_middel_of_edges(edges,bit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_verices_in_the_middel_of_edges(const Range &_edges,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  Range edges = _edges;
  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityType_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedMoFemEntities.get<Composite_EntityType_And_ParentEntityType_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) ref_parent_ents_view.insert(&*miit);
  // refine edges on the other side of the prism
  typedef BasicMoFEMEntityAdjacenctMap_multiIndex::index<MoABEnt_mi_tag2>::type BasicMoFEMEntityAdjacenctMap_by_adj;
  BasicMoFEMEntityAdjacenctMap_by_adj &basicEntAdjacencies_by_adj = basicEntAdjacencies.get<MoABEnt_mi_tag2>();
  Range::iterator eit = edges.begin();
  for(;eit!=edges.end();eit++) {
    BasicMoFEMEntityAdjacenctMap_by_adj::iterator adj_miit = basicEntAdjacencies_by_adj.find(*eit);
    if(adj_miit==basicEntAdjacencies_by_adj.end()) continue;
    EntityHandle prism = adj_miit->ent;
    RefMoFEMElement_multiIndex::iterator miit2 = refinedMoFemElements.get<MoABEnt_mi_tag>().find(prism);
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
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,node));
      if(!p_ent.second) SETERRQ(PETSC_COMM_SELF,1,"this entity is there");
      if(verbose>2) {
	ostringstream ss;
	ss << *(p_ent.first) << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    } else {
      const EntityHandle node = (*miit_view)->get_ref_ent();
      bool success = refinedMoFemEntities.modify(refinedMoFemEntities.get<MoABEnt_mi_tag>().find(node),RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsitency in data");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface) {
  PetscFunctionBegin;
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,false); CHKERR_PETSC(rval);
  ierr = refine_TET(tets,bit,respect_interface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_TET(const Range &_tets,const BitRefLevel &bit,const bool respect_interface) {
  PetscFunctionBegin;
  //FIXME: refinment is based on entity handlers, should work on global ids of nodes, this will allow parallelize agortihm in the future
  PetscFunctionBegin;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ents_by_ent;
  ref_ents_by_ent &ref_ents_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>();
  // find all verices which parent is edge
  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityType_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedMoFemEntities.get<Composite_EntityType_And_ParentEntityType_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) ref_parent_ents_view.insert(&*miit);
  typedef RefMoFEMElement_multiIndex::index<MoABEnt_mi_tag>::type ref_MoFEMFiniteElement_by_ent;
  ref_MoFEMFiniteElement_by_ent &ref_MoFEMFiniteElement = refinedMoFemElements.get<MoABEnt_mi_tag>();
  typedef RefMoFEMElement_multiIndex::index<Composite_of_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type ref_ent_by_composite;
  ref_ent_by_composite &by_composite = refinedMoFemElements.get<Composite_of_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
  // find oposite intrface nodes
  typedef BasicMoFEMEntityAdjacenctMap_multiIndex::index<EntType_mi_tag>::type AdjPrism_by_type;
  AdjPrism_by_type::iterator face_prism_miit = basicEntAdjacencies.get<EntType_mi_tag>().lower_bound(MBTRI);
  AdjPrism_by_type::iterator hi_face_prism_miit = basicEntAdjacencies.get<EntType_mi_tag>().upper_bound(MBTRI);
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
  Range tets = _tets.subset_by_type(MBTET);
  Range::iterator tit = tets.begin();
  for(;tit!=tets.end();tit++) {
    ref_MoFEMFiniteElement_by_ent::iterator miit2 = ref_MoFEMFiniteElement.find(*tit);
    if(miit2==ref_MoFEMFiniteElement.end()) SETERRQ(PETSC_COMM_SELF,1,"this tet is not in MoFEMFiniteElement");
    //connectivity
    const EntityHandle* conn; 
    int num_nodes; 
    moab.get_connectivity(*tit,conn,num_nodes,true); 
    assert(num_nodes==4);
    for(int nn = 0;nn<num_nodes;nn++) {
      bool success = refinedMoFemEntities.modify(refinedMoFemEntities.get<MoABEnt_mi_tag>().find(conn[nn]),RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"can not set refinment bit level to tet node");
    }
    //get edges
    BitRefEdges parent_edges_bit(0);
    EntityHandle edge_new_nodes[6];
    fill(&edge_new_nodes[0],&edge_new_nodes[6],no_handle); 
    int split_edges[6];  
    fill(&split_edges[0],&split_edges[6],-1); 
    //hash map of nodes (RefMoFEMEntity) by edges (EntityHandle)
    map<EntityHandle /*edge*/,const RefMoFEMEntity* /*node*/> map_ref_nodes_by_edges; 
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
    EntityHandle new_tets_conns[8*4];
    fill(&new_tets_conns[0],&new_tets_conns[8*4],no_handle);
    int sub_type = -1,nb_new_tets = 0;
    switch (parent_edges_bit.count()) {
      case 0: {
	  ref_ents_by_ent::iterator tit_miit;
	  tit_miit = ref_ents_ent.find(*tit);
	  if(tit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  bool success = refinedMoFemEntities.modify(tit_miit,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"imposible tet");
	  Range tit_conn;
	  rval = moab.get_connectivity(&*tit,1,tit_conn,true); CHKERR_PETSC(rval);
	  for(Range::iterator nit = tit_conn.begin();nit!=tit_conn.end();nit++) {
	    ref_ents_by_ent::iterator nit_miit = ref_ents_ent.find(*nit);
	    if(nit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find face in refinedMoFemEntities");
	    bool success = refinedMoFemEntities.modify(nit_miit,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"imposible node");
	  }
	  Range tit_edges;
	  rval = moab.get_adjacencies(&*tit,1,1,false,tit_edges); CHKERR_PETSC(rval);
	  for(Range::iterator eit = tit_edges.begin();eit!=tit_edges.end();eit++) {
	    ref_ents_by_ent::iterator eit_miit = ref_ents_ent.find(*eit);
	    if(eit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find face in refinedMoFemEntities");
	    bool success = refinedMoFemEntities.modify(eit_miit,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"imposible edge");
	  }
	  Range tit_faces;
	  rval = moab.get_adjacencies(&*tit,1,2,false,tit_faces); CHKERR_PETSC(rval);
	  if(tit_faces.size()!=4) SETERRQ(PETSC_COMM_SELF,1,"existing tet in mofem databsee should have 4 adjacent edges");
	  for(Range::iterator fit = tit_faces.begin();fit!=tit_faces.end();fit++) {
	    ref_ents_by_ent::iterator fit_miit = ref_ents_ent.find(*fit);
	    if(fit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find face in refinedMoFemEntities");
	    bool success = refinedMoFemEntities.modify(fit_miit,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"imposible face");
	  }
	  continue;
	}
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
      //add this tet if exist to this ref level
      EntityHandle tet = miit_composite2->get_ref_ent();
      refinedMoFemEntities.modify(refinedMoFemEntities.find(tet),RefMoFEMEntity_change_add_bit(bit));
      //set bit that this element is in databse - no need to create it
      ref_tets_bit.set(tt,1);
      if(verbose>2) {
	ostringstream ss;
	ss << miit_composite2->get_RefMoFEMElement() << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    }
    if(miit_composite!=hi_miit_composite) {
      //if that tet has the same pattern of splitted edges it has to have the same number of refined 
      //children elements - if not thorw an error
      if(ref_tets_bit.count()!=(unsigned int)nb_new_tets) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    } else {
      //if this element was not refined or was reffined with diffrent patterns of splitted edges create new elements
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
	  int ref_type[2];
	  ref_type[0] = parent_edges_bit.count();
	  ref_type[1] = sub_type; 
	  rval = moab.tag_set_data(th_RefType,&ref_tets[tt],1,ref_type); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefParentHandle,&ref_tets[tt],1,&*tit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitLevel,&ref_tets[tt],1,&bit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitEdge,&ref_tets[tt],1,&parent_edges_bit); CHKERR_PETSC(rval);
	  //add refined entity
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_MoFEMEntity = refinedMoFemEntities.insert(RefMoFEMEntity(moab,ref_tets[tt]));
	  //add refined element
	  pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement 
	    = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_MoFEMEntity.first)));
	  //set bit that this element is now in databse
	  ref_tets_bit.set(tt);
	  if(verbose>2) {
	    ostringstream ss;
	    ss << "add tet: " << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement()) << endl;
	    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }
	}
      }
      //find parents for new edges and faces
      //get tet edges and faces
      Range tit_edges,tit_faces;
      rval = moab.get_adjacencies(&*tit,1,1,false,tit_edges); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(&*tit,1,2,false,tit_faces); CHKERR_PETSC(rval);
      Range edges_nodes[6],faces_nodes[4];
      //for edges - add ref nodes
      //edges_nodes[ee] - contains all nodes on edge ee inluding mid nodes if exist
      Range::iterator eit = tit_edges.begin();
      for(int ee = 0;eit!=tit_edges.end();eit++,ee++) {
	rval = moab.get_connectivity(&*eit,1,edges_nodes[ee],true); CHKERR_PETSC(rval);
	map<EntityHandle,const RefMoFEMEntity*>::iterator map_miit = map_ref_nodes_by_edges.find(*eit);
	if(map_miit!=map_ref_nodes_by_edges.end()) {
	  edges_nodes[ee].insert(map_miit->second->get_ref_ent());
	}
      }
      //for faces - add ref nodes
      //faces_nodes[ff] - contains all nodes on face ff inluding mid nodes if exist
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
      //add ref nodes to tet
      //tet_nodes contains all nodes on tet inluding mid edge nodes
      Range tet_nodes;
      rval = moab.get_connectivity(&*tit,1,tet_nodes,true); CHKERR_PETSC(rval);
      for(map<EntityHandle,const RefMoFEMEntity*>::iterator map_miit = map_ref_nodes_by_edges.begin();
	map_miit != map_ref_nodes_by_edges.end();map_miit++) {
	tet_nodes.insert(map_miit->second->get_ref_ent());
      }
      Range ref_edges;
      //get all all edges of refined tets
      rval = moab.get_adjacencies(ref_tets,nb_new_tets,1,true,ref_edges,Interface::UNION); CHKERR_PETSC(rval);
      //check for all ref edge and set parents
      for(Range::iterator reit = ref_edges.begin();reit!=ref_edges.end();reit++) {
	Range ref_edges_nodes;
	rval = moab.get_connectivity(&*reit,1,ref_edges_nodes,true); CHKERR_PETSC(rval);
	//check if ref edge is an coarse edge
	int ee = 0;
	for(;ee<6;ee++) {
	  //two nodes are common (node[0],node[1],ref_node (if exist))
	  //this tests if given edge is contained by edge of refined tetrahedral
	  if(intersect(edges_nodes[ee],ref_edges_nodes).size()==2) {
	    EntityHandle edge = tit_edges[ee];
	    rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&edge); CHKERR_PETSC(rval);
	    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*reit));
	    bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"imposible to set edge pranet");
	    if(p_ent.second) {
	      if(verbose>2) {
		ostringstream ss;
		ss << "edge parent: " << *(p_ent.first) << endl;
		PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }}
	    break;
	  }  
	}
	if(ee<6) continue; //this refined edge is contined by edge of tetrahedral
	//check if ref edge is in coarse face
	int ff = 0;
	for(;ff<4;ff++) {
	  //two nodes are common (node[0],node[1],ref_node (if exist))
	  //thi tests if givem edge is contained by face of  tetrahedral
	  if(intersect(faces_nodes[ff],ref_edges_nodes).size()==2) {
	    EntityHandle face = tit_faces[ff];
	    rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&face); CHKERR_PETSC(rval);
	    //add edge to refinedMoFemEntities
	    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*reit));
	    bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"imposible to set edge parent");
	    if(p_ent.second) {
	      if(verbose>2) {
		ostringstream ss;
		ss << "face parent: " << *(p_ent.first) << endl;
		PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }}
	    break;
	  }
	}
	if(ff<4) continue; //this refibed egde is contained by face of tetrahedral
	// check if ref edge is in coarse tetrahedral (i.e. that is internal edge of refined tetrahedral)
	if(intersect(tet_nodes,ref_edges_nodes).size()==2) {
	  rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&*tit); CHKERR_PETSC(rval);
	  //add edge to refinedMoFemEntities
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*reit));
	  bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"imposible to set edge parent");
	  if(p_ent.second) {
	    if(verbose>2) {
	      ostringstream ss;
	      ss << "tet parent: " << *(p_ent.first) << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }}
	  continue;
	}
	//refined edge is not child of any edge, face or tetrahedral, this is imposible edge
	SETERRQ(PETSC_COMM_SELF,1,"imposible refined edge");
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
	  //check if refined edge is contained by face of tetrahedral
	  if(intersect(faces_nodes[ff],ref_faces_nodes).size()==3) {
	    EntityHandle face = tit_faces[ff];
	    rval = moab.tag_set_data(th_RefParentHandle,&*rfit,1,&face); CHKERR_PETSC(rval);
	    //add face to refinedMoFemEntities
	    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*rfit));
	    bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"imposible to set face parent");
	    if(p_ent.second) {
	      if(verbose>2) {
		ostringstream ss;
		ss << "face parent: " << *(p_ent.first) << endl;
		PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }}
	    break;
	  }
	}
	if(ff<4) continue; //this face is contained by one of tetrahedrals 
	//check if ref face is in coarse tetrahedral
	//this is ref face which is contained by tetrahedral volume
	if(intersect(tet_nodes,ref_faces_nodes).size()==3) {
	  rval = moab.tag_set_data(th_RefParentHandle,&*rfit,1,&*tit); CHKERR_PETSC(rval);
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*rfit));
	  //add face to refinedMoFemEntities
	  bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"imposible to set face parent");
	  if(p_ent.second) {
	    if(verbose>2) {
	      ostringstream ss;
	      ss << "tet parent: " << *(p_ent.first) << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }}
	  continue;
	}
	SETERRQ(PETSC_COMM_SELF,1,"imposible refined face");
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  //FIXME: refinment is based on entity handlers, should work on global ids of nodes, this will allow parallelize agortihm in the future
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ENTs_by_ent;
  typedef RefMoFEMElement_multiIndex::index<Composite_of_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type ref_fe_by_composite;
  ref_fe_by_composite &ref_fe_by_comp = refinedMoFemElements.get<Composite_of_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
  //find all verices which parent is edge
  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityType_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents_by_comp = refinedMoFemEntities.get<Composite_EntityType_And_ParentEntityType_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents_by_comp.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents_by_comp.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    ref_parent_ents_view.insert(&*miit);
  }
  Range prisms;
  rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,false); CHKERR_PETSC(rval);
  Range::iterator pit = prisms.begin();
  for(;pit!=prisms.end();pit++) {
    ref_ENTs_by_ent::iterator miit_prism = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*pit);   
    if(miit_prism==refinedMoFemEntities.end()) SETERRQ(PETSC_COMM_SELF,1,"this prism is not in ref database");
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
      refinedMoFemEntities.modify(miit_prism,RefMoFEMEntity_change_add_bit(bit));
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
      refinedMoFemEntities.modify(refinedMoFemEntities.find(miit_composite2->get_ref_ent()),RefMoFEMEntity_change_add_bit(bit));
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
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,ref_prisms[pp]));
	  pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
	  try {
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
	  } catch (const char* msg) {
	    SETERRQ(PETSC_COMM_SELF,1,msg);
	  }
	  ref_prism_bit.set(pp);
	  ierr = add_prism_to_basicEntAdjacencies(ref_prisms[pp]); CHKERRQ(ierr);
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
PetscErrorCode FieldCore::refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ENTs_by_ent;
  ref_ENTs_by_ent::iterator miit = refinedMoFemEntities.find(meshset);
  if(miit==refinedMoFemEntities.end()) SETERRQ(PETSC_COMM_SELF,1,"this meshset is not in ref database");
  ierr = refine_get_childern(meshset,bit,meshset,MBEDGE,recursive,verb); CHKERRQ(ierr);
  ierr = refine_get_childern(meshset,bit,meshset,MBTRI,recursive,verb); CHKERRQ(ierr);
  ierr = refine_get_childern(meshset,bit,meshset,MBTET,recursive,verb); CHKERRQ(ierr);
  refinedMoFemEntities.modify(miit,RefMoFEMEntity_change_add_bit(bit));
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_get_ents(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb) {
  PetscFunctionBegin;
  Range ents;
  ierr = refine_get_ents(bit,mask,type,ents,verb); CHKERRQ(ierr);
  rval = moab.add_entities(meshset,ents); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_get_ents(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  RefMoFEMEntity_multiIndex::index<EntType_mi_tag>::type::iterator miit = refinedMoFemEntities.get<EntType_mi_tag>().lower_bound(type);
  for(;miit!=refinedMoFemEntities.get<EntType_mi_tag>().upper_bound(type);miit++) {
    BitRefLevel bit2 = miit->get_BitRefLevel(); 
    if((bit2&mask) != bit2) continue;
    if(verb > 2) {
      ostringstream ss;
      ss << *miit << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    if(verb > 3) {
      ostringstream ss;
      ss << bit << endl;
      ss << mask << endl;
      ss << bit2 << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    if((bit2&bit).any()) {
      if(verb > 3) {
	ostringstream ss;
	ss << "add ent to meshset" << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
      EntityHandle ent = miit->get_ref_ent();
      ents.insert(ent);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_get_ents(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) {
  PetscFunctionBegin;
  RefMoFEMEntity_multiIndex::iterator miit = refinedMoFemEntities.begin();
  for(;miit!=refinedMoFemEntities.end();miit++) {
    BitRefLevel bit2 = miit->get_BitRefLevel(); 
    if((bit2&mask) != bit2) continue;
    if((bit2&bit).any()) {
      EntityHandle ent = miit->get_ref_ent();
      rval = moab.add_entities(meshset,&ent,1); CHKERR_PETSC(rval);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_get_ents(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents) {
  PetscFunctionBegin;
  RefMoFEMEntity_multiIndex::iterator miit = refinedMoFemEntities.begin();
  for(;miit!=refinedMoFemEntities.end();miit++) {
    BitRefLevel bit2 = miit->get_BitRefLevel(); 
    if((bit2&mask) != bit2) continue;
    if((bit2&bit).any()) {
      EntityHandle ent = miit->get_ref_ent();
      ents.insert(ent);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_get_childern(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityHandle_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedMoFemEntities.get<Composite_EntityHandle_And_ParentEntityType_mi_tag>();
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
PetscErrorCode FieldCore::problem_get_FE(const string &problem_name,const string &fe_name,const EntityHandle meshset) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no sach poblem like < %s >",problem_name.c_str());
  NumeredMoFEMFiniteElement_multiIndex &numeredFiniteElements = const_cast<NumeredMoFEMFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
  NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator miit = numeredFiniteElements.get<MoFEMFiniteElement_name_mi_tag>().lower_bound(fe_name);
  for(;miit!=numeredFiniteElements.get<MoFEMFiniteElement_name_mi_tag>().upper_bound(fe_name);miit++) {
    EntityHandle ent = miit->get_ent();
    rval = moab.add_entities(meshset,&ent,1); CHKERR_PETSC(rval);
    int part = miit->get_part();
    rval = moab.tag_set_data(th_Part,&ent,1,&part); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  moabCubitMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    ierr = miit->get_Cubit_msId_entities_by_dimension(moab,dimension,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  moabCubitMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    ierr = miit->get_Cubit_msId_entities_by_dimension(moab,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_Cubit_msId_entities_by_dimension(msId,Cubit_BC_bitset(CubitBCType),dimension,entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType,
  Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_Cubit_msId_entities_by_dimension(msId,Cubit_BC_bitset(CubitBCType),entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset) {
  PetscFunctionBegin;
  moabCubitMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    meshset = miit->meshset;
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_CubitBCType_meshsets(const unsigned int CubitBCType,Range &meshsets) {
  PetscFunctionBegin;
  moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<CubitMeshSets_mi_tag>().lower_bound(CubitBCType);
  moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator 
    hi_miit = cubit_meshsets.get<CubitMeshSets_mi_tag>().upper_bound(CubitBCType);
  for(;miit!=hi_miit;miit++) {
    meshsets.insert(miit->meshset);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::VecCreateGhost(const string &name,RowColData rc,Vec *V) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_local_idx;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spellig)",name.c_str());
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
PetscErrorCode FieldCore::VecScatterCreate(Vec xin,string &x_problem,RowColData x_rc,Vec yin,string &y_problem,RowColData y_rc,VecScatter *newctx,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_x = moFEMProblems_set.find(x_problem);
  if(p_x==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spellig)",x_problem.c_str());
  moFEMProblems_by_name::iterator p_y = moFEMProblems_set.find(y_problem);
  if(p_y==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spellig)",y_problem.c_str());
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_glob_idx;
  dofs_by_glob_idx::iterator y_dit,hi_y_dit;
  switch (y_rc) {
    case Row:
      y_dit = p_y->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(p_y->get_nb_local_dofs_row()); // should be lower
      break;
    case Col:
      y_dit = p_y->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().lower_bound(p_y->get_nb_local_dofs_col()); // should be lower
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type dofs_by_uid;
  const dofs_by_uid* x_numered_dofs_by_uid;
  switch (x_rc) {
    case Row:
      x_numered_dofs_by_uid = &(p_x->numered_dofs_rows.get<Unique_mi_tag>());
      break;
    case Col:
      x_numered_dofs_by_uid = &(p_x->numered_dofs_cols.get<Unique_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  vector<int> idx(0),idy(0);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  for(;y_dit!=hi_y_dit;y_dit++) {
    if(y_dit->get_part()!=pcomm->rank()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    dofs_by_uid::iterator x_dit;
    x_dit = x_numered_dofs_by_uid->find(y_dit->get_unique_id());
    if(x_dit==x_numered_dofs_by_uid->end()) continue;
    idx.push_back(x_dit->get_petsc_gloabl_dof_idx());
    idy.push_back(y_dit->get_petsc_gloabl_dof_idx());
  }
  IS ix,iy;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD,idx.size(),&idx[0],PETSC_USE_POINTER,&ix); CHKERRQ(ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD,idy.size(),&idy[0],PETSC_USE_POINTER,&iy); CHKERRQ(ierr);
  if(verb>3) {
    ISView(ix,PETSC_VIEWER_STDOUT_WORLD);
    ISView(iy,PETSC_VIEWER_STDOUT_WORLD);
  }
  ierr = ::VecScatterCreate(xin,ix,yin,iy,newctx); CHKERRQ(ierr);
  ierr = ISDestroy(&ix); CHKERRQ(ierr);
  ierr = ISDestroy(&iy); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::MatCreateMPIAIJWithArrays(const string &name,Mat *Aij,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = partition_create_Mat<Part_mi_tag>(name,Aij,MATMPIAIJ,false,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::MatCreateSeqAIJWithArrays(const string &name,Mat *Aij,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = partition_create_Mat<Part_mi_tag>(name,Aij,MATAIJ,false,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::set_local_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_local_idx;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str()); 
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
PetscErrorCode FieldCore::set_global_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type dofs_by_global_idx;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str());
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
PetscErrorCode FieldCore::set_other_global_VecCreateGhost(
  const string &name,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,
  int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type dofs_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found",name.c_str());
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
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit = moabFields.get<FieldName_mi_tag>().find(cpy_field_name);
  if(cpy_fit==moabFields.get<FieldName_mi_tag>().end()) {
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
	    DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx>::type::iterator diiiit;
	    diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	    if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx>().end()) {
	      EntityHandle ent = miit->get_ent();
	      rval = moab.add_entities(cpy_fit->get_meshset(),&ent,1); CHKERR_PETSC(rval);
	      //create field moabent
	      ApproximationOrder order = miit->get_max_order();
	      pair<MoFEMEntity_multiIndex::iterator,bool> p_e_miit;
	      try {
		MoFEMEntity moabent(moab,cpy_fit->get_MoFEMField_ptr(),miit->get_RefMoFEMEntity_ptr());
		p_e_miit = entsMoabField.insert(moabent);
	      } catch (const std::exception& ex) {
		ostringstream ss;
		ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
		SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	      }
	      if(p_e_miit.first->get_max_order()<order) {
		bool success = entsMoabField.modify(p_e_miit.first,MoFEMEntity_change_order(moab,order));
		if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
	      }
	      //create field moabdof
	      DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent>::type::iterator hi_diit,diit;
	      diit = dofsMoabField.get<Composite_Name_And_Ent>().lower_bound(boost::make_tuple(field_name,miit->get_ent()));
	      hi_diit = dofsMoabField.get<Composite_Name_And_Ent>().upper_bound(boost::make_tuple(field_name,miit->get_ent()));
	      for(;diit!=hi_diit;diit++) {
		DofMoFEMEntity mdof(&*(p_e_miit.first),diit->get_dof_order(),diit->get_dof_rank(),diit->get_EntDofIdx());
		pair<DofMoFEMEntity_multiIndex::iterator,bool> cpy_p_diit;
		cpy_p_diit = dofsMoabField.insert(mdof);
		if(cpy_p_diit.second) {
		  bool success = dofsMoabField.modify(cpy_p_diit.first,DofMoFEMEntity_active_change(true));
		  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
		}
	      }
	      diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	      if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
	    }
	    diiiit->get_FieldData() = array[miit->get_petsc_gloabl_dof_idx()];
	    if(verb > 1) {
	      ostringstream ss;
	      ss << *diiiit << "set " << array[miit->get_petsc_gloabl_dof_idx()] << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }
	  }
	    //if(verb > 0) {
	      //cerr << "AAAAAAAAAAAA\n";
	      //ierr = check_number_of_ents_in_ents_field(cpy_field_name); CHKERRQ(ierr);
	    //}
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      ierr = VecRestoreArray(V_glob,&array); CHKERRQ(ierr);
      ierr = VecDestroy(&V_glob); CHKERRQ(ierr);
      ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    }
    break;
    case SCATTER_FORWARD: {
	for(;miit!=hi_miit;miit++) {
	  if(pcomm->rank()!=miit->get_part()) continue;
	  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx>::type::iterator diiiit;
	  diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	  if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx>().end()) {
	    SETERRQ(PETSC_COMM_SELF,1,"no data to fill the vector (top tip: you want scatter forward of scatter reverse?)");
	  }
	  ierr = VecSetValue(V,miit->get_petsc_gloabl_dof_idx(),diiiit->get_FieldData(),mode); CHKERRQ(ierr);
	}
	ierr = VecAssemblyBegin(V); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(V); CHKERRQ(ierr);
      } 
      break;  
    default:
     SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_3dENTS_sides(const int msId,const Cubit_BC_bitset CubitBCType,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  moabCubitMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    ierr = FieldCore::get_msId_3dENTS_sides(miit->meshset,recursive,verb); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_3dENTS_sides(const EntityHandle SideSet,const bool recursive,int verb) {
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
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"adj. node if ents3d %u\n",nodes.size());
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
  skin_edges = subtract(skin_edges,skin_faces_edges); // from skin edges subtract edges from skin faces of 3d ents (only internal edges)
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"subtract skin_edges %u\n",skin_edges.size());
  Range skin_nodes;
  rval = moab.get_connectivity(skin_edges,skin_nodes,true); CHKERR_PETSC(rval);
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"subtract skin_nodes %u\n",skin_nodes.size());
  nodes = subtract(nodes,skin_nodes); // nodes adjacent to all splitted face edges except those on internal edge
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"adj. node if ents3d but not on the internal edge %u\n",nodes.size());
  // ents3 that are adjacent to nodes on splitted faces but not those which are on the nodes on internal edgea
  ents3d.clear();
  rval = moab.get_adjacencies(nodes,3,true,ents3d,Interface::UNION); CHKERR_PETSC(rval);
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"adj. ents3d to nodes %u\n",ents3d.size());
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
    PetscPrintf(PETSC_COMM_WORLD,"Nb. of other side ents3d in set %u\n",other_side.size());
  }
  if(verb>3) {
    ierr = moab.write_file("side.vtk","VTK","",&children[0],1); CHKERRQ(ierr);
    ierr = moab.write_file("other_side.vtk","VTK","",&children[1],1); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_3dENTS_split_sides(
  const EntityHandle meshset,const BitRefLevel &bit,
  const int msId,const Cubit_BC_bitset CubitBCType,const bool add_iterfece_entities,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  moabCubitMeshSet_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_mi_tag>().end()) {
    ierr = FieldCore::get_msId_3dENTS_split_sides(
      meshset,bit,miit->meshset,add_iterfece_entities,recursive,verb); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_3dENTS_split_sides(
  const EntityHandle meshset,const BitRefLevel &bit,
  const EntityHandle SideSet,const bool add_iterfece_entities,const bool recursive,
  int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  vector<EntityHandle> children;
  //get children meshsets
  rval = moab.get_child_meshsets(SideSet,children);  CHKERR_PETSC(rval);
  if(children.size()!=2) {
    SETERRQ(PETSC_COMM_SELF,1,"should be 2 child meshsets, each of them contains tets on two sides of interface");
  }
  //get child of child of the first meshset 
  //contains vertices on interface, excluding those on crack front
  vector<EntityHandle> children_nodes_and_skin_edges;
  rval = moab.get_child_meshsets(children[0],children_nodes_and_skin_edges);  CHKERR_PETSC(rval);
  if(children_nodes_and_skin_edges.size()!=1) {
    SETERRQ(PETSC_COMM_SELF,1,"should be 1 child of the child, containing vertices on the interface");
  }
  //faces of interface
  Range triangles;
  rval = moab.get_entities_by_type(SideSet,MBTRI,triangles,recursive);  CHKERR_PETSC(rval);
  //tetrahedrasl on "father" side
  Range side_ents3d;
  rval = moab.get_entities_by_type(children[0],MBTET,side_ents3d,false);  CHKERR_PETSC(rval);
  //tetrahedral on "mather" side
  Range other_ents3d;
  rval = moab.get_entities_by_type(children[1],MBTET,other_ents3d,false);  CHKERR_PETSC(rval);
  //nodes on interface but not on crack front (those should not be splitted)
  Range nodes;
  rval = moab.get_entities_by_type(children_nodes_and_skin_edges[0],MBVERTEX,nodes,false);  CHKERR_PETSC(rval);
  Range meshset_3d_ents,meshset_2d_ents;
  rval = moab.get_entities_by_dimension(meshset,3,meshset_3d_ents,true); CHKERR_PETSC(rval);
  Range meshset_tets = meshset_3d_ents.subset_by_type(MBTET);
  rval = moab.get_adjacencies(meshset_tets,2,false,meshset_2d_ents,moab::Interface::UNION); CHKERR_PETSC(rval);
  side_ents3d = intersect(meshset_3d_ents,side_ents3d);
  other_ents3d = intersect(meshset_3d_ents,other_ents3d); 
  triangles = intersect(meshset_2d_ents,triangles);
  if(verb>3) {
    PetscPrintf(PETSC_COMM_WORLD,"triangles %u\n",triangles.size());
    PetscPrintf(PETSC_COMM_WORLD,"side_ents3d %u\n",side_ents3d.size());
    PetscPrintf(PETSC_COMM_WORLD,"nodes %u\n",nodes.size());
  }
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ents_by_ent_type;
  ref_ents_by_ent_type &ref_ents_by_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>();
  //maps nodes on "father" and "mather" side
  map<
    EntityHandle /*node on "mather" side*/,
    EntityHandle /*node on "father" side*/
    > map_nodes;
  //add new nodes on interface and create map
  Range::iterator nit = nodes.begin();
  double coord[3];
  for(;nit!=nodes.end();nit++) {
    rval = moab.get_coords(&*nit,1,coord); CHKERR_PETSC(rval);	
    EntityHandle new_node;
    rval = moab.create_vertex(coord,new_node); CHKERR(rval);
    map_nodes[*nit] = new_node;
    ref_ents_by_ent_type::iterator miit_ref_ent = ref_ents_by_ent.find(*nit);
    if(miit_ref_ent == ref_ents_by_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find node in MoFEM database");
    //create new node on "father" side
    //parent is node on "mather" side
    rval = moab.tag_set_data(th_RefParentHandle,&new_node,1,&*nit); CHKERR_PETSC(rval);
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,new_node));
    //set ref bit level to node on "father" side
    refinedMoFemEntities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(bit));
    //set ref bit level to node on "mather" side
    refinedMoFemEntities.modify(miit_ref_ent,RefMoFEMEntity_change_add_bit(bit));
  }
  //crete meshset for new mesh bit level
  EntityHandle meshset_for_bit_level;
  rval = moab.create_meshset(MESHSET_SET,meshset_for_bit_level); CHKERR_PETSC(rval);
  //subtract those elements which will be refined, i.e. disconetcted form other side elements, and connected to new prisms, if they area created
  meshset_3d_ents = subtract(meshset_3d_ents,side_ents3d);
  rval = moab.add_entities(meshset_for_bit_level,meshset_3d_ents); CHKERR_PETSC(rval);
  for(int dd = 0;dd<3;dd++) {
    Range ents_dd;
    rval = moab.get_adjacencies(meshset_3d_ents,dd,false,ents_dd,moab::Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(meshset_for_bit_level,ents_dd); CHKERR_PETSC(rval);
  }
  //create new tets on "father" side
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
	nb_new_conn++;
	if(verb>6) {
	  PetscPrintf(PETSC_COMM_WORLD,"nodes %u -> %d\n",conn[ii],new_conn[ii]);
	}
      } else {
	new_conn[ii] = conn[ii];
      }
    }
    if(nb_new_conn==0) {
      if(verb>3) {
	EntityHandle meshset_error_out;
	rval = moab.create_meshset(MESHSET_SET,meshset_error_out); CHKERR_PETSC(rval);
	rval = moab.add_entities(meshset_error_out,&*tit,1); CHKERR_PETSC(rval);
	ierr = moab.write_file("error_out.vtk","VTK","",&meshset_error_out,1); CHKERRQ(ierr);
      }
      SETERRQ(PETSC_COMM_SELF,1,"database insonistency, in side_ent3 is a tet which has no common node with interface");
    }
    //here is created new tet
    EntityHandle tet;
    rval = moab.create_element(MBTET,new_conn,4,tet); CHKERR_PETSC(rval);
    rval = moab.tag_set_data(th_RefParentHandle,&tet,1,&*tit); CHKERR_PETSC(rval);
    rval = moab.add_entities(meshset_for_bit_level,&tet,1); CHKERR_PETSC(rval);
    rval = moab.add_entities(meshset_for_bit_level,new_conn,4); CHKERR_PETSC(rval);
    new_tets.insert(tet);
  }
  Range new_ents; 
  //create new entities by adjecies form new tets
  rval = moab.get_adjacencies(new_tets,1,true,new_ents,Interface::UNION); CHKERR_PETSC(rval);
  rval = moab.get_adjacencies(new_tets,2,true,new_ents,Interface::UNION); CHKERR_PETSC(rval);
  Range ents; 
  //add new edges and triangles to mofem database
  rval = moab.get_adjacencies(triangles,1,false,ents,Interface::UNION); CHKERR_PETSC(rval);
  Range new_ents_in_database; //this range contains all new entities
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
    if(miit_ref_ent == ref_ents_by_ent.end()) {
      SETERRQ(PETSC_COMM_SELF,1,"this entity (edge or tri) should be already in database");
    }
    Range new_ent; //contains all entities (edges or triangles) added to mofem database
    switch (moab.type_from_handle(*eit)) {
      case MBTRI: {
	  //get entity based on its connectivity
	  rval = moab.get_adjacencies(new_conn,3,2,false,new_ent); CHKERR_PETSC(rval);
	  if(new_ent.size() != 1) SETERRQ(PETSC_COMM_SELF,1,"this tri should be in moab database"); 
	  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"new_ent %u\n",new_ent.size());
	  //add prism element
	  if(add_iterfece_entities) {
	    //set prism connectivity
	    EntityHandle prism_conn[6] = { 
	      conn[0],conn[1],conn[2],
	      new_conn[0],new_conn[1],new_conn[2] 
	    };
	    //cerr << 
	    //  conn[0] << " " << conn[1] << " " << conn[2] << " ::: " 
	    //  << new_conn[0] << " " << new_conn[1] << " " << new_conn[2] << endl;
	    EntityHandle prism = no_handle;
	    rval = moab.create_element(MBPRISM,prism_conn,6,prism); CHKERR_PETSC(rval);
	    ierr = add_prism_to_basicEntAdjacencies(prism,verb/*nb_new_conn < 3 ? 1 : 0*/); CHKERRQ(ierr);
	    rval = moab.add_entities(meshset_for_bit_level,&prism,1); CHKERR_PETSC(rval);
	  }
	} break;
      case MBEDGE: {
	  rval = moab.get_adjacencies(new_conn,2,1,false,new_ent); CHKERR_PETSC(rval);
	  if(new_ent.size()!=1) {
	    SETERRQ(PETSC_COMM_SELF,1,"this edge should be in moab database");
	  }
	} break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"huston we have problem !!!");
    }
    if(new_ent.size()!=1) {
      SETERRQ1(PETSC_COMM_SELF,1,"new_ent.size() = %u, size always should be 1",new_ent.size());
    }
    //set parent 
    rval = moab.tag_set_data(th_RefParentHandle,&*new_ent.begin(),1,&*eit); CHKERR_PETSC(rval);
    //add to database
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,new_ent[0]));
    refinedMoFemEntities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(bit));
    new_ents_in_database.insert(new_ent.begin(),new_ent.end());
  }
  //all other entities, some ents like triangles and faces on the side of tets
  Range side_adj_faces_and_edges;
  rval = moab.get_adjacencies(side_ents3d,1,true,side_adj_faces_and_edges,Interface::UNION); CHKERR_PETSC(rval);
  rval = moab.get_adjacencies(side_ents3d,2,true,side_adj_faces_and_edges,Interface::UNION); CHKERR_PETSC(rval);
  //subtract entities already added to mofem database
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
    if(miit_ref_ent == ref_ents_by_ent.end()) {
      SETERRQ(PETSC_COMM_SELF,1,"entity should be in MoFem database");
    }
    Range new_ent;
    switch (moab.type_from_handle(*eit)) {
      case MBTRI: {
	  rval = moab.get_adjacencies(new_conn,3,2,false,new_ent); CHKERR_PETSC(rval);
	}
	break;
      case MBEDGE: {
	  rval = moab.get_adjacencies(new_conn,2,1,false,new_ent); CHKERR_PETSC(rval);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"huston we have problem");   
    }
    if(new_ent.size()!=1) {
      SETERRQ1(PETSC_COMM_SELF,1,"database insonistency, new_ent.size() = %u",new_ent.size());
    }
    //add entity to mofem database
    rval = moab.tag_set_data(th_RefParentHandle,&*new_ent.begin(),1,&*eit); CHKERR_PETSC(rval);
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent 
      = refinedMoFemEntities.insert(RefMoFEMEntity(moab,new_ent[0]));
    refinedMoFemEntities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(bit));
    if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"new_ent %u\n",new_ent.size());
    new_ents_in_database.insert(new_ent.begin(),new_ent.end());
  }
  //finalise by adding new tets and prism ti bitlelvel
  ierr = seed_ref_level_3D(meshset_for_bit_level,bit); CHKERRQ(ierr);
  rval = moab.delete_entities(&meshset_for_bit_level,1); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_prism_to_basicEntAdjacencies(const EntityHandle prism,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  vector<EntityHandle> Ents(8,no_handle);
  try {
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,prism));
    pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
    if(p_ent.second) {
      p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
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
      //set bit common for faces with side number 3 and 4
      RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*face_side3.begin());
      if(miit_ref_ent!=refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) {
	BitRefLevel bit = miit_ref_ent->get_BitRefLevel();
	if(face_side4.empty()) SETERRQ(PETSC_COMM_SELF,1,"database insonistency");
	miit_ref_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*face_side4.begin());
	if(miit_ref_ent!=refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) {
	  bit &= miit_ref_ent->get_BitRefLevel();
	  refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	}
      }
    } 
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::problem_basic_method_preProcess(const string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  // find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in databse %s",problem_name.c_str());
  // finite element
  typedef NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type FEs_by_composite;
  ierr = method.set_problem(&*p_miit); CHKERRQ(ierr);
  ierr = method.set_fields(&moabFields); CHKERRQ(ierr);
  ierr = method.set_ents_multiIndex(&entsMoabField); CHKERRQ(ierr);
  ierr = method.set_dofs_multiIndex(&dofsMoabField); CHKERRQ(ierr);
  ierr = method.set_fes_multiIndex(&finiteElements); CHKERRQ(ierr);
  ierr = method.set_fes_data_multiIndex(&finiteElementsMoFEMEnts); CHKERRQ(ierr);
  ierr = method.set_adjacencies(&entFEAdjacencies); CHKERRQ(ierr);
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::problem_basic_method_postProcess(const string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  // find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in databse %s",problem_name.c_str());
  // finite element
  typedef NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type FEs_by_composite;
  ierr = method.set_problem(&*p_miit); CHKERRQ(ierr);
  ierr = method.set_fields(&moabFields); CHKERRQ(ierr);
  ierr = method.set_ents_multiIndex(&entsMoabField); CHKERRQ(ierr);
  ierr = method.set_dofs_multiIndex(&dofsMoabField); CHKERRQ(ierr);
  ierr = method.set_fes_multiIndex(&finiteElements); CHKERRQ(ierr);
  ierr = method.set_fes_data_multiIndex(&finiteElementsMoFEMEnts); CHKERRQ(ierr);
  ierr = method.set_adjacencies(&entFEAdjacencies); CHKERRQ(ierr);
  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::loop_finite_elements(const string &problem_name,const string &fe_name,FEMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  //
  ierr = loop_finite_elements(problem_name,fe_name,method,pcomm->rank(),pcomm->rank(),verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::loop_finite_elements(
  const string &problem_name,const string &fe_name,FEMethod &method,int lower_rank,int upper_rank,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  // find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in databse %s",problem_name.c_str());
  // finite element
  typedef NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type FEs_by_composite;
  method.fe_name = fe_name;
  ierr = method.set_problem(&*p_miit); CHKERRQ(ierr);
  ierr = method.set_fields(&moabFields); CHKERRQ(ierr);
  ierr = method.set_ents_multiIndex(&entsMoabField); CHKERRQ(ierr);
  ierr = method.set_dofs_multiIndex(&dofsMoabField); CHKERRQ(ierr);
  ierr = method.set_fes_multiIndex(&finiteElements); CHKERRQ(ierr);
  ierr = method.set_fes_data_multiIndex(&finiteElementsMoFEMEnts); CHKERRQ(ierr);
  ierr = method.set_adjacencies(&entFEAdjacencies); CHKERRQ(ierr);
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);
  FEs_by_composite &numeredFiniteElements = 
    (const_cast<NumeredMoFEMFiniteElement_multiIndex&>(p_miit->numeredFiniteElements)).get<Composite_mi_tag>();
  FEs_by_composite::iterator miit = numeredFiniteElements.lower_bound(boost::make_tuple(fe_name,lower_rank));
  FEs_by_composite::iterator hi_miit = numeredFiniteElements.upper_bound(boost::make_tuple(fe_name,upper_rank));
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
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  }
  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type dofs_by_name;
  // find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem not in databse %s",problem_name.c_str());
  ierr = method.set_fields(&moabFields); CHKERRQ(ierr);
  ierr = method.set_ents_multiIndex(&entsMoabField); CHKERRQ(ierr);
  ierr = method.set_dofs_multiIndex(&dofsMoabField); CHKERRQ(ierr);
  ierr = method.set_fes_multiIndex(&finiteElements); CHKERRQ(ierr);
  ierr = method.set_fes_data_multiIndex(&finiteElementsMoFEMEnts); CHKERRQ(ierr);
  ierr = method.set_adjacencies(&entFEAdjacencies); CHKERRQ(ierr);
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
PetscErrorCode FieldCore::get_problem(const string &problem_name,const MoFEMProblem **problem_ptr) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found, (top tip: check spelling)",problem_name.c_str());
  *problem_ptr = &*p_miit;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_dofs(const DofMoFEMEntity_multiIndex **dofsMoabField_ptr) {
  PetscFunctionBegin;
  *dofsMoabField_ptr = &dofsMoabField;
  PetscFunctionReturn(0);
}
MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator FieldCore::get_ent_moabfield_by_name_begin(const string &field_name) {
  return entsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
}
MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator FieldCore::get_ent_moabfield_by_name_end(const string &field_name) {
  return entsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator FieldCore::get_dofs_by_name_begin(const string &field_name) {
  return dofsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
}
DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator FieldCore::get_dofs_by_name_end(const string &field_name) {
  return dofsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent>::type::iterator FieldCore::get_dofs_by_name_and_ent_begin(const string &field_name,const EntityHandle ent) {
  return dofsMoabField.get<Composite_Name_And_Ent>().lower_bound(boost::make_tuple(field_name,ent));
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent>::type::iterator FieldCore::get_dofs_by_name_and_ent_end(const string &field_name,const EntityHandle ent) {
  return dofsMoabField.get<Composite_Name_And_Ent>().upper_bound(boost::make_tuple(field_name,ent));
}
PetscErrorCode FieldCore::get_finite_elements(const MoFEMFiniteElement_multiIndex **finiteElements_ptr) {
  PetscFunctionBegin;
  *finiteElements_ptr = &finiteElements;
  PetscFunctionReturn(0);
}
EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator FieldCore::get_fes_moabfield_by_name_begin(const string &fe_name) {
  return finiteElementsMoFEMEnts.get<MoFEMFiniteElement_name_mi_tag>().lower_bound(fe_name);
}
EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator FieldCore::get_fes_moabfield_by_name_end(const string &fe_name) {
  return finiteElementsMoFEMEnts.get<MoFEMFiniteElement_name_mi_tag>().upper_bound(fe_name);
}
PetscErrorCode FieldCore::check_number_of_ents_in_ents_field(const string& name) {
  PetscFunctionBegin;
  MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator it = entsMoabField.get<FieldName_mi_tag>().find(name);
  EntityHandle meshset = it->get_meshset();
  int num_entities;
  rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
  if(num_entities != distance(entsMoabField.get<FieldName_mi_tag>().lower_bound(name),entsMoabField.get<FieldName_mi_tag>().upper_bound(name))) {
    SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and multiindex < %s >",name.c_str());
  }
  PetscFunctionReturn(0);
}

}
