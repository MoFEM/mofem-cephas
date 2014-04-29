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
#include<version.h>

namespace MoFEM {

const static int debug = 1;

PetscErrorCode print_MoFem_verison() {
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"version %d.%d.%d\n",MoFEM_VERSION_MAJOR,MoFEM_VERSION_MINOR,MoFEM_VERSION_BUILD);
  PetscPrintf(PETSC_COMM_WORLD,"git commit id %s\n",GIT_SHA1_NAME);
  PetscFunctionReturn(0);
}

FieldCore::FieldCore(Interface& _moab,int _verbose): 
  moab(_moab),verbose(_verbose) {

  const EntityHandle root_meshset = moab.get_root_set();
  if(verbose>0) {
    print_MoFem_verison();
  }

  // Version
  Tag th_version;
  stringstream strs_version;
  strs_version << "MoFEM_version_" << MoFEM_VERSION_MAJOR << "." << MoFEM_VERSION_MINOR << "." << MoFEM_VERSION_BUILD;
  string version = strs_version.str();
  rval = moab.tag_get_handle("_MoFEM_VERSION",version.size()*sizeof(char),MB_TYPE_OPAQUE,
    th_version,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,NULL); CHKERR_THROW(rval);
  const char *ptr_version = version.c_str();
  rval = moab.tag_set_data(th_version,&root_meshset,1,ptr_version); CHKERR_THROW(rval);

  //tags saved in vtk-files
  const int def_part = -1;
  rval = moab.tag_get_handle("PARTITION",1,MB_TYPE_INTEGER,th_Part,MB_TAG_CREAT|MB_TAG_SPARSE,&def_part); CHKERR_THROW(rval);

  //Tags Ref
  EntityHandle def_handle = 0;
  rval = moab.tag_get_handle("_RefParentHandle",1,MB_TYPE_HANDLE,th_RefParentHandle,MB_TAG_CREAT|MB_TAG_SPARSE,&def_handle); CHKERR_THROW(rval);
  const int def_type[] = {0,0};
  rval = moab.tag_get_handle("_RefType",2,MB_TYPE_INTEGER,th_RefType,MB_TAG_CREAT|MB_TAG_SPARSE,def_type); CHKERR_THROW(rval);
  BitRefLevel def_bit_level = 0;
  rval = moab.tag_get_handle("_RefBitLevel",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_RefBitLevel,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); CHKERR_THROW(rval);
  BitRefEdges def_bit_egde = 0;
  rval = moab.tag_get_handle("_RefBitEdge",sizeof(BitRefEdges),MB_TYPE_OPAQUE,
    th_RefBitEdge,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_egde); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_RefFEMeshset",1,MB_TYPE_HANDLE,
    th_RefFEMeshset,MB_TAG_CREAT|MB_TAG_SPARSE,&root_meshset);
    
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
  clear_map();
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

bool FieldCore::check_field(const string &name) const {
  typedef MoFEMField_multiIndex::index<FieldName_mi_tag>::type field_set_by_name;
  const field_set_by_name &set = moabFields.get<FieldName_mi_tag>();
  field_set_by_name::iterator miit = set.find(name);
  if(miit==set.end()) return false;
  return true;
}
const MoFEMField* FieldCore::get_field_structure(const string& name) {
  typedef MoFEMField_multiIndex::index<FieldName_mi_tag>::type field_set_by_name;
  const field_set_by_name &set = moabFields.get<FieldName_mi_tag>();
  field_set_by_name::iterator miit = set.find(name);
  if(miit==set.end()) return NULL;
  return &*miit;
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
  refinedMoFemEntities.clear();
  refinedMoFemElements.clear();
  moabFields.clear();
  entsMoabField.clear();
  dofsMoabField.clear();
  finiteElements.clear();
  finiteElementsMoFEMEnts.clear();
  entFEAdjacencies.clear();
  moFEMProblems.clear();
  cubit_meshsets.clear();
  PetscFunctionReturn(0);
} 
PetscErrorCode FieldCore::add_field(const string& name,const BitFieldId id,const FieldSpace space,const ApproximationRank rank,enum MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fit;
  fit = moabFields.get<FieldName_mi_tag>().find(name);
  if(fit != moabFields.get<FieldName_mi_tag>().end() ) {
    if(bh == MF_EXCL) {
      SETERRQ1(PETSC_COMM_SELF,1,"field is <%s> in database",name.c_str());
    }
  } else {
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
  	"field not inserted %s (top tip, it could be already there)",
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
PetscErrorCode FieldCore::rebuild_database(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_map(); CHKERRQ(ierr);
  ierr = initialiseDatabseInformationFromMesh(verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::initialiseDatabseInformationFromMesh(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
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
	if(verb > 0) {
	  ostringstream ss;
	  ss << "read cubit" << base_meshset << endl;
	  PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	}
	//PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	//PetscSynchronizedFlush(PETSC_COMM_WORLD); 
	//ierr = seed_ref_level_MESHSET(*mit,0); CHKERRQ(ierr);
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
	try {
	switch (moab.type_from_handle(*eit)) {
	  case MBVERTEX:
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_VERTEX(moab,&*p_ref_ent.first)));
	    break;
	  case MBEDGE:
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_EDGE(moab,&*p_ref_ent.first)));
	    break;
	  case MBTRI:
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TRI(moab,&*p_ref_ent.first)));
	    break;
	  case MBTET:
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_ref_ent.first)));
	    break;
	  case MBPRISM:
	    ierr = add_prism_to_mofem_database(*eit,verb); CHKERRQ(ierr);
	    break;
	  case MBENTITYSET:
  	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_MESHSET(moab,&*p_ref_ent.first)));
	    break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"Only finite elements of type MBTET, MBPRISM and MBENTITYSET are implemented");
	}
	if(p_MoFEMFiniteElement.second) {
	  //PetscPrintf(PETSC_COMM_WORLD,"Warrning: this entity should be already in refined finite elements database");
	  //SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, this entity should be already in refined finite elements database");
	}
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
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
  //build ref entities meshset
  for(int dd = 0;dd<=3;dd++) {
    Range ents;
    rval = moab.get_entities_by_dimension(0,dd,ents,false); CHKERR_PETSC(rval);
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();eit++) {
      switch (moab.type_from_handle(*eit)) {
	case MBVERTEX:
	case MBEDGE:
	case MBTRI:
	case MBTET:
	case MBPRISM:
	  break;
	default: 
	 continue;
      }
      RefMoFEMEntity mofem_ent(moab,*eit);
      BitRefLevel bit = mofem_ent.get_BitRefLevel();
      if(bit.none()) {
	continue;
      }
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p;
      p = refinedMoFemEntities.insert(mofem_ent);
    }
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
      //rval = moab.get_connectivity(edges,nodes,true); CHKERR_PETSC(rval);
      //use get adjacencies, this will allow take in account adjacencies set user
      rval = moab.get_adjacencies(edges,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
      {
	Range topo_nodes;
	rval = moab.get_connectivity(edges,topo_nodes,true); CHKERR_PETSC(rval);
	Range mid_nodes;
	rval = moab.get_connectivity(edges,mid_nodes,false); CHKERR_PETSC(rval);
	mid_nodes = subtract(mid_nodes,topo_nodes);
	nodes = subtract(nodes,mid_nodes);
      }
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
      SETERRQ(PETSC_COMM_SELF,1,"sorry, not implemented");
      break;
    case Hdiv:
      SETERRQ(PETSC_COMM_SELF,1,"sorry, not implemented, Hdiv not implemented for EDGEs");
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"add_ents_to_field_by_EDGEs this field not work for EDGEs");
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
      //rval = moab.get_connectivity(tris,nodes,true); CHKERR_PETSC(rval);
      //use get adjacencies, this will allow take in account adjacencies set user
      rval = moab.get_adjacencies(tris,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
      {
	Range topo_nodes;
	rval = moab.get_connectivity(tris,topo_nodes,true); CHKERR_PETSC(rval);
	Range mid_nodes;
	rval = moab.get_connectivity(tris,mid_nodes,false); CHKERR_PETSC(rval);
	mid_nodes = subtract(mid_nodes,topo_nodes);
	nodes = subtract(nodes,mid_nodes);
      }
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
PetscErrorCode FieldCore::add_ents_to_field_by_VERTICEs(const Range &nodes,const BitFieldId id,int verb) {
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
PetscErrorCode FieldCore::add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  Range nodes;
  rval = moab.get_entities_by_type(meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
  ierr = add_ents_to_field_by_VERTICEs(nodes,id,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_field_by_VERTICEs(const Range &nodes,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_VERTICEs(nodes,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
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
PetscErrorCode FieldCore::add_ents_to_field_by_TETs(const Range &tets,const BitFieldId id,int verb) {
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
      //rval = moab.get_connectivity(tets,nodes,true); CHKERR_PETSC(rval);
      //use get adjacencies, this will allow take in account adjacencies set user
      rval = moab.get_adjacencies(tets,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
      {
	Range topo_nodes;
	rval = moab.get_connectivity(tets,topo_nodes,true); CHKERR_PETSC(rval);
	Range mid_nodes;
	rval = moab.get_connectivity(tets,mid_nodes,false); CHKERR_PETSC(rval);
	mid_nodes = subtract(mid_nodes,topo_nodes);
	nodes = subtract(nodes,mid_nodes);
      }
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
PetscErrorCode FieldCore::add_ents_to_field_by_TETs(const EntityHandle meshset,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,true); CHKERR_PETSC(rval);
  ierr = add_ents_to_field_by_TETs(tets,id,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_field_by_TETs(const Range &tets,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TETs(tets,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
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
PetscErrorCode FieldCore::set_field_order(const Range &ents,const BitFieldId id,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
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
  Range ents_ = intersect(ents,ents_of_id_meshset);
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents for order change in the field %d\n",ents_.size());
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
  Range::iterator eit = ents_.begin();
  for(unsigned int ee = 0;ee<ents_.size();ee++,eit++) {
    //sanity check
    switch(miit->get_space()) {
      case H1:
	if(moab.type_from_handle(*eit)==MBVERTEX) {
	  if(order!=1) {
	    SETERRQ(PETSC_COMM_SELF,1,"approximation order for H1 sapce and vertex diffrent than 1 makes not sense"); 
	  }
	}
	break;
      case Hdiv:
	if(moab.type_from_handle(*eit)==MBVERTEX) {
	  SETERRQ(PETSC_COMM_SELF,1,"Hdiv space on vertices makes no sense"); 
	} 
	if(moab.type_from_handle(*eit)==MBEDGE) {
	  SETERRQ(PETSC_COMM_SELF,1,"Hdiv space on edges makes no sense"); 
	} 
	break;
      default:
	break;
    }
    //
    MoFEMEntity_multiIndex_ent_view::iterator miit3 = ents_id_view.find(*eit);
    if(miit3!=ents_id_view.end()) {
      const ApproximationOrder old_ApproximationOrder = (*miit3)->get_max_order();
      if(old_ApproximationOrder==order) continue;
      MoFEMEntity_multiIndex::iterator miit4 = entsMoabField.get<Unique_mi_tag>().find((*miit3)->get_unique_id());
      assert(miit4!=entsMoabField.end());
      if(miit4->get_max_order()<order) nb_ents_set_order_up++;
      if(miit4->get_max_order()>order) nb_ents_set_order_down++;
      typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type dof_set_type;
      dof_set_type& set_set = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>();
      dof_set_type::iterator miit5 = set_set.lower_bound(boost::make_tuple(miit4->get_name_ref(),miit4->get_ent()));
      dof_set_type::iterator hi_miit6 = set_set.upper_bound(boost::make_tuple(miit4->get_name_ref(),miit4->get_ent()));
      for(;miit5!=hi_miit6;miit5++) {
	if(miit5->get_dof_order()<=order) continue;
	bool success = dofsMoabField.modify(dofsMoabField.project<0>(miit5),DofMoFEMEntity_active_change(false));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      }
      bool success = entsMoabField.modify(entsMoabField.project<0>(miit4),MoFEMEntity_change_order(moab,order));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    } else {
      *(ApproximationOrder*)tag_data_order[ee] = order;
      RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*eit);
      if(miit_ref_ent==refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"database inconsistency");
      try { 
	MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
	//if(moabent.get_order_nb_dofs(moabent.get_max_order())==0) continue; 
	pair<MoFEMEntity_multiIndex::iterator,bool> e_miit = entsMoabField.insert(moabent);
	bool success = entsMoabField.modify(e_miit.first,MoFEMEntity_change_order(moab,order));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents for which order was increased %d (order %d)\n",nb_ents_set_order_up,order);
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents for which order was reduced %d (order %d)\n",nb_ents_set_order_down,order);
    PetscPrintf(PETSC_COMM_WORLD,"nb. of ents for which order set %d (order %d)\n",nb_ents_set_order_new,order);
  }
  if(verb>4) {
    list_ent_by_field_id(id);
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
  try{
    ierr = set_field_order(ents,id,order,verb); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
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
PetscErrorCode FieldCore::set_field_order(const Range &ents,const string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try{
    ierr = set_field_order(ents,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,id,order,verb); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::dofs_NoField(const BitFieldId id,map<EntityType,int> &dof_counter) {
  PetscFunctionBegin;
  //field it
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  //find fiels
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) SETERRQ(PETSC_COMM_SELF,1,"field not found");
  //serch if field meshset is in database
  RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(miit->meshset);
  if(miit_ref_ent==refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"database inconsistency");
  pair<MoFEMEntity_multiIndex::iterator,bool> e_miit;
  try {
    //create database entity
    MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
    e_miit = entsMoabField.insert(moabent);
    //this is nor real field in space (set order to zero)
    bool success = entsMoabField.modify(e_miit.first,MoFEMEntity_change_order(moab,0));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
      if(d_miit.second) {
	dof_counter[d_miit.first->get_ent_type()]++;
      }
      bool success = dofsMoabField.modify(d_miit.first,DofMoFEMEntity_active_change(true));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
PetscErrorCode FieldCore::dofs_L2H1HcurlHdiv(const BitFieldId id,map<EntityType,int> &dof_counter,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ents_by_ents;
  //find field
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) SETERRQ(PETSC_COMM_SELF,1,"field not found");
  //ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(miit->meshset,ents_of_id_meshset,false); CHKERR_PETSC(rval);
  //create dofsMoabField
  Range::iterator eit = ents_of_id_meshset.begin();
  for(;eit!=ents_of_id_meshset.end();eit++) {
    // check if ent is in ref meshset
    RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*eit);
    if(miit_ref_ent==refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) {
      SETERRQ(PETSC_COMM_SELF,1,"database inconsistency");
    }
    // create mofem entity linked to ref ent
    MoFEMEntity_multiIndex::iterator e_miit;
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
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
	      dof_counter[d_miit.first->get_ent_type()]++;
	      bool success = dofsMoabField.modify(d_miit.first,DofMoFEMEntity_active_change(true));
	      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
	    } 
	    //check ent
	    if(d_miit.first->get_ent()!=e_miit->get_ent()) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    }
	    if(d_miit.first->get_ent_type()!=e_miit->get_ent_type()) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    }
	    if(d_miit.first->get_id()!=e_miit->get_id()) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    }
	    //check dof
	    if(d_miit.first->get_dof_order()!=oo) {
	      ostringstream ss;
	      ss << "data inconsistency!" << endl;
	      ss << "should be " << mdof << endl;
	      ss << "but is " << *d_miit.first << endl;
	      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	    }
	    if(d_miit.first->get_max_order()!=e_miit->get_max_order()) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    }
	  } catch (const char* msg) {
	    SETERRQ(PETSC_COMM_SELF,1,msg);
	  }
	}
      }
    }
    if(DD != e_miit->get_max_rank()*e_miit->get_order_nb_dofs(e_miit->get_max_order())) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
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
    map<EntityType,int> dof_counter;
    if(verbose>0) {
      PetscPrintf(PETSC_COMM_WORLD,"Build Field %s\n",miit->get_name().c_str());
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
      int _dof_counter_ = 0;
      for(map<EntityType,int>::iterator it = dof_counter.begin();
	it!=dof_counter.end();it++) {
	switch (it->first) {
	  case MBVERTEX:
	    PetscPrintf(PETSC_COMM_WORLD,"nb added dofs (vertices) %d\n",it->second);
	  break;
	  case MBEDGE:
	    PetscPrintf(PETSC_COMM_WORLD,"nb added dofs (edges) %d\n",it->second);
	  break;
	  case MBTRI:
	    PetscPrintf(PETSC_COMM_WORLD,"nb added dofs (triangles) %d\n",it->second);
	  break;
	  case MBTET:
	    PetscPrintf(PETSC_COMM_WORLD,"nb added dofs (tets) %d\n",it->second);
	  break;
	  case MBENTITYSET:
	    PetscPrintf(PETSC_COMM_WORLD,"nb added dofs (meshsets) %d\n",it->second);
	  break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	}
	_dof_counter_ += it->second;
      }
      PetscPrintf(PETSC_COMM_WORLD,"nb added dofs %d\n",_dof_counter_);
    }
    if(verb>1) {
      list_ent_by_field_id(miit->get_id());
      list_dof_by_field_id(miit->get_id());
    }
  }
  PetscPrintf(PETSC_COMM_WORLD,"Nb. dofs %u\n",dofsMoabField.size());
  *build_MoFEM = 1<<0;
  //PetscFunctionReturn(0);
  return 0;
}
PetscErrorCode FieldCore::list_dof_by_field_name(const string &name) const {
  PetscFunctionBegin;
  PetscErrorCode _ierr;
  _ierr = list_dof_by_field_id(get_BitFieldId(name)); CHKERRQ(_ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::list_ent_by_field_name(const string &name) const {
  PetscFunctionBegin;
  PetscErrorCode _ierr;
  _ierr = list_ent_by_field_id(get_BitFieldId(name)); CHKERRQ(_ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::list_dof_by_field_id(const BitFieldId id) const {
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
PetscErrorCode FieldCore::list_ent_by_field_id(const BitFieldId id) const {
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
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"this < %s > is not there",MoFEMFiniteElement_name.c_str());
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_row_change_bit_add(get_BitFieldId(name_row)));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"this < %s > is not there",MoFEMFiniteElement_name.c_str());
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_row_change_bit_off(get_BitFieldId(name_row)));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
EntityHandle FieldCore::get_finite_element_meshset(const BitFEId id) const {
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
  const finiteElements_by_id& set = finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = set.find(id);
  if(miit==set.end()) THROW_AT_LINE("finite element not found");
  return miit->meshset;
}
EntityHandle FieldCore::get_finite_element_meshset(const string& name) const {	
  return get_finite_element_meshset(get_BitFEId(name));
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
PetscErrorCode FieldCore::add_problem(const string& name,enum MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  const moFEMProblems_by_name& set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator miit = set.find(name);
  if(miit==set.end()) {
    BitProblemId id = get_problem_shift();
    ierr = add_problem(id,name); CHKERRQ(ierr);
  } else if(bh == MF_EXCL) {
    SETERRQ1(PETSC_COMM_SELF,1,"problem is in database %s",name.c_str());
  }
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
  const EntityHandle idm = get_finite_element_meshset(id);
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
  const EntityHandle idm = get_finite_element_meshset(id);
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
  const EntityHandle idm = get_finite_element_meshset(id);
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
    idm = get_finite_element_meshset(id);
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
  const EntityHandle idm = get_finite_element_meshset(id);
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

PetscErrorCode FieldCore::add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const BitFEId id,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  EntityHandle idm = no_handle;
  try {
    idm = get_finite_element_meshset(id);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  Range prisms;
  rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,recursive); CHKERR_PETSC(rval);
  rval = moab.add_entities(idm,prisms); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_PRISMs(const Range& tets,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_finite_element_meshset(id);
  rval = moab.add_entities(idm,tets.subset_by_type(MBPRISM)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_PRISMs(const Range& prims,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_PRISMs(prims,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const string &name,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_PRISMs(meshset,get_BitFEId(name),recursive);  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const string &name,EntityType type,int verb) {
  PetscFunctionBegin;
  ierr = add_ents_to_finite_element_EntType_by_bit_ref(bit,BitRefLevel().set(),name,type,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const string &name,EntityType type,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);
  typedef RefMoFEMElement_multiIndex::index<EntType_mi_tag>::type refMoabFE_by_type;
  refMoabFE_by_type &ref_MoFEMFiniteElement = refinedMoFemElements.get<EntType_mi_tag>();
  refMoabFE_by_type::iterator miit = ref_MoFEMFiniteElement.lower_bound(type);
  refMoabFE_by_type::iterator hi_miit = ref_MoFEMFiniteElement.upper_bound(type);
  if(verb > 1) {
    PetscPrintf(PETSC_COMM_WORLD,"nb. ref elements in database %d\n",distance(miit,hi_miit));
  }
  int nb_add_FEs = 0;
  for(;miit!=hi_miit;miit++) {
    BitRefLevel bit2 = miit->get_BitRefLevel(); 
    //check if all bits in mask are ib fe bit2
    //if((miit->get_BitRefLevel()&bit)!=bit) continue;
    if((bit2&mask) != bit2) continue;
    if((bit2&bit).any()) {
      EntityHandle ent = miit->get_ref_ent();
      rval = moab.add_entities(idm,&ent,1); CHKERR_PETSC(rval);
      nb_add_FEs++;
    }
  }
  if(verb > 0) {
    ostringstream ss;
    ss << "Add Nb. FEs " << nb_add_FEs << " form BitRef " << bit << endl;
    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const string& name,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);
  if(recursive==false){
    rval = moab.add_entities(idm,&meshset,1); CHKERR_PETSC(rval);
  } else {
    Range meshsets;
    rval = moab.get_entities_by_type(meshset,MBENTITYSET,meshsets,false); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,meshsets); CHKERR_PETSC(rval);
  }
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
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::modify_problem_ref_level_set_bit(const string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name& set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator miit = set.find(name_problem);
  ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,1,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,problem_change_ref_level_bit_set(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::build_finite_element_data_dofs(EntMoFEMFiniteElement &EntFe,int verb) {
  PetscFunctionBegin;
  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  FEDofMoFEMEntity_multiIndex &data_dofs = const_cast<FEDofMoFEMEntity_multiIndex&>(EntFe.data_dofs);
  //clear data dofs multiindex //FIXME should be cleaned when dofs are cleaned form datasets
  data_dofs.clear();
  DofMoFEMEntity_multiIndex_active_view data_view;
  ierr = EntFe.get_MoFEMFiniteElement_data_dof_uid_view(dofsMoabField,data_view,Interface::UNION); CHKERRQ(ierr);
  DofMoFEMEntity_multiIndex_active_view::iterator viit_data,hi_viit_data;
  //loops over active dofs only
  viit_data = data_view.lower_bound(1); 
  hi_viit_data = data_view.upper_bound(1);
  for(;viit_data!=hi_viit_data;viit_data++) {
    try {
      switch((*viit_data)->get_space()) {
	case H1:
	case Hdiv:
	case Hcurl:
	case L2: 
	case NoField:
	{
	  SideNumber *side_number_ptr = EntFe.get_side_number_ptr(moab,(*viit_data)->get_ent());
	  FEDofMoFEMEntity FEDof(side_number_ptr,&**viit_data);
	  //add dofs to finite element multi_index database
	  pair<FEDofMoFEMEntity_multiIndex::iterator,bool> p;
	  p = data_dofs.insert(FEDof);
	  if(!p.second) {
	    //SETERRQ(PETSC_COMM_SELF,1,"insertion unsucessfull");
	  }
	}
	break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,1,msg);
    }
  }
  viit_data = data_view.lower_bound(1); 
  if(data_dofs.size()!=(unsigned int)distance(viit_data,hi_viit_data)) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::build_finite_element_uids_view(EntMoFEMFiniteElement &EntFe,int verb) {
  PetscFunctionBegin;
  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,1,"fields not build");
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_by_id;
  typedef MoFEMField_multiIndex::index<Meshset_mi_tag>::type field_by_meshset;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ent_by_ent;
  typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type dof_set_type;
  field_by_id &moabFields_by_id = moabFields.get<BitFieldId_mi_tag>();
  field_by_meshset &moabFields_by_meshset = moabFields.get<Meshset_mi_tag>();
  dof_set_type& dof_set = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>();
  EntityHandle fe_ent = EntFe.get_ent();
  //get id of mofem fields for row, col and data
  enum IntLoop { Row = 0,Col,Data,Last };
  BitFieldId FEAdj_fields[Last] = { 
    EntFe.get_BitFieldId_row(),
    EntFe.get_BitFieldId_col(),
    EntFe.get_BitFieldId_data() 
  };
  //get refinment level
  const BitRefLevel& bit_ref_MoFEMFiniteElement = EntFe.get_BitRefLevel();
  Range tets,faces,edges,nodes,meshsets,adj_ents,ent_ents;
  Range::iterator eit_eit;
  DofMoFEMEntity_multiIndex_uid_view* MoFEMFiniteElement_dof_uid_view[Last] = {
    &EntFe.row_dof_view, 
    &EntFe.col_dof_view, 
    &EntFe.data_dof_view
  };
  for(int ss = 0;ss<Last;ss++) {
    MoFEMFiniteElement_dof_uid_view[ss]->clear();
  }
  //lopp over all fields in database
  for(unsigned int ii = 0;ii<BitFieldId().size();ii++) {
    // common field id for Row, Col and Data
    BitFieldId id_common = 0;
    //check if the field (ii) is added to finite element
    for(int ss = 0;ss<Last;ss++) id_common |= FEAdj_fields[ss]&BitFieldId().set(ii);
    if( id_common.none() ) continue;
    //find in database data associated with the field (ii)
    field_by_id::iterator miit = moabFields_by_id.find(BitFieldId().set(ii));
    if(miit==moabFields_by_id.end()) {
      SETERRQ(PETSC_COMM_SELF,1,"data incosistency");
    }
    //get field (ii) space
    FieldSpace space = miit->get_space();
    //resolve entities on element, those entities are used to build
    //tag with dof uids on finite element tag
    switch (moab.type_from_handle(fe_ent)) {
      case MBVERTEX:
	switch (space) {
	  case H1: 
	    adj_ents.insert(fe_ent);
	    break;
	  case NoField: {
	    EntityHandle field_meshset = miit->get_meshset();
	    adj_ents.insert(field_meshset);
	  }
	  break;
      	  default:
  	   SETERRQ(PETSC_COMM_SELF,1,"this field is not implemented for VERTEX finite element");
	}
	break;
      case MBEDGE:
	switch (space) {
	  case H1: if(nodes.empty()) {
	      //moab.get_connectivity(&fe_ent,1,nodes,true);
	      //use get adjacencies, this will allow take in account adjacencies set user
	      rval = moab.get_adjacencies(&fe_ent,1,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
	      {
		Range topo_nodes;
		rval = moab.get_connectivity(&fe_ent,1,topo_nodes,true); CHKERR_PETSC(rval);
		Range mid_nodes;
		rval = moab.get_connectivity(&fe_ent,1,mid_nodes,false); CHKERR_PETSC(rval);
		mid_nodes = subtract(mid_nodes,topo_nodes);
		nodes = subtract(nodes,mid_nodes);
	      }
	    }
	    adj_ents.insert(nodes.begin(),nodes.end());
	    adj_ents.insert(fe_ent);
	    break;
	  case NoField: {
	    EntityHandle field_meshset = miit->get_meshset();
	    adj_ents.insert(field_meshset);
	  }
	  break;
      	  default:
  	   SETERRQ(PETSC_COMM_SELF,1,"this field is not implemented for EDGE finite element");
	}
	break;
      case MBTRI: 
	switch (space) {
	  case H1: 
	    //add nodes
	    if(nodes.empty()) {
	      //moab.get_connectivity(&fe_ent,1,nodes,true);
	      //use get adjacencies, this will allow take in account adjacencies set user
	      rval = moab.get_adjacencies(&fe_ent,1,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
	      {
		Range topo_nodes;
		rval = moab.get_connectivity(&fe_ent,1,topo_nodes,true); CHKERR_PETSC(rval);
		Range mid_nodes;
		rval = moab.get_connectivity(&fe_ent,1,mid_nodes,false); CHKERR_PETSC(rval);
		mid_nodes = subtract(mid_nodes,topo_nodes);
		nodes = subtract(nodes,mid_nodes);
	      }
	    }
	    adj_ents.insert(nodes.begin(),nodes.end());
	    //add edges
	    if(edges.empty()) {
	      moab.get_adjacencies(&fe_ent,1,1,false,edges);
	    }
	    adj_ents.insert(edges.begin(),edges.end());
	    for(Range::iterator eeit = edges.begin();eeit!=edges.end();eeit++) {
	      EntFe.get_side_number_ptr(moab,*eeit);
	    }
	    //add faces
	    adj_ents.insert(fe_ent);
	    break;
	    case NoField: {
	      EntityHandle field_meshset = miit->get_meshset();
	      adj_ents.insert(field_meshset);
	    }
	   break;
      	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"this field is not implemented for TRI finite element");
	}
	break;
      case MBTET:
	 switch (space) {
	  case H1: if(nodes.empty()) {
	    //moab.get_connectivity(&fe_ent,1,nodes,true);
	    //use get adjacencies, this will allow take in account adjacencies set user
	    rval = moab.get_adjacencies(&fe_ent,1,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
	    {
	      Range topo_nodes;
	      rval = moab.get_connectivity(&fe_ent,1,topo_nodes,true); CHKERR_PETSC(rval);
	      Range mid_nodes;
	      rval = moab.get_connectivity(&fe_ent,1,mid_nodes,false); CHKERR_PETSC(rval);
	      mid_nodes = subtract(mid_nodes,topo_nodes);
	      nodes = subtract(nodes,mid_nodes);
	    }
	    if(nodes.size()<4) {
	      SETERRQ(PETSC_COMM_SELF,1,"TET has at least 4 adjacent nodes; it can has more if user add more adjacencies");
	    }
	  }
  	   adj_ents.insert(nodes.begin(),nodes.end());
  	  case Hcurl: if(edges.empty()) moab.get_adjacencies(&fe_ent,1,1,false,edges);
  	   adj_ents.insert(edges.begin(),edges.end());
	   for(Range::iterator eeit = edges.begin();eeit!=edges.end();eeit++) {
	      EntFe.get_side_number_ptr(moab,*eeit);
	    }
  	  case Hdiv: if(faces.empty()) moab.get_adjacencies(&fe_ent,1,2,false,faces);
  	   adj_ents.insert(faces.begin(),faces.end());
	   for(Range::iterator fit = faces.begin();fit!=faces.end();fit++) {
	      EntFe.get_side_number_ptr(moab,*fit);
	    }
  	  case L2:
  	   adj_ents.insert(fe_ent);
  	   break;
	  case NoField: {
	      EntityHandle field_meshset = miit->get_meshset();
	      adj_ents.insert(field_meshset);
	   }
	   break;
      	  default:
  	   SETERRQ(PETSC_COMM_SELF,1,"this field is not implemented for TET finite element");
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
	      if(side_ptr->side_number!=ee) SETERRQ1(PETSC_COMM_SELF,1,"data insonsitency for edge %d",ee);
	      rval = moab.side_element(prism,1,6+ee,edge); CHKERR_PETSC(rval);
	      side_ptr = EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,edge);
	      if(side_ptr->side_number!=ee+6) {
		if(side_ptr->side_number!=ee) {
		  SETERRQ1(PETSC_COMM_SELF,1,"data insonsitency for edge %d",ee);
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
	      if(side_ptr->side_number!=nn) SETERRQ1(PETSC_COMM_SELF,1,"data insonsitency for node %d",nn);
	      rval = moab.side_element(prism,0,nn+3,node); CHKERR_PETSC(rval);
	      side_ptr = EntFe.get_RefMoFEMElement()->get_side_number_ptr(moab,node);
	      if(side_ptr->side_number!=nn+3) {
		if(side_ptr->side_number!=nn) {
		  SETERRQ1(PETSC_COMM_SELF,1,"data insonsitency for node %d",nn);
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
	    case H1: if(nodes.empty()) {
		//moab.get_connectivity(&fe_ent,1,nodes,true);
		//use get adjacencies, this will allow take in account adjacencies set user
		rval = moab.get_adjacencies(&fe_ent,1,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
		{
		  Range topo_nodes;
		  rval = moab.get_connectivity(&fe_ent,1,topo_nodes,true); CHKERR_PETSC(rval);
		  Range mid_nodes;
		  rval = moab.get_connectivity(&fe_ent,1,mid_nodes,false); CHKERR_PETSC(rval);
		  mid_nodes = subtract(mid_nodes,topo_nodes);
		  nodes = subtract(nodes,mid_nodes);
		}
	      }
	      adj_ents.insert(nodes.begin(),nodes.end());
	    case Hcurl: {
	      SideNumber_multiIndex::nth_index<2>::type::iterator
		siit = side_table.get<2>().lower_bound(MBEDGE), hi_siit = side_table.get<2>().upper_bound(MBEDGE);
	      for(;siit!=hi_siit;siit++) adj_ents.insert(siit->ent);
	    }
	    case Hdiv: {
	      SideNumber_multiIndex::nth_index<2>::type::iterator
		siit = side_table.get<2>().lower_bound(MBTRI), hi_siit = side_table.get<2>().upper_bound(MBTRI);
	      for(;siit!=hi_siit;siit++) adj_ents.insert(siit->ent);
	    }
	    case L2:
	      adj_ents.insert(fe_ent); 
	      break;
	    case NoField: {
	      EntityHandle field_meshset = miit->get_meshset();
	      adj_ents.insert(field_meshset);
	    }
	    break;
	    default:
	      SETERRQ(PETSC_COMM_SELF,1,"this field is not implemented for PRISM finite element");
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
    //loop over adjacent to finite entities, and find dofs on those entities
    //this part is to build MoFEMFiniteElement_dof_uid_view
    Range::iterator eit2 = adj_ents.begin();
    for(;eit2!=adj_ents.end();eit2++) {
      ref_ent_by_ent::iterator ref_ent_miit = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*eit2);
      if(ref_ent_miit==refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) {
	SETERRQ(PETSC_COMM_SELF,1,"ref ent not in database"); 
      }
      const BitRefLevel& bit_ref_ent = ref_ent_miit->get_BitRefLevel();
      if(!(bit_ref_MoFEMFiniteElement&bit_ref_ent).any()) {
	ostringstream ss;
	ss << "top tip: check if you seed mesh with the elements for bit ref level1" << endl;
	ss << "inconsitency in database entity" << " type " 
	  << moab.type_from_handle(*eit2) << " bits ENT " << bit_ref_ent << endl;
	ss << "inconsitency in database entity" << " type " 
	  << moab.type_from_handle(EntFe.get_ent()) << " bits FE  " << bit_ref_MoFEMFiniteElement << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      dof_set_type::iterator ents_miit2 = dof_set.lower_bound(boost::make_tuple(miit->get_name_ref(),ref_ent_miit->get_ref_ent()));
      dof_set_type::iterator ents_hi_miit2 = dof_set.upper_bound(boost::make_tuple(miit->get_name_ref(),ref_ent_miit->get_ref_ent()));
      for(int ss = 0;ss<Last;ss++) {
	if( !(FEAdj_fields[ss].test(ii)) ) continue;
	dof_set_type::iterator ents_miit3 = ents_miit2;
	for(;ents_miit3!=ents_hi_miit2;ents_miit3++) {
	  MoFEMFiniteElement_dof_uid_view[ss]->insert(&*ents_miit3);
	}
      }
    }
  }
  if(verb>2) {
    ostringstream ss;
    ss << "add: FE data"  << endl << EntFe << endl;
    //rows
    DofMoFEMEntity_multiIndex_uid_view MoFEMFiniteElement_row_dof_uid_view;
    ierr = EntFe.get_MoFEMFiniteElement_row_dof_uid_view(dofsMoabField,MoFEMFiniteElement_row_dof_uid_view); CHKERRQ(ierr);
    DofMoFEMEntity_multiIndex_uid_view::iterator miit_row = MoFEMFiniteElement_row_dof_uid_view.begin();
    ss << "rows dofs" << endl;
    for(;miit_row!=MoFEMFiniteElement_row_dof_uid_view.end();miit_row++) ss << **miit_row << endl;
    //cols
    DofMoFEMEntity_multiIndex_uid_view MoFEMFiniteElement_col_dof_uid_view;
    ierr = EntFe.get_MoFEMFiniteElement_col_dof_uid_view(dofsMoabField,MoFEMFiniteElement_col_dof_uid_view); CHKERRQ(ierr);
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
    EntityHandle meshset = get_finite_element_meshset(MoFEMFiniteElement_miit->get_id());
    // get entities from finite element meshset // if meshset 
    Range MoFEMFiniteElement_ents;
    rval = moab.get_entities_by_handle(meshset,MoFEMFiniteElement_ents,false); CHKERR_PETSC(rval);
    //loop meshset Ents and add finite elements
    Range::iterator eit = MoFEMFiniteElement_ents.begin();
    for(;eit!=MoFEMFiniteElement_ents.end();eit++) {
      // note: iterator is a wrapper
      // check if is in refinedMoFemElements database
      ref_MoFEMFiniteElement_by_ent::iterator ref_MoFEMFiniteElement_miit; 
      ref_MoFEMFiniteElement_miit = refinedMoFemElements.get<MoABEnt_mi_tag>().find(*eit);
      if(ref_MoFEMFiniteElement_miit == refinedMoFemElements.get<MoABEnt_mi_tag>().end()) {
	ostringstream ss;
	ss << "ref MoFEMFiniteElement not in database ent = " << *eit;
	ss << " type " << moab.type_from_handle(*eit);
	ss << " " << *MoFEMFiniteElement_miit;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      EntMoFEMFiniteElement EntFe(moab,ref_MoFEMFiniteElement_miit->get_RefMoFEMElement(),&*MoFEMFiniteElement_miit);
      pair<EntMoFEMFiniteElement_multiIndex::iterator,bool> p = finiteElementsMoFEMEnts.insert(EntFe);
      ierr = build_finite_element_uids_view(const_cast<EntMoFEMFiniteElement&>(*p.first),verb); CHKERRQ(ierr);
      ierr = build_finite_element_data_dofs(const_cast<EntMoFEMFiniteElement&>(*p.first),verb); CHKERRQ(ierr);
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
      int count = distance(miit,hi_miit);
      ostringstream ss;
      ss << *id_MoFEMFiniteElement << " Nb. FEs " << count << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
  }
  *build_MoFEM |= 1<<1;
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::build_adjacencies(const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,1,"field not build");
  if(!(*build_MoFEM)&(1<<1)) SETERRQ(PETSC_COMM_SELF,1,"fe not build");
  typedef MoFEMEntity_multiIndex::index<Unique_mi_tag>::type ents_by_uid;
  EntMoFEMFiniteElement_multiIndex::iterator fit = finiteElementsMoFEMEnts.begin();
  for(;fit!=finiteElementsMoFEMEnts.end();fit++) {
    if(!ents.empty()) {
      if(ents.find(fit->get_ent())==ents.end()) {
	continue;
      }
    }
    UId ent_uid;
    ent_uid = 0;
    DofMoFEMEntity_multiIndex_uid_view::iterator rvit;
    rvit = fit->row_dof_view.begin();
    for(;rvit!=fit->row_dof_view.end();rvit++) {
      if( ent_uid == (*rvit)->get_MoFEMEntity_ptr()->get_unique_id()) continue;
      ent_uid = (*rvit)->get_MoFEMEntity_ptr()->get_unique_id();
      pair<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
      p = entFEAdjacencies.insert(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap((*rvit)->get_MoFEMEntity_ptr(),&*fit));
      bool success = entFEAdjacencies.modify(
	p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_by_what(by_row));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    }
    ent_uid = 0;
    DofMoFEMEntity_multiIndex_uid_view::iterator cvit;
    cvit = fit->col_dof_view.begin();
    for(;cvit!=fit->col_dof_view.end();cvit++) {
      if( ent_uid == (*cvit)->get_MoFEMEntity_ptr()->get_unique_id()) continue;
      ent_uid = (*cvit)->get_MoFEMEntity_ptr()->get_unique_id();
      pair<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
      p = entFEAdjacencies.insert(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap((*cvit)->get_MoFEMEntity_ptr(),&*fit));
      bool success = entFEAdjacencies.modify(
	p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_by_what(by_col));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    }
    ent_uid = 0;
    DofMoFEMEntity_multiIndex_uid_view::iterator dvit;
    dvit = fit->data_dof_view.begin();
    for(;dvit!=fit->data_dof_view.end();dvit++) {
      if( ent_uid == (*dvit)->get_MoFEMEntity_ptr()->get_unique_id()) continue;
      ent_uid = (*dvit)->get_MoFEMEntity_ptr()->get_unique_id();
      pair<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
      p = entFEAdjacencies.insert(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap((*dvit)->get_MoFEMEntity_ptr(),&*fit));
      bool success = entFEAdjacencies.modify(
	p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_by_what(by_data));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
PetscErrorCode FieldCore::build_adjacencies(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range ents;
  ierr = get_entities_by_ref_level(bit,mask,ents); CHKERRQ(ierr);
  ierr = build_adjacencies(ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::build_adjacencies(const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = build_adjacencies(bit,BitRefLevel().set(),verb); CHKERRQ(ierr);
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
  //iterate problems
  for(;p_miit!=moFEMProblems.end();p_miit++) {
    if(p_miit->get_BitRefLevel().none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",p_miit->get_name().c_str());
    }
    //zero finite elements
    bool success = moFEMProblems.modify(p_miit,problem_clear_numered_finite_elements_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    //miit2 iterator for finite elements
    EntMoFEMFiniteElement_multiIndex::iterator miit2 = finiteElementsMoFEMEnts.begin();
    EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finiteElementsMoFEMEnts.end();
    DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
    EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
    //iterate all finite elemen entities in database
    for(;miit3!=hi_miit2;miit3++) {
      //if element is in problem
      if((miit3->get_id()&p_miit->get_BitFEId()).any()) {
	//if finite element bit level has all refined bits sets
	if((miit3->get_BitRefLevel()&p_miit->get_BitRefLevel())==p_miit->get_BitRefLevel()) {
	  //get dof uids for rows and columns
	  ierr = miit3->get_MoFEMFiniteElement_row_dof_uid_view(dofsMoabField,dofs_rows); CHKERRQ(ierr);
	  ierr = miit3->get_MoFEMFiniteElement_col_dof_uid_view(dofsMoabField,dofs_cols); CHKERRQ(ierr);
	}
      }
    }
    //zero rows
    success = moFEMProblems.modify(p_miit,problem_zero_nb_rows_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    //zero cols
    success = moFEMProblems.modify(p_miit,problem_zero_nb_cols_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    //add dofs for rows
    DofMoFEMEntity_multiIndex_active_view::iterator miit4,hi_miit4;
    miit4 = dofs_rows.lower_bound(1);
    hi_miit4 = dofs_rows.upper_bound(1);
    for(;miit4!=hi_miit4;miit4++) {
      success = moFEMProblems.modify(p_miit,problem_row_change(&**miit4));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    }
    //add dofs for cols
    DofMoFEMEntity_multiIndex_active_view::iterator miit5,hi_miit5;
    miit5 = dofs_cols.lower_bound(1);
    hi_miit5 = dofs_cols.upper_bound(1);
    for(;miit5!=hi_miit5;miit5++) {
      success = moFEMProblems.modify(p_miit,problem_col_change(&**miit5));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    }
    //number dofs on rows and collumns
    success = moFEMProblems.modify(p_miit,problem_row_number_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    success = moFEMProblems.modify(p_miit,problem_col_number_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    //job done, some debugging and postprocessing
    if(verbose>0) {
      PetscPrintf(PETSC_COMM_WORLD,"Problem %s Nb. rows %u Nb. cols %u\n",
	p_miit->get_name().c_str(),
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
PetscErrorCode FieldCore::clear_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  MoFEMProblem_multiIndex::iterator p_miit = moFEMProblems.begin();
  //iterate problems
  for(;p_miit!=moFEMProblems.end();p_miit++) {
    //zero rows
    bool success = moFEMProblems.modify(p_miit,problem_zero_nb_rows_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    //zero cols
    success = moFEMProblems.modify(p_miit,problem_zero_nb_cols_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    //clear finite elements
    success = moFEMProblems.modify(p_miit,problem_clear_numered_finite_elements_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
  }
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
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Simple partition problem %s\n",name.c_str());
  }
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
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      if(miit_row->part == pcomm->rank()) {
	success = dofs_row_by_idx.modify(miit_row,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      if(miit_col->part == pcomm->rank()) {
	success = dofs_col_by_idx.modify(miit_col,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"Partition problem %s\n",name.c_str());
  }
  typedef NumeredDofMoFEMEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",name.c_str());
  DofIdx nb_dofs_row = p_miit->get_nb_dofs_row();
  Mat Adj;
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"\tcreate Adj matrix\n");
  }
  try {
    ierr = partition_create_Mat<Idx_mi_tag>(name,&Adj,MATMPIADJ,true,verb); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"\t<- done\n");
  }
  //PetscBarrier(PETSC_NULL);
  PetscInt m,n;
  ierr = MatGetSize(Adj,&m,&n); CHKERRQ(ierr);
  if(m!=p_miit->get_nb_dofs_row()) {
    SETERRQ3(PETSC_COMM_SELF,1,"row number inconsistency %d != %d nb cols. %d",m,p_miit->get_nb_dofs_row(),p_miit->get_nb_dofs_col());
  }
  if(n!=p_miit->get_nb_dofs_col()) {
    SETERRQ(PETSC_COMM_SELF,1,"col number inconsistency");
  }
  if(m != n) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  if(verb>2) {
    MatView(Adj,PETSC_VIEWER_STDOUT_WORLD);
  }
  //partitioning
  MatPartitioning part;
  IS is;
  ierr = MatPartitioningCreate(MPI_COMM_WORLD,&part); CHKERRQ(ierr);
  ierr = MatPartitioningSetAdjacency(part,Adj); CHKERRQ(ierr);
  ierr = MatPartitioningSetFromOptions(part); CHKERRQ(ierr);
  ierr = PetscBarrier(PETSC_NULL); CHKERRQ(ierr);
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
  if(size_is_gather != (int)nb_dofs_row) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  ISGetSize(is_num,&size_is_num);
  if(size_is_num != (int)nb_dofs_row) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
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
  try {
  for(;miit_dofs_row!=dofs_row_by_idx_no_const.end();miit_dofs_row++,miit_dofs_col++) {
    if(miit_dofs_col==dofs_col_by_idx_no_const.end()) {
      SETERRQ(PETSC_COMM_SELF,1,"check finite element definition, nb. of rows is not equal to number for columns");
    }
    if(miit_dofs_row->get_unique_id()!=miit_dofs_col->get_unique_id()) {
      SETERRQ(PETSC_COMM_SELF,1,"check finite element definition, nb. of rows is not equal to columns");
    }
    if(miit_dofs_row->dof_idx!=miit_dofs_col->dof_idx) {
      SETERRQ(PETSC_COMM_SELF,1,"check finite element definition, nb. of rows is not equal to columns");
    }
    assert(petsc_idx[miit_dofs_row->dof_idx]>=0);
    assert(petsc_idx[miit_dofs_row->dof_idx]<(int)p_miit->get_nb_dofs_row());
    bool success = dofs_row_by_idx_no_const.modify(miit_dofs_row,NumeredDofMoFEMEntity_part_change(part_number[miit_dofs_row->dof_idx],petsc_idx[miit_dofs_row->dof_idx]));
    if(!success) {
      SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    }
    success = dofs_col_by_idx_no_const.modify(miit_dofs_col,NumeredDofMoFEMEntity_part_change(part_number[miit_dofs_col->dof_idx],petsc_idx[miit_dofs_col->dof_idx]));
    if(!success) {
      SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    }
    if(miit_dofs_row->part == pcomm->rank()) {
      assert(miit_dofs_row->part==miit_dofs_col->part);
      assert(miit_dofs_row->petsc_gloabl_dof_idx==miit_dofs_col->petsc_gloabl_dof_idx);
      success = dofs_row_by_idx_no_const.modify(miit_dofs_row,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
      if(!success) {
	SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      }
      success = dofs_col_by_idx_no_const.modify(miit_dofs_col,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
      if(!success) {
	SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      }
    }
  }
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
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
  //PetscBarrier(PETSC_NULL);
  if(debug>0) {
    try {
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
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,1,msg);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
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
    PetscPrintf(PETSC_COMM_WORLD,"Compose problem %s from rows of %s and columns of %s\n",
      p_miit->get_name().c_str(),problem_for_rows.c_str(),problem_for_cols.c_str());
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
	  if(!(adj_miit->by_other&by_row)) {
	    // if it is not row if element
	    continue; 
	  }
	  if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) {
	    // if element is not part of prblem
	    continue; 
	  }
	  if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&miit_row->get_BitRefLevel()).none()) {
	    // if entity is not problem refinment level
	    continue; 
	  }
	  NumeredDofMoFEMEntity_multiIndex_uid_view row_dof_view;
	  ierr = adj_miit->EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_row_dof_uid_view( 
	    dofs_row,row_dof_view,Interface::UNION); CHKERRQ(ierr);
	  NumeredDofMoFEMEntity_multiIndex_uid_view::iterator rdvit;
	  rdvit = row_dof_view.begin();
	  for(;rdvit!=row_dof_view.end();rdvit++) {
	    DofIdx petsc_global_idx = (*rdvit)->get_petsc_gloabl_dof_idx();
	    rows_problem_map[petsc_global_idx] = (*rdvit);
	  }
	}
      }
    }
    if(rows_problem_map.size() != (unsigned int)p_miit->get_nb_dofs_row()) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
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
	  if(!(adj_miit->by_other&by_col)) {
	    // if it is not row if element
	    continue; 
	  }
	  if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) {
	    // if element is not part of prblem
	    continue; 
	  }
	  if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&miit_col->get_BitRefLevel()).none()) {
	    // if entity is not problem refinment level
	    continue; 
	  }
	  NumeredDofMoFEMEntity_multiIndex_uid_view col_dof_view;
	  ierr = adj_miit->EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_col_dof_uid_view( 
	    dofs_col,col_dof_view,Interface::UNION); CHKERRQ(ierr);
	  NumeredDofMoFEMEntity_multiIndex_uid_view::iterator cdvit;
	  cdvit = col_dof_view.begin();
	  for(;cdvit!=col_dof_view.end();cdvit++) {
	    DofIdx petsc_global_idx = (*cdvit)->get_petsc_gloabl_dof_idx();
	    cols_problem_map[petsc_global_idx] = (*cdvit);
	  }
	}
      }
    }
    if(cols_problem_map.size() != (unsigned int)p_miit->get_nb_dofs_col()) {
      SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency, %u != %u",cols_problem_map.size(),(unsigned int)p_miit->get_nb_dofs_col());
    }
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
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    if(pr_dof->get_part() == pcomm->rank()) {
      success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    }
  }
  } else {
    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(p_miit_row,diit)) {
      bool success = moFEMProblems.modify(moFEMProblems.project<0>(p_miit),problem_row_change(diit->get_DofMoFEMEntity_ptr()));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      int part_number = diit->get_part();
      int petsc_global_dof = diit->get_petsc_gloabl_dof_idx();
      NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_row_by_uid.find(diit->get_unique_id());
      if(pr_dof == dofs_row_by_uid.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      if(pr_dof->get_part() == pcomm->rank()) {
	success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    if(pr_dof->get_part() == pcomm->rank()) {
      success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
    }
  }
  } else {
    for(_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(p_miit_col,diit)) {
      bool success = moFEMProblems.modify(moFEMProblems.project<0>(p_miit),problem_col_change(diit->get_DofMoFEMEntity_ptr()));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      int part_number = diit->get_part();
      int petsc_global_dof = diit->get_petsc_gloabl_dof_idx();
      NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_col_by_uid.find(diit->get_unique_id());
      success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      if(pr_dof->get_part() == pcomm->rank()) {
	success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
    if((miit3->get_id()&p_miit->get_BitFEId()).none()) continue; // if element is not part of prblem
    if((miit3->get_BitRefLevel()&p_miit->get_BitRefLevel())!=p_miit->get_BitRefLevel()) continue; // if entity is not problem refinment level
    {
      NumeredDofMoFEMEntity_multiIndex_uid_view rows_view,cols_view;
      //rows_view
      ierr = miit3->get_MoFEMFiniteElement_row_dof_uid_view(p_miit->numered_dofs_rows,rows_view,Interface::UNION); CHKERRQ(ierr);
      if(rows_view.empty()) continue;
      //cols_vies
      ierr = miit3->get_MoFEMFiniteElement_col_dof_uid_view(p_miit->numered_dofs_cols,cols_view,Interface::UNION); CHKERRQ(ierr);
      if(cols_view.empty()) continue;
      pair<NumeredMoFEMFiniteElement_multiIndex::iterator,bool> p;
      p = numeredFiniteElements.insert(NumeredMoFEMFiniteElement(&*miit3));
      if(!p.second) {
	SETERRQ(PETSC_COMM_SELF,1,"element is there");
      }
      //rows element dof multiindices
      FENumeredDofMoFEMEntity_multiIndex &rows_dofs = const_cast<FENumeredDofMoFEMEntity_multiIndex&>(p.first->rows_dofs);
      rows_dofs.clear();
      NumeredDofMoFEMEntity_multiIndex_uid_view::iterator viit_rows = rows_view.begin();
      vector<int> parts(pcomm->size(),0);
      for(;viit_rows!=rows_view.end();viit_rows++) {
	try {
	  SideNumber *side_number_ptr = p.first->get_side_number_ptr(moab,(*viit_rows)->get_ent());
	  FENumeredDofMoFEMEntity FEDof(side_number_ptr,&**viit_rows);
	  pair<FENumeredDofMoFEMEntity_multiIndex::iterator,bool> pp;
	  pp = rows_dofs.insert(FEDof);
	  if(!p.second) {
	    SETERRQ(PETSC_COMM_SELF,1,"element is there");
	  }
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
	}
	parts[(*viit_rows)->part]++;
      }
      vector<int>::iterator pos = max_element(parts.begin(),parts.end());
      unsigned int max_part = distance(parts.begin(),pos);
      bool success = numeredFiniteElements.modify(p.first,NumeredMoFEMFiniteElement_change_part(max_part));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
      if(do_skip) if(max_part!=pcomm->rank()) continue; 
      //cols elmeny dof multiindices
      FENumeredDofMoFEMEntity_multiIndex &cols_dofs = const_cast<FENumeredDofMoFEMEntity_multiIndex&>(p.first->cols_dofs);
      cols_dofs.clear();
      NumeredDofMoFEMEntity_multiIndex_uid_view::iterator viit_cols = cols_view.begin();
      for(;viit_cols!=cols_view.end();viit_cols++) {
	try {
	  SideNumber *side_number_ptr = p.first->get_side_number_ptr(moab,(*viit_cols)->get_ent());
	  FENumeredDofMoFEMEntity FEDof(side_number_ptr,&**viit_cols);
	  pair<FENumeredDofMoFEMEntity_multiIndex::iterator,bool> pp;
	  pp = cols_dofs.insert(FEDof);
	  if(!p.second) {
	    SETERRQ(PETSC_COMM_SELF,1,"element is there");
	  }
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
	}
      }
      if(verb>1) {
	ostringstream ss;
	ss << *p_miit << endl;
	ss << *p.first << endl;
	typedef FENumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type FENumeredDofMoFEMEntity_multiIndex_by_Unique_mi_tag;
	FENumeredDofMoFEMEntity_multiIndex_by_Unique_mi_tag::iterator miit = p.first->rows_dofs.get<Unique_mi_tag>().begin();
	for(;miit!= p.first->rows_dofs.get<Unique_mi_tag>().end();miit++) ss << "rows: " << *miit << endl;
	miit = p.first->cols_dofs.get<Unique_mi_tag>().begin();
	for(;miit!=p.first->cols_dofs.get<Unique_mi_tag>().end();miit++) ss << "cols: " << *miit << endl;
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	PetscSynchronizedFlush(PETSC_COMM_WORLD); 
      }
    }
  }
  if(verb>0) {
    typedef NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_Part_mi_tag>::type NumeredMoFEMFiniteElement_multiIndex_by_part;
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator MoFEMFiniteElement_miit = numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().lower_bound(pcomm->rank());
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator hi_MoMoFEMFiniteElement_miitFEMFE_miit = numeredFiniteElements.get<MoFEMFiniteElement_Part_mi_tag>().upper_bound(pcomm->rank());
    int count = distance(MoFEMFiniteElement_miit,hi_MoMoFEMFiniteElement_miitFEMFE_miit);
    ostringstream ss;
    ss << *p_miit;
    ss << " Nb. elems " << count << " on proc " << pcomm->rank() << endl;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }
  *build_MoFEM |= 1<<5;  
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::partition_check_matrix_fill_in(const string &problem_name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  struct test_matrix_fill_in:public FEMethod {
    Interface& moab;
    Mat A;
    PetscErrorCode ierr;
    ErrorCode rval;

    test_matrix_fill_in(Interface& _moab,Mat _A): moab(_moab),A(_A) {};
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(refinedMoFemElements->find(fe_ptr->get_ent())==refinedMoFemElements->end()) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      FENumeredDofMoFEMEntity_multiIndex::iterator rit = row_multiIndex->begin();
      for(;rit!=row_multiIndex->end();rit++) {
	if(refinedMoFemEntities->find(rit->get_ent())==refinedMoFemEntities->end()) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	if(!rit->get_active()) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Composite_unique_mi_tag>::type::iterator ait;
	ait = fem_adjacencies->get<Composite_unique_mi_tag>().find(boost::make_tuple(
	      rit->get_MoFEMEntity_ptr()->get_unique_id(),fe_ptr->get_unique_id()));
	if(ait==fem_adjacencies->end()) {
	    ostringstream ss;
	    ss << *rit << endl;
	    ss << *fe_ptr << endl;
	    ss << "dof: " << rit->get_BitRefLevel() << endl;
	    ss << "fe: " << fe_ptr->get_BitRefLevel() << endl;
	    ss << "problem: " << problem_ptr->get_BitRefLevel() << endl;
	    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	    SETERRQ(PETSC_COMM_SELF,1,"adjacencies data inconsistency");
	} else {
	  UId uid = ait->get_ent_unique_id();
	  if(ents_moabfield->find(uid) == ents_moabfield->end()) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  } 
	  if(dofs_moabfield->find(rit->get_unique_id())==dofs_moabfield->end()) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	}
	int row = rit->get_petsc_gloabl_dof_idx();
	FENumeredDofMoFEMEntity_multiIndex::iterator cit = col_multiIndex->begin();
	for(;cit!=col_multiIndex->end();cit++) {
	  if(refinedMoFemEntities->find(cit->get_ent())==refinedMoFemEntities->end()) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  if(!cit->get_active()) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  int col = cit->get_petsc_gloabl_dof_idx();
	  ait = fem_adjacencies->get<Composite_unique_mi_tag>().find(boost::make_tuple(
	      cit->get_MoFEMEntity_ptr()->get_unique_id(),fe_ptr->get_unique_id()));
	  if(ait==fem_adjacencies->end()) {
	    SETERRQ(PETSC_COMM_SELF,1,"adjacencies data inconsistency");
	  } else {
	    UId uid = ait->get_ent_unique_id();
	    if(ents_moabfield->find(uid) == ents_moabfield->end()) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    } 
	    if(dofs_moabfield->find(cit->get_unique_id())==dofs_moabfield->end()) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    }
	  }
  	  ierr = MatSetValue(A,row,col,1,INSERT_VALUES);

	  if(ierr!=0) {
	  //if(row == 87 && col == 909) {
	  
	    EntityHandle ent = fe_ptr->get_ent();
      
	    ostringstream ss;
	    ss << "fe:\n" << *fe_ptr << endl;
	    ss << "row:\n" << *rit << endl;
	    ss << "col:\n" << *cit << endl;

	    ss << "fe:\n" << fe_ptr->get_BitRefLevel() << endl;
	    ss << "row:\n" << rit->get_BitRefLevel() << endl;
	    ss << "col:\n" << cit->get_BitRefLevel() << endl;

	    ss << "edges:\n";
	    for(int ee = 0;ee<6;ee++) {
	      EntityHandle edge;
	      rval = moab.side_element(ent,1,ee,edge); CHKERR_THROW(rval);
	      ss << edge << " ";
	    }
	    ss << endl;

	    PetscPrintf(PETSC_COMM_WORLD,"%s\n",ss.str().c_str());
	  //}
	  }
	}
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


  };

  Mat A;
  ierr = MatCreateMPIAIJWithArrays(problem_name,&A); CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);  CHKERRQ(ierr);
  test_matrix_fill_in method(moab,A);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",problem_name.c_str());
  }
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"check problem < %s >\n",problem_name.c_str());
  }
  //MoFEMFiniteElement set
  MoFEMFiniteElement_multiIndex::iterator fe = finiteElements.begin();
  MoFEMFiniteElement_multiIndex::iterator hi_fe = finiteElements.end();
  for(;fe!=hi_fe;fe++) {
    if(verb>0) {
      PetscPrintf(PETSC_COMM_WORLD,"\tcheck element %s\n",fe->get_name().c_str());
    }

    ierr = loop_finite_elements(problem_name,fe->get_name(),method,0,pcomm->size(),verb);  CHKERRQ(ierr);

  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatDestroy(&A); CHKERRQ(ierr);
 
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
      if(fe_it->rows_dofs.empty()) SETERRQ(PETSC_COMM_SELF,1,"no row dofs on this element why is it here?");
      if(fe_it->cols_dofs.empty()) SETERRQ(PETSC_COMM_SELF,1,"no col dofs on this element why is it here?");
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
        if(diit->petsc_local_dof_idx!=(DofIdx)-1) SETERRQ(PETSC_COMM_SELF,1,"inconsistent data, ghost dof already set");
        bool success = dof_by_uid_no_const[ss]->modify(diit,NumeredDofMoFEMEntity_local_idx_change(nb_local_dofs[ss]++));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
    if(eiit == refinedMoFemEntities.get<MoABEnt_mi_tag>().end())  {
      SETERRQ(PETSC_COMM_SELF,1,"entity is not in database");
    }
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
      case MBTET: 
	p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*eiit)));	
	break;
      case MBPRISM: 
	p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*eiit)));	
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::seed_ref_level_2D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  Range ents2d;
  rval = moab.get_entities_by_dimension(meshset,2,ents2d,false); CHKERR_PETSC(rval);
  ierr = seed_ref_level_2D(ents2d,bit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::seed_ref_level_2D(const Range &ents2d,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin; 
  if(verb==-1) verb = verbose;
  try {
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
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
      Range ents;
      rval = moab.get_adjacencies(ents2d,dd,true,ents,Interface::UNION); CHKERR_PETSC(rval);
      Range::iterator eit = ents.begin();
      for(;eit!=ents.end();eit++) {
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*eit));
        if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
	  bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
PetscErrorCode FieldCore::seed_ref_level_3D(const Range &ents3d,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  try {
    if(verb > 1) {
      PetscPrintf(PETSC_COMM_WORLD,"nb. 3d entities for seed %d\n",ents3d.size());
    }
    Range::iterator tit = ents3d.begin();
    for(;tit!=ents3d.end();tit++) {
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*tit));
      if(debug > 0) {
	ierr = test_moab(moab,*tit); CHKERRQ(ierr);
      }
      bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
    for(int dd = 0;dd<3;dd++) {
      Range ents;
      rval = moab.get_adjacencies(ents3d.subset_by_type(MBTET),dd,true,ents,Interface::UNION); CHKERR_PETSC(rval);
      Range::iterator eit = ents.begin();
      for(;eit!=ents.end();eit++) {
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*eit));
	bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
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
  Range ents3d;
  rval = moab.get_entities_by_dimension(meshset,3,ents3d,false); CHKERR_PETSC(rval);
  ierr = seed_ref_level_3D(ents3d,bit,verb); CHKERRQ(ierr);
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
PetscErrorCode FieldCore::shift_left_bit_ref(const int shift,int verb) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::shift_right_bit_ref(const int shift,int verb) {
  PetscFunctionBegin;	
  if(verb==-1) verb = verbose;
  BitRefLevel delete_bits;
  for(int ii = 0;ii<shift;ii++) {
    delete_bits.set(0);
    ierr = delete_ents_by_bit_ref(delete_bits,delete_bits,verb); CHKERRQ(ierr);
  }
  RefMoFEMEntity_multiIndex::iterator ent_it = refinedMoFemEntities.begin();
  for(;ent_it!=refinedMoFemEntities.end();ent_it++) {
    if(verb>5) {
      cout << ent_it->get_BitRefLevel() << " : ";
    }
    bool success = refinedMoFemEntities.modify(ent_it,RefMoFEMEntity_change_right_shift(shift));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsistency in data");
    if(verb>5) {
      cout << ent_it->get_BitRefLevel() << endl;
    }
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb) {
  PetscFunctionBegin;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  rval = moab.add_entities(meshset,ents); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = moab.get_entities_by_type(0,type,ents,false); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();) {
    BitRefLevel bit2;
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
    if(mask.any()&&bit2.none()) {
      eit = ents.erase(eit);
      continue;
    }
    if((bit2&mask) != bit2) {
      eit = ents.erase(eit);
      continue;
    }
    if((bit2&bit).none()) {
      eit = ents.erase(eit);
      continue;
    }
    eit++;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) {
  PetscFunctionBegin;
  Range ents;
  ierr = get_entities_by_ref_level(bit,mask,ents); CHKERRQ(ierr);
  rval = moab.add_entities(meshset,ents); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents) {
  PetscFunctionBegin;
  ierr = moab.get_entities_by_handle(0,ents,false); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();) {
    switch (moab.type_from_handle(*eit)) {
      case MBVERTEX:
      case MBEDGE:
      case MBTRI:
      case MBTET:
      case MBPRISM:
      break;
      case MBENTITYSET:
      default:
	eit = ents.erase(eit);
	continue;
    }
    BitRefLevel bit2;
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
    if(mask.any()&&bit2.none()) {
      eit = ents.erase(eit);
      continue;
    }
    if((bit2&mask) != bit2) {
      eit = ents.erase(eit);
      continue;
    }
    if((bit2&bit).none()) {
      eit = ents.erase(eit);
      continue;
    }
    eit++;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::update_meshset_by_entities_children(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  Range ents;
  rval = moab.get_entities_by_handle(parent,ents,recursive);  
  if(rval != MB_SUCCESS) {
    cerr << parent << endl; 
    cerr << moab.type_from_handle(parent) <<  " " << MBENTITYSET << endl;
  } CHKERR_PETSC(rval);

  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityHandle_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedMoFemEntities.get<Composite_EntityHandle_And_ParentEntityType_mi_tag>();
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
	if(ref_ent == *eit) continue;
	if(ref_ent == 0) {
	  SETERRQ(PETSC_COMM_SELF,1,"this should not happen");
	}
	if(moab.type_from_handle(*eit)==MBENTITYSET) {
	  SETERRQ(PETSC_COMM_SELF,1,"this should not happen");
	}
	rval = moab.add_entities(child,&ref_ent,1); CHKERR_PETSC(rval);
	if(verb>1) {
	  ostringstream ss;
	  ss << "good bit " << *miit << endl;
	  PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::iterator fit = moabFields.begin();
  for(;fit!=moabFields.end();fit++) {
    EntityHandle meshset = fit->get_meshset();
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTET,false,verb);  CHKERRQ(ierr);
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTRI,false,verb);  CHKERRQ(ierr);
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBEDGE,false,verb);  CHKERRQ(ierr);
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBVERTEX,false,verb);  CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::update_field_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,int verb) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator miit;
  miit = moabFields.get<FieldName_mi_tag>().find(name);
  EntityHandle meshset = miit->get_meshset();
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTET,false,verb);  CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTRI,false,verb);  CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBEDGE,false,verb);  CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBVERTEX,false,verb);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::update_finite_element_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb) {
  PetscFunctionBegin;
  typedef MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type finiteElements_by_name;
  const finiteElements_by_name& set = finiteElements.get<MoFEMFiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator miit = set.find(name);
  if(miit==set.end()) THROW_AT_LINE(("finite element < "+name+" > not found (top tip: check spelling)").c_str());
  EntityHandle meshset = miit->get_meshset();
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,fe_ent_type,false,verb);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::problem_get_FE(const string &problem_name,const string &fe_name,const EntityHandle meshset) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<MoFEMProblem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<MoFEMProblem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem like < %s >",problem_name.c_str());
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
bool FieldCore::check_msId_meshset(const int msId,const Cubit_BC_bitset CubitBCType) {
  moabCubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    return true;
  } 
  return false;
}
PetscErrorCode FieldCore::add_Cubit_msId(const Cubit_BC_bitset CubitBCType,const int msId) {
  PetscFunctionBegin;
  if(check_msId_meshset(msId,CubitBCType)) {
    SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",msId);
  } 
  try {
    CubitMeshSets cmeshset(moab,CubitBCType,msId);
    if((cmeshset.CubitBCType&Cubit_BC_bitset(NodeSet|SideSet|BlockSet)).any()) {
      pair<moabCubitMeshSet_multiIndex::iterator,bool> p = cubit_meshsets.insert(cmeshset);
      if(!p.second) SETERRQ(PETSC_COMM_SELF,1,"meshset not inserted");
    }
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::delete_Cubit_msId(const Cubit_BC_bitset CubitBCType,const int msId) {
  PetscFunctionBegin;
  moabCubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit==cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",msId);
  }
  EntityHandle meshset = miit->get_meshset();
  cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().erase(miit);
  rval = moab.delete_entities(&meshset,1); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  moabCubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    ierr = miit->get_Cubit_msId_entities_by_dimension(moab,dimension,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_Cubit_msId_entities_by_dimension(const int msId,const Cubit_BC_bitset CubitBCType,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  moabCubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    ierr = miit->get_Cubit_msId_entities_by_dimension(moab,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
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
PetscErrorCode FieldCore::get_Cubit_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset) {
  PetscFunctionBegin;
  moabCubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    meshset = miit->meshset;
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_Cubit_meshsets(const unsigned int CubitBCType,Range &meshsets) {
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
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",name.c_str());
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
  int count = distance(miit,hi_miit);
  if(count != nb_ghost_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
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
  if(p_x==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",x_problem.c_str());
  moFEMProblems_by_name::iterator p_y = moFEMProblems_set.find(y_problem);
  if(p_y==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",y_problem.c_str());
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
  if(size!=nb_local_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency: check ghost vector, problem with nb. of local nodes");
  ierr = VecGetLocalSize(Vlocal,&size); CHKERRQ(ierr);
  if(size!=nb_local_dofs+nb_ghost_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency: check ghost vector, problem with nb. of ghost nodes");
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
      if(size!=nb_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency: nb. of dofs and declared nb. dofs in database");
      if(size!=distance(miit,hi_miit)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency: nb. of dofs and declared nb. dofs in database");
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
      if(size!=nb_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency: nb. of dofs and declared nb. dofs in database");
      PetscScalar *array;
      VecGetArray(V_glob,&array);
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++) {
	    if(miit->get_petsc_gloabl_dof_idx()>=size) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency: nb. of dofs and declared nb. dofs in database");
	    }
	    DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator diiiit;
	    diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	    if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) {
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
		ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
		SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	      }
	      if(p_e_miit.first->get_max_order()<order) {
		bool success = entsMoabField.modify(p_e_miit.first,MoFEMEntity_change_order(moab,order));
		if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
	      }
	      //create field moabdof
	      DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator hi_diit,diit;
	      diit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_name,miit->get_ent()));
	      hi_diit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_name,miit->get_ent()));
	      for(;diit!=hi_diit;diit++) {
		DofMoFEMEntity mdof(&*(p_e_miit.first),diit->get_dof_order(),diit->get_dof_rank(),diit->get_EntDofIdx());
		pair<DofMoFEMEntity_multiIndex::iterator,bool> cpy_p_diit;
		cpy_p_diit = dofsMoabField.insert(mdof);
		if(cpy_p_diit.second) {
		  bool success = dofsMoabField.modify(cpy_p_diit.first,DofMoFEMEntity_active_change(true));
		  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
		}
	      }
	      diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	      if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
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
	  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator diiiit;
	  diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	  if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) {
	    SETERRQ(PETSC_COMM_SELF,1,"no data to fill the vector (top tip: you want scatter forward or scatter reverse?)");
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
PetscErrorCode FieldCore::field_axpy(const double alpha,const string& field_name_x,const string& field_name_y,
  bool error_if_missing,bool creat_if_missing) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator x_fit = moabFields.get<FieldName_mi_tag>().find(field_name_x);
  if(x_fit==moabFields.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"x field < %s > not found, (top tip: check spelling)",field_name_x.c_str());
  }
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator y_fit = moabFields.get<FieldName_mi_tag>().find(field_name_y);
  if(y_fit==moabFields.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"y field < %s > not found, (top tip: check spelling)",field_name_y.c_str());
  }
  if(x_fit->get_space() != y_fit->get_space()) {
    SETERRQ2(PETSC_COMM_SELF,1,"space for field < %s > and field <%s> are not compatible",field_name_x.c_str(),field_name_y.c_str());
  }
  if(x_fit->get_max_rank() != y_fit->get_max_rank()) {
    SETERRQ2(PETSC_COMM_SELF,1,"rank for field < %s > and field <%s> are not compatible",field_name_x.c_str(),field_name_y.c_str());
  }
  MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator x_eit;
  x_eit = entsMoabField.get<FieldName_mi_tag>().lower_bound(field_name_x.c_str());
  for(;x_eit!=entsMoabField.get<FieldName_mi_tag>().upper_bound(field_name_x.c_str());x_eit++) {
    int nb_dofs_on_x_entity = x_eit->tag_FieldData_size/sizeof(FieldData);
    for(int dd = 0;dd<nb_dofs_on_x_entity;dd++) {
      ApproximationOrder dof_order = x_eit->tag_dof_order_data[dd];
      ApproximationRank dof_rank = x_eit->tag_dof_rank_data[dd];
      FieldData data = x_eit->tag_FieldData[dd];
      DofMoFEMEntity_multiIndex::index<Composite_Name_Ent_Order_And_Rank_mi_tag>::type::iterator dit;
      dit = dofsMoabField.get<Composite_Name_Ent_Order_And_Rank_mi_tag>().find(
	boost::make_tuple(field_name_y.c_str(),x_eit->get_ent(),dof_order,dof_rank));
      if(dit == dofsMoabField.get<Composite_Name_Ent_Order_And_Rank_mi_tag>().end()) {
	if(creat_if_missing) {
	  SETERRQ(PETSC_COMM_SELF,1,"not yet implemented");
	} else {
	  if(error_if_missing) {
	    ostringstream ss;
	    ss << "dof on ent " << x_eit->get_ent() << " order " << dof_order << " rank " << dof_rank << " does not exist";
	    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	  } else {
	    continue;
	  }
	}
      }
      dit->get_FieldData() += alpha*data;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::set_field(const double val,const EntityType type,const string& field_name) {
  PetscFunctionBegin;
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag >::type::iterator dit,hi_dit;
  dit = dofsMoabField.get<Composite_Name_And_Type_mi_tag >().lower_bound(boost::make_tuple(field_name,type));
  hi_dit = dofsMoabField.get<Composite_Name_And_Type_mi_tag >().upper_bound(boost::make_tuple(field_name,type));
  for(;dit!=hi_dit;dit++) {
    dit->get_FieldData() = val;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::field_scale(const double alpha,const string& field_name) {
  PetscFunctionBegin;
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
  for(;dit!=hi_dit;dit++) {
    dit->get_FieldData() *= alpha;
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
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());
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
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());
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
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());
  // finite element
  typedef NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type FEs_by_composite;
  method.fe_name = fe_name;
  method.refinedMoFemEntities = &refinedMoFemEntities;
  method.refinedMoFemElements = &refinedMoFemElements;
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
      ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
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
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem not in database %s",problem_name.c_str());
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
    ierr = method.set_dof(miit->get_DofMoFEMEntity_ptr()); CHKERRQ(ierr);
    ierr = method.set_numered_dof(&*miit); CHKERRQ(ierr);
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::loop_dofs(const string &field_name,EntMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = method.set_fields(&moabFields); CHKERRQ(ierr);
  ierr = method.set_ents_multiIndex(&entsMoabField); CHKERRQ(ierr);
  ierr = method.set_dofs_multiIndex(&dofsMoabField); CHKERRQ(ierr);
  ierr = method.set_fes_multiIndex(&finiteElements); CHKERRQ(ierr);
  ierr = method.set_fes_data_multiIndex(&finiteElementsMoFEMEnts); CHKERRQ(ierr);
  ierr = method.set_adjacencies(&entFEAdjacencies); CHKERRQ(ierr);
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator miit,hi_miit;
  miit = dofsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_miit = dofsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(;miit!=hi_miit;miit++) {
    ierr = method.set_dof(&*miit); CHKERRQ(ierr);
    ierr = method.set_numered_dof(NULL); CHKERRQ(ierr);
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_ref_ents(const RefMoFEMEntity_multiIndex **refinedMoFemEntities_ptr) {
  PetscFunctionBegin;
  *refinedMoFemEntities_ptr = &refinedMoFemEntities;
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
DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator FieldCore::get_dofs_by_name_begin(const string &field_name) const {
  return dofsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
}
DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator FieldCore::get_dofs_by_name_end(const string &field_name) const {
  return dofsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator FieldCore::get_dofs_by_name_and_ent_begin(const string &field_name,const EntityHandle ent) {
  return dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_name,ent));
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator FieldCore::get_dofs_by_name_and_ent_end(const string &field_name,const EntityHandle ent) {
  return dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_name,ent));
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator FieldCore::get_dofs_by_name_and_type_begin(const string &field_name,const EntityType type) {
  return dofsMoabField.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,type));
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator FieldCore::get_dofs_by_name_and_type_end(const string &field_name,const EntityType type) {
  return dofsMoabField.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,type));
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
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator it = moabFields.get<FieldName_mi_tag>().find(name);
  if(it == moabFields.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"field not found < %s >",name.c_str());
  }
  EntityHandle meshset = it->get_meshset();
  int num_entities;
  rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
  if(entsMoabField.get<FieldName_mi_tag>().count(it->get_name()) 
    != (unsigned int)num_entities) {
    SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and field multiindex < %s >",name.c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::check_number_of_ents_in_ents_field() {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator it = moabFields.get<FieldName_mi_tag>().begin();
  for(;it!=moabFields.get<FieldName_mi_tag>().end();it++) {
    if(it->get_space() == NoField) continue; //FIXME: should be treated proprly, not test is just skiped for this NoField space
    EntityHandle meshset = it->get_meshset();
    int num_entities;
    rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
    if(entsMoabField.get<FieldName_mi_tag>().count(it->get_name()) != (unsigned int)num_entities) {
      SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and field multiindex < %s >",it->get_name().c_str());
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::check_number_of_ents_in_ents_finite_element(const string& name) {
  PetscFunctionBegin;
  MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator it;
  it = finiteElements.get<MoFEMFiniteElement_name_mi_tag>().find(name);
  if(it == finiteElements.get<MoFEMFiniteElement_name_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"finite element not found < %s >",name.c_str());
  }
  EntityHandle meshset = it->get_meshset();
  int num_entities;
  rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
  if(finiteElementsMoFEMEnts.get<MoFEMFiniteElement_name_mi_tag>().count(it->get_name().c_str()) 
    != (unsigned int)num_entities) {
    SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and finite elements multiindex < %s >",it->get_name().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::check_number_of_ents_in_ents_finite_element() {
  PetscFunctionBegin;
  MoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator it;
  it = finiteElements.get<MoFEMFiniteElement_name_mi_tag>().begin();
  for(;it!=finiteElements.get<MoFEMFiniteElement_name_mi_tag>().end();it++) {
    EntityHandle meshset = it->get_meshset();
    int num_entities;
    rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
    if(finiteElementsMoFEMEnts.get<MoFEMFiniteElement_name_mi_tag>().count(it->get_name().c_str()) 
      != (unsigned int)num_entities) {
      SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and finite elements multiindex < %s >",it->get_name().c_str());
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) {
  PetscFunctionBegin;
  RefMoFEMEntity from_ref_entiti(moab,from_entiti);
  //cerr << "from:\n";
  //cerr << from_ref_entiti << endl;
  rval = moab.get_adjacencies(&from_entiti,1,to_dimension,false,adj_entities); CHKERR_PETSC(rval);
  Range::iterator eit = adj_entities.begin();
  //cerr << "to:\n";
  for(;eit!=adj_entities.end();) {
    RefMoFEMEntity adj_entiti(moab,*eit);
    //cerr << "\t" << adj_entiti << endl;
    if(from_ref_entiti.get_BitRefLevel() != adj_entiti.get_BitRefLevel()) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) {
  PetscFunctionBegin;
  RefMoFEMEntity from_ref_entiti(moab,from_entiti);
  //cerr << "from:\n";
  //cerr << from_ref_entiti << endl;
  rval = moab.get_adjacencies(&from_entiti,1,to_dimension,false,adj_entities); CHKERR_PETSC(rval);
  Range::iterator eit = adj_entities.begin();
  //cerr << "to:\n";
  for(;eit!=adj_entities.end();) {
    RefMoFEMEntity adj_entiti(moab,*eit);
    //cerr << "\t" << adj_entiti << endl;
    if(!(from_ref_entiti.get_BitRefLevel()&adj_entiti.get_BitRefLevel()).any()) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_adjacencies(
    const MoFEMProblem *problem_ptr,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type,
    const int verb) {
  PetscFunctionBegin;
  BitRefLevel bit = problem_ptr->get_BitRefLevel();
  ierr = get_adjacencies(bit,from_entities,num_netities,to_dimension,adj_entities,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type,const int verb) {
  PetscFunctionBegin;
  if(verb>0) {
    ostringstream ss;
    ss << "from: " << bit << endl << "to: " << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  rval = moab.get_adjacencies(from_entities,num_netities,to_dimension,false,adj_entities,operation_type); CHKERR_PETSC(rval);
  Range::iterator eit = adj_entities.begin();
  //cerr << "to:\n";
  for(;eit!=adj_entities.end();) {
    RefMoFEMEntity adj_entiti(moab,*eit);
    if(verb>0) {
      ostringstream ss;
      ss << "\t" << adj_entiti << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    if(!(adj_entiti.get_BitRefLevel()&bit).any() ) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  PetscFunctionReturn(0);
}

//clear,remove and delete

PetscErrorCode FieldCore::clear_dofs_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  {
    DofMoFEMEntity_multiIndex::iterator dit;
    dit = dofsMoabField.begin();
    for(;dit!=dofsMoabField.end();) {
      BitRefLevel bit2 = dit->get_BitRefLevel(); 
      if(dit->get_ent_type()==MBENTITYSET) {
	dit++;
	continue;
      }
      if((bit2&mask)!=bit2) {
	dit++;
	continue;
      }
      if((bit2&bit).none()) {
	dit++;
	continue;
      }
      dit = dofsMoabField.erase(dit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::clear_dofs_fields(const string &name,const Range ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
    dit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
    hi_dit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
    for(;dit!=hi_dit;) {
      dit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().erase(dit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::clear_ents_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_dofs_fields(bit,mask,verb); CHKERRQ(ierr);
  ierr = clear_adjacencies_entities(bit,mask,verb); CHKERRQ(ierr);
  MoFEMEntity_multiIndex::iterator eit;
  eit = entsMoabField.begin();
  for(;eit!=entsMoabField.end();) {
    if(eit->get_ent_type()==MBENTITYSET) {
      eit++;
      continue;
    }
    BitRefLevel bit2 = eit->get_BitRefLevel(); 
    if((bit2&mask)!=bit2) {
      eit++;
      continue;
    }
    if((bit2&bit).none()) {
      eit++;
      continue;
    }
    EntityHandle ent = eit->get_ent();
    rval = moab.tag_delete_data(eit->field_ptr->th_AppOrder,&ent,1); CHKERR_PETSC(rval);
    if(eit->tag_FieldData_size>0) {
      rval = moab.tag_delete_data(eit->field_ptr->th_FieldData,&ent,1); CHKERR_PETSC(rval);
    }
    eit = entsMoabField.erase(eit);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::clear_ents_fields(const string &name,const Range ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_dofs_fields(name,ents,verb); CHKERRQ(ierr);
  ierr = clear_adjacencies_entities(name,ents,verb); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
    dit = entsMoabField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
    hi_dit = entsMoabField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
    for(;dit!=hi_dit;) {
      dit = entsMoabField.get<Composite_Name_And_Ent_mi_tag>().erase(dit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::clear_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_adjacencies_finite_elements(bit,mask,verb); CHKERRQ(ierr);
  EntMoFEMFiniteElement_multiIndex::iterator fe_it = finiteElementsMoFEMEnts.begin();
  for(;fe_it!=finiteElementsMoFEMEnts.end();) {
    BitRefLevel bit2 = fe_it->get_BitRefLevel(); 
    if(fe_it->get_ent_type()==MBENTITYSET) {
      fe_it++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      fe_it++;
      continue;
    }
    if((bit2&bit).none()) {
      fe_it++;
      continue;
    }
    fe_it = finiteElementsMoFEMEnts.erase(fe_it);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::clear_finite_elements(const string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_adjacencies_finite_elements(name,ents,verb); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    EntMoFEMFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator fit,hi_fit;
    fit = finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
    hi_fit = finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
    for(;fit!=hi_fit;) {
      fit = finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().erase(fit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::clear_adjacencies_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator ait;
  ait = entFEAdjacencies.begin();
  for(;ait!=entFEAdjacencies.end();) {
    BitRefLevel bit2 = ait->EntMoFEMFiniteElement_ptr->get_BitRefLevel(); 
    if(ait->EntMoFEMFiniteElement_ptr->get_ent_type()==MBENTITYSET) {
      ait++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      ait++;
      continue;
    }
    if((bit2&bit).none()) {
      ait++;
      continue;
    }
    ait = entFEAdjacencies.erase(ait);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::clear_adjacencies_finite_elements(const string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<MoABFEEnt_mi_tag>::type::iterator ait,hi_ait;
    ait = entFEAdjacencies.get<MoABFEEnt_mi_tag>().lower_bound(*eit);
    hi_ait = entFEAdjacencies.get<MoABFEEnt_mi_tag>().upper_bound(*eit);
    for(;ait!=hi_ait;) {
      if(ait->EntMoFEMFiniteElement_ptr->get_name() == name) {
	ait = entFEAdjacencies.get<MoABFEEnt_mi_tag>().erase(ait);
      } else {
	ait++;
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::clear_adjacencies_entities(const BitRefLevel &bit,const BitRefLevel &mask,int verb ) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator ait;
  ait = entFEAdjacencies.begin();
  for(;ait!=entFEAdjacencies.end();) {
    BitRefLevel bit2 = ait->MoFEMEntity_ptr->get_BitRefLevel(); 
    if(ait->MoFEMEntity_ptr->get_ent_type()==MBENTITYSET) {
      ait++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      ait++;
      continue;
    }
    if((bit2&bit).none()) {
      ait++;
      continue;
    }
    ait = entFEAdjacencies.erase(ait);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::clear_adjacencies_entities(const string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<MoABEnt_mi_tag>::type::iterator ait,hi_ait;
    ait = entFEAdjacencies.get<MoABEnt_mi_tag>().lower_bound(*eit);
    hi_ait = entFEAdjacencies.get<MoABEnt_mi_tag>().upper_bound(*eit);
    for(;ait!=hi_ait;) {
      if(ait->MoFEMEntity_ptr->get_name() == name) {
	ait = entFEAdjacencies.get<MoABEnt_mi_tag>().erase(ait);
      } else {
	ait++;
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::remove_ents_from_field_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_ents_fields(bit,mask,verb); CHKERRQ(ierr);
  MoFEMField_multiIndex::iterator f_it = moabFields.begin();
  for(;f_it!=moabFields.end();f_it++) {
    EntityHandle meshset = f_it->get_meshset();
    Range ents_to_remove;
    rval = moab.get_entities_by_handle(
      meshset,ents_to_remove,false); CHKERR_PETSC(rval);
    Range::iterator eit = ents_to_remove.begin();
    for(;eit!=ents_to_remove.end();) {
      if(moab.type_from_handle(*eit)==MBENTITYSET) {
	eit = ents_to_remove.erase(eit);
	continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
      if((bit2&mask)!=bit2) {
	eit = ents_to_remove.erase(eit);
	continue;
      }
      if((bit2&bit).none()) {
	eit = ents_to_remove.erase(eit);
	continue;
      }
      MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator iit;
      iit = entsMoabField.get<Composite_Name_And_Ent_mi_tag>().find(boost::make_tuple(f_it->get_name(),*eit));
      if(iit != entsMoabField.get<Composite_Name_And_Ent_mi_tag>().end()) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      eit++;
    }
    rval = moab.remove_entities(meshset,ents_to_remove); CHKERR_PETSC(rval);
    if(verb>0) {
      PetscPrintf(PETSC_COMM_WORLD,"number of removed entities = %u from field %s\n",ents_to_remove.size(),f_it->get_name().c_str());
      if(verb>1) {
	int num_entities;
	rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
	PetscPrintf(PETSC_COMM_WORLD,"\tnumber of entities in database = %u and meshset = %u\n",
	    entsMoabField.get<BitFieldId_mi_tag>().count(f_it->get_id()),num_entities);
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::remove_ents_from_field(const string& name,const EntityHandle meshset,const EntityType type,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents); CHKERR_PETSC(rval);
  ierr = remove_ents_from_field(name,ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::remove_ents_from_field(const string& name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  EntityHandle meshset;
  try {
    meshset = get_field_meshset(name);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  rval = moab.remove_entities(meshset,ents); CHKERR_PETSC(rval);
  ierr = clear_ents_fields(name,ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::remove_ents_from_finite_element_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_finite_elements(bit,mask,verb); CHKERRQ(ierr);
  MoFEMFiniteElement_multiIndex::iterator fe_it = finiteElements.begin();
  for(;fe_it!=finiteElements.end();fe_it++) {
    EntityHandle meshset = fe_it->get_meshset();
    Range ents_to_remove;
    rval = moab.get_entities_by_handle(
      meshset,ents_to_remove,false); CHKERR_PETSC(rval);
    Range::iterator eit = ents_to_remove.begin();
    for(;eit!=ents_to_remove.end();) {
      if(moab.type_from_handle(*eit)==MBENTITYSET) {
	eit = ents_to_remove.erase(eit);
	continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
      if((bit2&mask)!=bit2) {
	eit = ents_to_remove.erase(eit);
	continue;
      }
      if((bit2&bit).none()) {
	eit = ents_to_remove.erase(eit);
	continue;
      }
      EntMoFEMFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator iit;
      iit = finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().find(
	boost::make_tuple(fe_it->get_name(),*eit));
      if(iit != finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().end()) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      eit++;
    }
    rval = moab.remove_entities(meshset,ents_to_remove); CHKERR_PETSC(rval);
    if(verb>0) {
      PetscPrintf(PETSC_COMM_WORLD,"number of removed entities = %u from finite element %s\n",ents_to_remove.size(),fe_it->get_name().c_str());
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::remove_ents_from_finite_element(const string &name,const EntityHandle meshset,const EntityType type,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents,false); CHKERR_PETSC(rval);
  ierr = remove_ents_from_finite_element(name,ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::remove_ents_from_finite_element(const string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_finite_elements(name,ents,verb); CHKERRQ(ierr);
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);
  rval = moab.remove_entities(idm,ents); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = delete_finite_elements_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  ierr = remove_ents_from_field_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  RefMoFEMEntity_multiIndex::iterator ent_it = refinedMoFemEntities.begin();
  for(;ent_it!=refinedMoFemEntities.end();) {
    BitRefLevel bit2 = ent_it->get_BitRefLevel(); 
    if(ent_it->get_ent_type()==MBENTITYSET) {
      ent_it++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      ent_it++;
      continue;
    }
    if((bit2&bit).none()) {
      ent_it++;
      continue;
    }
    ent_it = refinedMoFemEntities.erase(ent_it);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::delete_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  Range ents_to_delete;
  rval = moab.get_entities_by_handle(0,ents_to_delete,false); CHKERR_PETSC(rval);
  {
    Range::iterator eit = ents_to_delete.begin();
    for(;eit!=ents_to_delete.end();) {
      if(moab.type_from_handle(*eit)==MBENTITYSET) {
	eit = ents_to_delete.erase(eit);
	continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
      if((bit2&mask)!=bit2) {
	eit = ents_to_delete.erase(eit);
	continue;
      }
      if((bit2&bit).none()) {
	eit = ents_to_delete.erase(eit);
	continue;
      }
      eit++;
    }
  }
  { //remove parent
    Range::iterator eit = ents_to_delete.begin();
    for(;eit != ents_to_delete.end();eit++) {
      RefMoFEMEntity_multiIndex::index<MoABEnt_MoABEnt_mi_tag>::type::iterator pit,hi_pit;
      pit = refinedMoFemEntities.get<MoABEnt_MoABEnt_mi_tag>().lower_bound(*eit);
      hi_pit = refinedMoFemEntities.get<MoABEnt_MoABEnt_mi_tag>().upper_bound(*eit);
      for(;pit!=hi_pit;pit++) {
	EntityHandle ent = pit->get_ref_ent();
	if(ents_to_delete.find(ent) != ents_to_delete.end()) {
	  continue;
	}
	ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
	if(pcomm->rank()==0) {
	  EntityHandle out_meshset;
	  rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	  rval = moab.add_entities(out_meshset,&ent,1); CHKERR_PETSC(rval);
	  rval = moab.add_entities(out_meshset,&*eit,1); CHKERR_PETSC(rval);
	  rval = moab.write_file("error.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	  rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
	}
	ostringstream ss;
	ss << "child:\n" << *pit << endl;
	ss << "parent:\n" << RefMoFEMEntity(moab,*eit) << endl;
	SETERRQ1(PETSC_COMM_SELF,1,
	  "entity can not be removed, it is parent for some other entity\n%s",ss.str().c_str());
	/*bool success = refinedMoFemEntities.modify(
	  refinedMoFemEntities.project<0>(pit),RefMoFEMEntity_change_remove_parent(moab));
	if(!success) {
	  SETERRQ(PETSC_COMM_SELF,1,"mofification unsucessfull");
	}*/
      }
    }
  }
  { //remove deleted entities form cubit meshsets
    moabCubitMeshSet_multiIndex::iterator cubit_it;
    cubit_it = cubit_meshsets.begin();
    for(;cubit_it!=cubit_meshsets.end();cubit_it++) {
      EntityHandle cubit_meshset = cubit_it->meshset; 
      rval = moab.remove_entities(cubit_meshset,ents_to_delete); CHKERR_PETSC(rval);
    }
  }
  ierr = remove_ents_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"number of deleted entities = %u\n",ents_to_delete.size());

  }
  if(verb>2) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    rval = moab.add_entities(out_meshset,ents_to_delete.subset_by_type(MBTET)); CHKERR_PETSC(rval);
    rval = moab.write_file("debug_ents_to_delete.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  //delete entities form moab
  rval = moab.delete_entities(ents_to_delete); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::delete_finite_elements_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = remove_ents_from_finite_element_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  RefMoFEMElement_multiIndex::iterator fe_it = refinedMoFemElements.begin();
  for(;fe_it!=refinedMoFemElements.end();) {
    BitRefLevel bit2 = fe_it->get_BitRefLevel(); 
    if(fe_it->get_ent_type()==MBENTITYSET) {
      fe_it++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      fe_it++;
      continue;
    }
    if((bit2&bit).none()) {
      fe_it++;
      continue;
    }
    fe_it = refinedMoFemElements.erase(fe_it);
  }
  PetscFunctionReturn(0);
}

}
