/** \file Core_dataStructures.cpp
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


#include <Core_dataStructures.hpp>

namespace MoFEM {

const bool Idx_mi_tag::IamNotPartitioned = true;
const bool Part_mi_tag::IamNotPartitioned = false;

//moab base meshsets
CubitMeshSets::CubitMeshSets(Interface &moab,const EntityHandle _meshset): 
  meshset(_meshset),CubitBCType(UnknownSet),msId(NULL),tag_bc_data(NULL),tag_bc_size(0),tag_block_header_data(NULL) {
  ErrorCode rval;
  Tag nsTag,ssTag,nsTag_data,ssTag_data,bhTag,bhTag_header;
  rval = moab.tag_get_handle(DIRICHLET_SET_TAG_NAME,nsTag); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle(NEUMANN_SET_TAG_NAME,ssTag); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle((string(DIRICHLET_SET_TAG_NAME)+"__BC_DATA").c_str(),nsTag_data); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle((string(NEUMANN_SET_TAG_NAME)+"__BC_DATA").c_str(),ssTag_data); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle(MATERIAL_SET_TAG_NAME,bhTag); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle("BLOCK_HEADER",bhTag_header); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_tags_on_entity(meshset,tag_handles); CHKERR(rval);CHKERR_THROW(rval);
  vector<Tag>::iterator tit = tag_handles.begin();
  for(;tit!=tag_handles.end();tit++) {
    if(
      *tit == nsTag ||
      *tit == ssTag ||
      *tit == bhTag) {
      rval = moab.tag_get_by_ptr(*tit,&meshset,1,(const void **)&msId); CHKERR(rval);CHKERR_THROW(rval);
    }
    if(
      (*tit == nsTag_data)||
      (*tit == ssTag_data)) {
    }
    if(*tit == nsTag) {
      if(*msId != -1) {
	CubitBCType = NodeSet;
      }
    }
    if(*tit == ssTag) {
      if(*msId != -1) {
	CubitBCType = SideSet;
      }
    }
    if(*tit == bhTag) {
      if(*msId != -1) {
	CubitBCType = BlockSet;
      }
    }
    if(
      (*tit == nsTag_data) ||
      (*tit == ssTag_data)) {
      rval = moab.tag_get_by_ptr(*tit,&meshset,1,(const void **)&tag_bc_data,&tag_bc_size); CHKERR(rval);CHKERR_THROW(rval);
    }
    if(*tit == bhTag_header) {
      rval = moab.tag_get_by_ptr(*tit,&meshset,1,(const void **)&tag_block_header_data); CHKERR(rval);CHKERR_THROW(rval);
      if(tag_block_header_data[9]>0) CubitBCType |= MaterialSet;
    }
  }
}
PetscErrorCode CubitMeshSets::get_Cubit_msId_entities_by_dimension(Interface &moab,const int dimension,Range &entities,const bool recursive)  const {
  PetscFunctionBegin;
  ErrorCode rval;
  rval = moab.get_entities_by_dimension(meshset,dimension,entities,recursive); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode CubitMeshSets::get_Cubit_msId_entities_by_dimension(Interface &moab,Range &entities,const bool recursive)  const {
  PetscFunctionBegin;
  if((CubitBCType&Cubit_BC_bitset(BlockSet)).any()) {
    if(tag_block_header_data!=NULL) {
      return get_Cubit_msId_entities_by_dimension(moab,tag_block_header_data[11],entities,recursive);
    } else {
      SETERRQ(PETSC_COMM_SELF,1,"dimension unknown");
    }
  }
  if((CubitBCType&Cubit_BC_bitset(NodeSet)).any()) {
    return get_Cubit_msId_entities_by_dimension(moab,1,entities,recursive);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode CubitMeshSets::get_Cubit_bc_data(vector<char>& bc_data) const {
  PetscFunctionBegin;
  bc_data.resize(tag_bc_size);
  copy(&tag_bc_data[0],&tag_bc_data[tag_bc_size],bc_data.begin());
  PetscFunctionReturn(0);
}
PetscErrorCode CubitMeshSets::print_Cubit_bc_data(ostream& os) const {
  PetscFunctionBegin;
  vector<char> bc_data;
  get_Cubit_bc_data(bc_data);
  os << "bc_data = ";
  std::vector<char>::iterator vit = bc_data.begin();
  for(;vit!=bc_data.end();vit++) {
    os << std::hex << (int)((unsigned char)*vit) << " ";
  }
  os << ": ";
  vit = bc_data.begin();
  for(;vit!=bc_data.end();vit++) {
    os << *vit;
  }
  os << std::endl;
  PetscFunctionReturn(0);
}
ostream& operator<<(ostream& os,const CubitMeshSets& e) {
  os << "meshset " << e.meshset << " type " << e.CubitBCType;
  if(e.msId != NULL) os << " msId " << *(e.msId);
  if(e.tag_block_header_data != NULL) {
    os << " block header: ";
    os << " blockID = " << e.tag_block_header_data[0];
    os << " blockElemType = " << e.tag_block_header_data[1];
    os << " blockMat = " << e.tag_block_header_data[9];
    os << " blockAttributeOrder = " << e.tag_block_header_data[5];
    os << " blockDimension = " << e.tag_block_header_data[11];
  }
  return os;
}

//basic moab ent
BasicMoFEMEntity::BasicMoFEMEntity(const EntityHandle _ent): ent(_ent) {
  switch (get_ent_type()) {
    case MBVERTEX:
    case MBEDGE:
    case MBTRI:
    case MBTET:
    case MBPRISM:
    case MBENTITYSET:
    break;
    default:
      THROW_AT_LINE("this entity type is currently not implemented");
  }
}
//ref moab ent
RefMoFEMEntity::RefMoFEMEntity(
  Interface &moab,const EntityHandle _ent): 
    BasicMoFEMEntity(_ent),tag_parent_ent(NULL),tag_BitRefLevel(NULL) {
  ErrorCode rval;
  Tag th_RefParentHandle,th_RefBitLevel;
  rval = moab.tag_get_handle("_RefParentHandle",th_RefParentHandle); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_RefBitLevel",th_RefBitLevel); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefParentHandle,&ent,1,(const void **)&tag_parent_ent); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel,&ent,1,(const void **)&tag_BitRefLevel); CHKERR_THROW(rval);
}
ostream& operator<<(ostream& os,const RefMoFEMEntity& e) {
  os << "ent " << e.ent;
  os << " parent ent " << e.get_parent_ent();
  os << " BitRefLevel " << e.get_BitRefLevel();
  os << " ent type " << e.get_ent_type();
  os << " ent parent type " << e.get_parent_ent_type();
  return os;
}

//ref moab MoFEMFiniteElement
RefMoFEMElement::RefMoFEMElement(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr):
  interface_RefMoFEMEntity<RefMoFEMEntity>(_RefMoFEMEntity_ptr) {
  ErrorCode rval;
  Tag th_RefBitEdge;
  rval = moab.tag_get_handle("_RefBitEdge",th_RefBitEdge); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitEdge,&ref_ptr->ent,1,(const void **)&tag_BitRefEdges); CHKERR_THROW(rval);
}
ostream& operator<<(ostream& os,const RefMoFEMElement& e) {
  os << " ref egdes " << e.get_BitRefEdges();
  os << " " << *e.ref_ptr;
  return os;
}
RefMoFEMElement_MESHSET::RefMoFEMElement_MESHSET(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr): RefMoFEMElement(moab,_RefMoFEMEntity_ptr) {
  switch (ref_ptr->get_ent_type()) {
    case MBENTITYSET:
    break;
    default:
      THROW_AT_LINE("this work only for MESHSETs");
  }
}
SideNumber* RefMoFEMElement_MESHSET::get_side_number_ptr(Interface &moab,EntityHandle ent) const { 
  NOT_USED(moab);
  NOT_USED(ent);
  SideNumber_multiIndex::iterator miit;
  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,-1,0,-1)).first;
  return const_cast<SideNumber*>(&*miit);
  THROW_AT_LINE("not implemented");
  return NULL;
}
RefMoFEMElement_PRISM::RefMoFEMElement_PRISM(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr): RefMoFEMElement(moab,_RefMoFEMEntity_ptr) {
  switch (ref_ptr->get_ent_type()) {
    case MBPRISM:
    break;
    default:
      THROW_AT_LINE("this work only for PRISMs");
  }
  ErrorCode rval;
  EntityHandle prism = get_ref_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = moab.get_connectivity(prism,conn,num_nodes,true); CHKERR_THROW(rval);
  assert(num_nodes == 6);
  for(int nn = 0;nn<6; nn++) {
    const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(conn[nn],nn,0,-1));
  }
}
SideNumber* RefMoFEMElement_PRISM::get_side_number_ptr(Interface &moab,EntityHandle ent) const {
  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  if(miit!=side_number_table.end()) return const_cast<SideNumber*>(&*miit);
  if(ref_ptr->ent == ent) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,0,0,0)).first;
    return const_cast<SideNumber*>(&*miit);
  }
  ErrorCode rval;
  int side_number,sense,offset;
  rval = moab.side_number(ref_ptr->ent,ent,side_number,sense,offset); CHKERR_THROW(rval);
  if(side_number==-1) {
    //get prism connectivity
    int num_nodes;
    const EntityHandle* conn;
    rval = moab.get_connectivity(ref_ptr->ent,conn,num_nodes,true); CHKERR_THROW(rval);
    assert(num_nodes==6);
    //get ent connectivity	
    const EntityHandle* conn_ent;
    rval = moab.get_connectivity(ent,conn_ent,num_nodes,true); CHKERR_THROW(rval);
    /*
    for(int nn = 0; nn<6;nn++) {
      cerr << conn[nn] << " ";
    };
    cerr << endl;
    for(int nn = 0; nn<num_nodes;nn++) {
      cerr << conn_ent[nn] << " ";
    }
    cerr << endl;
    */
    //buttom face
    EntityHandle face3[3] = { conn[0], conn[1], conn[2] };
    //top face
    EntityHandle face4[3] = { conn[3], conn[4], conn[5] };
    if(num_nodes == 3) {
      int sense_p1_map[3][3] = { {0,1,2}, {1,2,0}, {2,0,1} };
      int sense_m1_map[3][3] = { {0,2,1}, {2,1,0}, {1,0,2} };
      EntityHandle* conn0_3_ptr = find( face3, &face3[3], conn_ent[0] );
      if( conn0_3_ptr != &face3[3] ) {
	offset = distance( face3, conn0_3_ptr );
	if( 
	  face3[ sense_p1_map[offset][0] ] == conn_ent[0] &&
	  face3[ sense_p1_map[offset][1] ] == conn_ent[1] &&
	  face3[ sense_p1_map[offset][2] ] == conn_ent[2] ) {
	  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,3,1,offset)).first;
	  return const_cast<SideNumber*>(&*miit);
	} else if (
	  face3[ sense_m1_map[offset][0] ] == conn_ent[0] &&
	  face3[ sense_m1_map[offset][1] ] == conn_ent[1] &&
	  face3[ sense_m1_map[offset][2] ] == conn_ent[2] ) {
	  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,3,-1,offset)).first;
	  return const_cast<SideNumber*>(&*miit);
	} 
      }	
      EntityHandle* conn0_4_ptr = find( face4, &face4[3], conn_ent[0] );
      if( conn0_4_ptr != &face4[3] ) {
	  offset = distance( face4, conn0_4_ptr );
	  if( 
	    face4[ sense_p1_map[offset][0] ] == conn_ent[0] &&
	    face4[ sense_p1_map[offset][1] ] == conn_ent[1] &&
	    face4[ sense_p1_map[offset][2] ] == conn_ent[2] ) {
	    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,4,1,3+offset)).first;
	    return const_cast<SideNumber*>(&*miit);
	  } else if (
	    face4[ sense_m1_map[offset][0] ] == conn_ent[0] &&
	    face4[ sense_m1_map[offset][1] ] == conn_ent[1] &&
	    face4[ sense_m1_map[offset][2] ] == conn_ent[2] ) {
	    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,4,-1,3+offset)).first;
	    return const_cast<SideNumber*>(&*miit);
	  } else THROW_AT_LINE("Huston we have problem");
      } THROW_AT_LINE("Huston we have problem");
    }
    if(num_nodes == 2) {
      EntityHandle edges[6][2] = {
	{ conn[0], conn[1] }, { conn[1], conn[2] }, { conn[2], conn[0] },
	{ conn[3], conn[4] }, { conn[4], conn[5] }, { conn[5], conn[3] } };
      for(int ee = 0;ee<6;ee++) {
	if(
	  ( conn_ent[0] == edges[ee][0] )&&( conn_ent[1] == edges[ee][1] )||
	  ( conn_ent[0] == edges[ee][1] )&&( conn_ent[1] == edges[ee][0] ) ) {
	  side_number = ee;
	  if(ee>=3) {
	    side_number += 6;
	    EntityHandle* conn0_4_ptr = find( face4, &face4[3], conn_ent[0] );
	    offset = distance( face4, conn0_4_ptr ) + 3;
	  } else {
	    EntityHandle* conn0_3_ptr = find( face3, &face3[3], conn_ent[0] );
	    offset = distance( face3, conn0_3_ptr );
	  }
	  sense = 1;
	  if(( conn_ent[0] == edges[ee][1] )&&( conn_ent[1] == edges[ee][0] ))  sense = -1;
	  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,side_number,sense,offset)).first;
	  return const_cast<SideNumber*>(&*miit);
	}
      }
      THROW_AT_LINE("Huston we have problem");
    }
    ostringstream sss;
    sss << "this not working: " << ent << " type: " << moab.type_from_handle(ent) << " " << MBEDGE << " " << MBTRI << endl;
    THROW_AT_LINE(sss.str().c_str());
  }
  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,side_number,sense,offset)).first;
  return const_cast<SideNumber*>(&*miit);
  THROW_AT_LINE("not implemented");
  return NULL;
}
RefMoFEMElement_TET::RefMoFEMElement_TET(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr): RefMoFEMElement(moab,_RefMoFEMEntity_ptr) {
  ErrorCode rval;
  Tag th_RefType;
  switch (ref_ptr->get_ent_type()) {
    case MBTET:
    break;
    default:
      THROW_AT_LINE("this work only for TETs");
  }
  rval = moab.tag_get_handle("_RefType",th_RefType); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefType,&ref_ptr->ent,1,(const void **)&tag_type_data); CHKERR_THROW(rval);
}
SideNumber* RefMoFEMElement_TET::get_side_number_ptr(Interface &moab,EntityHandle ent) const {
  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  if(miit!=side_number_table.end()) return const_cast<SideNumber*>(&*miit);
  if(ref_ptr->ent == ent) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,0,0,0)).first;
    return const_cast<SideNumber*>(&*miit);
  }
  if(moab.type_from_handle(ent)==MBENTITYSET) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,-1,0,0)).first;
    return const_cast<SideNumber*>(&*miit);
  }
  ErrorCode rval;
  int side_number,sense,offset;
  rval = moab.side_number(ref_ptr->ent,ent,side_number,sense,offset); CHKERR_THROW(rval);
  if(side_number==-1) THROW_AT_LINE("this not working");
  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,side_number,sense,offset)).first;
  //cerr << side_number << " " << sense << " " << offset << endl;
  return const_cast<SideNumber*>(&*miit);
}
ostream& operator<<(ostream& os,const RefMoFEMElement_TET& e) {
  os << "ref type " << e.tag_type_data[0] << " ref sub type " << e.tag_type_data[1];
  os << " ref egdes " << e.get_BitRefEdges();
  os << " " << *e.ref_ptr;
  return os;
}

MoFEMField::MoFEMField(Interface &moab,const EntityHandle _meshset): meshset(_meshset),
  tag_id_data(NULL),tag_space_data(NULL),tag_rank_data(NULL),tag_name_data(NULL),tag_name_size(0) { 
  //Change those tags only by modyfiers
  ErrorCode rval;
  //id
  Tag th_FieldId;
  rval = moab.tag_get_handle("_FieldId",th_FieldId); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FieldId,&meshset,1,(const void **)&tag_id_data); CHKERR(rval);
  //space
  Tag th_FieldSpace;
  rval = moab.tag_get_handle("_FieldSpace",th_FieldSpace); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FieldSpace,&meshset,1,(const void **)&tag_space_data); CHKERR(rval);
  //name
  Tag th_FieldName;
  rval = moab.tag_get_handle("_FieldName",th_FieldName); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FieldName,&meshset,1,(const void **)&tag_name_data,&tag_name_size); CHKERR(rval);
  //data
  string Tag_data_name = "_App_Data_"+get_name();
  rval = moab.tag_get_handle(Tag_data_name.c_str(),th_FieldData); CHKERR(rval);
  //order
  string Tag_ApproximationOrder_name = "_App_Order_"+get_name();
  rval = moab.tag_get_handle(Tag_ApproximationOrder_name.c_str(),th_AppOrder); CHKERR(rval);
  //dof order
  string Tag_dof_ApproximationOrder_name = "_App_Dof_Order"+get_name();
  rval = moab.tag_get_handle(Tag_dof_ApproximationOrder_name.c_str(),th_AppDofOrder); CHKERR(rval);
  //rank
  Tag th_Rank;
  string Tag_rank_name = "_Field_Rank_"+get_name();
  rval = moab.tag_get_handle(Tag_rank_name.c_str(),th_Rank); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_Rank,&meshset,1,(const void **)&tag_rank_data); CHKERR(rval);
  //dof rank
  string Tag_dof_rank_name = "_Field_Dof_Rank_"+get_name();
  rval = moab.tag_get_handle(Tag_dof_rank_name.c_str(),th_DofRank); CHKERR(rval);
  switch (*tag_space_data) {
    case H1:
      forder_entityset = NULL;
      forder_vertex = fNBVERTEX_H1;
      forder_edge = fNBEDGE_H1;
      forder_face = fNBFACE_H1;
      forder_elem = fNBVOLUME_H1;
      break;
    case Hdiv:
      forder_entityset = NULL;
      forder_vertex = fNBVERTEX_Hdiv;
      forder_edge = fNBEDGE_Hdiv;
      forder_face = fNBFACE_Hdiv;
      forder_elem = fNBVOLUME_Hdiv;
      break;
    case Hcurl:
      forder_entityset = NULL;
      forder_vertex = fNBVERTEX_Hcurl;
      forder_edge = fNBEDGE_Hcurl;
      forder_face = fNBFACE_Hcurl;
      forder_elem = fNBVOLUME_Hcurl;
      break;
    case L2:
      forder_entityset = NULL;
      forder_vertex = fNBVERTEX_L2;
      forder_edge = fNBEDGE_L2;
      forder_face = fNBFACE_L2;
      forder_elem = fNBVOLUME_L2;
      break;
    case NoField:
      forder_entityset = fNBENTITYSET_nofield;
      forder_vertex = NULL;
      forder_edge = NULL;
      forder_face = NULL;
      forder_elem = NULL;
      break;
    default:
      THROW_AT_LINE("not implemented");
  }
}
ostream& operator<<(ostream& os,const MoFEMField& e) {
  os << "name "<<e.get_name()<<" BitFieldId "<< e.get_id().to_ulong() << " bit number " << e.get_bit_number() 
    << " space " << e.get_space() << " rank " << e.get_max_rank() << " meshset " << e.meshset;
  return os;
}

//moab problem
MoFEMProblem::MoFEMProblem(Interface &moab,const EntityHandle _meshset): meshset(_meshset) {
  ErrorCode rval;
  Tag th_ProblemId;
  rval = moab.tag_get_handle("_ProblemId",th_ProblemId); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemId,&meshset,1,(const void **)&tag_id_data); CHKERR(rval);
  Tag th_ProblemName;
  rval = moab.tag_get_handle("_ProblemName",th_ProblemName); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemName,&meshset,1,(const void **)&tag_name_data,&tag_name_size); CHKERR(rval);
  Tag th_ProblemNbDofsRow;
  rval = moab.tag_get_handle("_ProblemNbDofsRow",th_ProblemNbDofsRow); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemNbDofsRow,&meshset,1,(const void **)&tag_nbdof_data_row); CHKERR(rval);
  Tag th_ProblemNbDofsCol;
  rval = moab.tag_get_handle("_ProblemNbDofsCol",th_ProblemNbDofsCol); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemNbDofsCol,&meshset,1,(const void **)&tag_nbdof_data_col); CHKERR(rval);
  Tag th_ProblemLocalNbDofRow;
  rval = moab.tag_get_handle("_ProblemLocalNbDofsRow",th_ProblemLocalNbDofRow); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemLocalNbDofRow,&meshset,1,(const void **)&tag_local_nbdof_data_row); CHKERR(rval);
  Tag th_ProblemGhostNbDofRow;
  rval = moab.tag_get_handle("_ProblemGhostNbDofsRow",th_ProblemGhostNbDofRow); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemGhostNbDofRow,&meshset,1,(const void **)&tag_ghost_nbdof_data_row); CHKERR(rval);
  Tag th_ProblemLocalNbDofCol;
  rval = moab.tag_get_handle("_ProblemLocalNbDofsCol",th_ProblemLocalNbDofCol); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemLocalNbDofCol,&meshset,1,(const void **)&tag_local_nbdof_data_col); CHKERR(rval);
  Tag th_ProblemGhostNbDofCol;
  rval = moab.tag_get_handle("_ProblemGhostNbDofsCol",th_ProblemGhostNbDofCol); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemGhostNbDofCol,&meshset,1,(const void **)&tag_ghost_nbdof_data_col); CHKERR(rval);
  Tag th_ProblemFEId;
  rval = moab.tag_get_handle("_ProblemFEId",th_ProblemFEId); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_ProblemFEId,&meshset,1,(const void **)&tag_BitFEId_data); CHKERR(rval);
  Tag th_RefBitLevel;
  rval = moab.tag_get_handle("_RefBitLevel",th_RefBitLevel); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_RefBitLevel,&meshset,1,(const void **)&tag_BitRefLevel); CHKERR(rval);
}
ostream& operator<<(ostream& os,const MoFEMProblem& e) {
  os << "problem id " << e.get_id()
    << " MoFEMFiniteElement id " << e.get_BitFEId() 
    << " name "<<e.get_name();
  return os;
}
BitFEId MoFEMProblem::get_BitFEId() const {
  return *tag_BitFEId_data;
}
void problem_MoFEMFiniteElement_change_bit_add::operator()(MoFEMProblem &p) {
  *(p.tag_BitFEId_data) |= f_id;
}
problem_row_change::problem_row_change(const DofMoFEMEntity *_dof_ptr): dof_ptr(_dof_ptr) {
  assert(dof_ptr->active);
}
void problem_row_change::operator()(_MoFEMProblem_ &e) { 
  pair<NumeredDofMoFEMEntity_multiIndex::iterator,bool> p 
    = e.numered_dofs_rows.insert(NumeredDofMoFEMEntity((*(DofIdx*)e.tag_nbdof_data_row),dof_ptr)); 
  if(p.second) {
    (*(DofIdx*)e.tag_nbdof_data_row)++;
  }
}
problem_col_change::problem_col_change(const DofMoFEMEntity *_dof_ptr): dof_ptr(_dof_ptr) {}
void problem_col_change::operator()(_MoFEMProblem_ &e) { 
  pair<NumeredDofMoFEMEntity_multiIndex::iterator,bool> p 
    = e.numered_dofs_cols.insert(NumeredDofMoFEMEntity((*(DofIdx*)e.tag_nbdof_data_col),dof_ptr)); 
  if(p.second) {
    (*(DofIdx*)e.tag_nbdof_data_col)++;
  }
}
void problem_zero_nb_rows_change::operator()(_MoFEMProblem_ &e) { 
  (*(DofIdx*)e.tag_nbdof_data_row) = 0;
  e.numered_dofs_rows.clear();
}
void problem_zero_nb_cols_change::operator()(_MoFEMProblem_ &e) { 
  (*(DofIdx*)e.tag_nbdof_data_col) = 0;
  e.numered_dofs_cols.clear();
}

//moab ent
MoFEMEntity::MoFEMEntity(Interface &moab,const MoFEMField *_FieldData,const RefMoFEMEntity *_ref_ent_ptr): 
  interface_MoFEMField<MoFEMField>(_FieldData),interface_RefMoFEMEntity<RefMoFEMEntity>(_ref_ent_ptr),ref_mab_ent_ptr(_ref_ent_ptr),
  tag_order_data(NULL),tag_FieldData(NULL),tag_FieldData_size(0),tag_dof_order_data(NULL),tag_dof_rank_data(NULL) {
  const EntityType type = get_ent_type(); 
  switch (type) {
    case MBVERTEX:
      forder = field_ptr->forder_vertex;
      break;
    case MBEDGE:
      forder = field_ptr->forder_edge;
      break;
    case MBTRI:
      forder = field_ptr->forder_face;
      break;
    case MBTET:
      forder = field_ptr->forder_elem;
      break;
    case MBPRISM:
      forder = field_ptr->forder_elem;
      break;
    case MBENTITYSET:
      forder = field_ptr->forder_entityset;
      break;
    default: {
      THROW_AT_LINE("data inconsistency");
    }
  }
  ErrorCode rval;
  EntityHandle ent = get_ent();
  rval = moab.tag_get_by_ptr(field_ptr->th_AppOrder,&ent,1,(const void **)&tag_order_data); CHKERR(rval);
  uid = get_unique_id_calculate();
  rval = moab.tag_get_by_ptr(field_ptr->th_FieldData,&ent,1,(const void **)&tag_FieldData,&tag_FieldData_size); 
  if(rval == MB_SUCCESS) {
    if( (unsigned int)tag_FieldData_size != 0 ) {
      int tag_size[1];
      rval = moab.tag_get_by_ptr(field_ptr->th_AppDofOrder,&ent,1,(const void **)&tag_dof_order_data,tag_size); CHKERR(rval);
      assert(tag_size[0]/sizeof(ApproximationOrder) == tag_FieldData_size/sizeof(FieldData));
      rval = moab.tag_get_by_ptr(field_ptr->th_DofRank,&ent,1,(const void **)&tag_dof_rank_data,tag_size); CHKERR(rval);
      assert(tag_size[0]/sizeof(ApproximationRank) == tag_FieldData_size/sizeof(FieldData));
    }
  }
}
ostream& operator<<(ostream& os,const MoFEMEntity& e) {
  os << "ent_uid " << e.get_unique_id() << " entity "<< e.get_ent() << " type " << e.get_ent_type()
    <<" order "<<e.get_max_order()<<" "<<*e.field_ptr;
  return os;
}
void MoFEMEntity_change_order::operator()(MoFEMEntity &e) {
  ErrorCode rval;
  int nb_dofs = e.forder(order)*e.get_max_rank();
  ApproximationOrder& ent_order = *((ApproximationOrder*)e.tag_order_data);
  ent_order = order;
  EntityHandle ent = e.get_ent();
  rval = moab.tag_get_by_ptr(e.field_ptr->th_FieldData,&ent,1,(const void **)&e.tag_FieldData,&e.tag_FieldData_size); 
  if(rval == MB_SUCCESS) {
    //data
    if( nb_dofs*sizeof(FieldData) == (unsigned int)e.tag_FieldData_size ) {
      rval = moab.tag_get_by_ptr(e.field_ptr->th_FieldData,&ent,1,(const void **)&e.tag_FieldData,&e.tag_FieldData_size); CHKERR(rval);
      int tag_size[1];
      rval = moab.tag_get_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,(const void **)&e.tag_dof_order_data,tag_size); CHKERR(rval);
      assert(tag_size[0]/sizeof(ApproximationOrder) == e.tag_FieldData_size/sizeof(FieldData));
      rval = moab.tag_get_by_ptr(e.field_ptr->th_DofRank,&ent,1,(const void **)&e.tag_dof_rank_data,tag_size); CHKERR(rval);
      assert(tag_size[0]/sizeof(ApproximationRank) == e.tag_FieldData_size/sizeof(FieldData));
      return;
    }
    data.resize(e.tag_FieldData_size/sizeof(FieldData));
    FieldData *ptr_begin = (FieldData*)e.tag_FieldData;
    FieldData *ptr_end = (FieldData*)e.tag_FieldData + e.tag_FieldData_size/sizeof(FieldData);
    copy(ptr_begin,ptr_end,data.begin());
    int tag_size[1];
    //order
    rval = moab.tag_get_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,(const void **)&e.tag_dof_order_data,tag_size); 
    assert(rval == MB_SUCCESS);
    assert(tag_size[0]/sizeof(ApproximationOrder) == e.tag_FieldData_size/sizeof(FieldData));
    data_dof_order.resize(e.tag_FieldData_size/sizeof(FieldData));
    ApproximationOrder *ptr_dof_order_begin = (ApproximationOrder*)e.tag_dof_order_data;
    ApproximationOrder *ptr_dof_order_end = (ApproximationOrder*)e.tag_dof_order_data + e.tag_FieldData_size/sizeof(FieldData);
    copy(ptr_dof_order_begin,ptr_dof_order_end,data_dof_order.begin());
    //rank
    rval = moab.tag_get_by_ptr(e.field_ptr->th_DofRank,&ent,1,(const void **)&e.tag_dof_rank_data,tag_size); 
    assert(rval == MB_SUCCESS);
    assert(tag_size[0]/sizeof(ApproximationRank) == e.tag_FieldData_size/sizeof(FieldData));
    data_dof_rank.resize(e.tag_FieldData_size/sizeof(FieldData));
    ApproximationRank *ptr_dof_rank_begin = (ApproximationRank*)e.tag_dof_rank_data;
    ApproximationRank *ptr_dof_rank_end = (ApproximationRank*)e.tag_dof_rank_data + e.tag_FieldData_size/sizeof(FieldData);
    copy(ptr_dof_rank_begin,ptr_dof_rank_end,data_dof_rank.begin());
  }
  if(nb_dofs>0) {
    data.resize(nb_dofs,0);
    int tag_size[] = { data.size()*sizeof(FieldData) };
    void const* tag_data[] = { &data[0] };
    rval = moab.tag_set_by_ptr(e.field_ptr->th_FieldData,&ent,1,tag_data,tag_size); CHKERR(rval); assert(rval==MB_SUCCESS);
    rval = moab.tag_get_by_ptr(e.field_ptr->th_FieldData,&ent,1,(const void **)&e.tag_FieldData,&e.tag_FieldData_size); CHKERR(rval);
    //order
    data_dof_order.resize(nb_dofs,0);
    tag_size[0] = data_dof_order.size()*sizeof(ApproximationOrder);
    tag_data[0] = &data_dof_order[0];
    rval = moab.tag_set_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,tag_data,tag_size); CHKERR(rval); assert(rval==MB_SUCCESS);
    rval = moab.tag_get_by_ptr(e.field_ptr->th_AppDofOrder,&ent,1,(const void **)&e.tag_dof_order_data,tag_size); CHKERR(rval);
    //rank
    data_dof_rank.resize(nb_dofs,0);
    tag_size[0] = data_dof_rank.size()*sizeof(ApproximationRank);
    tag_data[0] = &data_dof_rank[0];
    rval = moab.tag_set_by_ptr(e.field_ptr->th_DofRank,&ent,1,tag_data,tag_size); CHKERR(rval); assert(rval==MB_SUCCESS);
    rval = moab.tag_get_by_ptr(e.field_ptr->th_DofRank,&ent,1,(const void **)&e.tag_dof_rank_data,tag_size); CHKERR(rval);
  }
}

//moab dof
DofMoFEMEntity::DofMoFEMEntity(const MoFEMEntity *_MoFEMEntity_ptr,const ApproximationOrder _dof_order,const ApproximationRank _dof_rank,const DofIdx _dof): 
    interface_MoFEMEntity<MoFEMEntity>(_MoFEMEntity_ptr), dof(_dof),active(false) {
  assert(field_ptr->tag_dof_order_data!=NULL);
  assert(field_ptr->tag_dof_rank_data!=NULL);
  ((ApproximationOrder*)field_ptr->tag_dof_order_data)[dof] = _dof_order;
  ((ApproximationRank*)field_ptr->tag_dof_rank_data)[dof] = _dof_rank;
  uid = get_unique_id_calculate();
}
ostream& operator<<(ostream& os,const DofMoFEMEntity& e) {
  os << "dof_uid " << e.get_unique_id() 
    << " dof_order " << e.get_dof_order()
    << " dof_rank " << e.get_dof_rank()
    << " dof " << e.dof
    << " active " << e.active 
    << " " << *(e.field_ptr);
  return os;
}
DofMoFEMEntity_active_change::DofMoFEMEntity_active_change(bool _active): active(_active) {}
void DofMoFEMEntity_active_change::operator()(DofMoFEMEntity &_dof_) {
  _dof_.active = active;
  assert((_dof_.get_dof_order()<=_dof_.get_max_order()));
}


//numered dof
NumeredDofMoFEMEntity::NumeredDofMoFEMEntity(const DofIdx idx,const DofMoFEMEntity* _DofMoFEMEntity_ptr): 
    interface_DofMoFEMEntity<DofMoFEMEntity>(_DofMoFEMEntity_ptr),
    dof_idx(idx),petsc_gloabl_dof_idx(-1),petsc_local_dof_idx(-1),part(-1) {}
ostream& operator<<(ostream& os,const NumeredDofMoFEMEntity& e) {
  os << "idx " << e.dof_idx << " part " << e.part 
    << " petsc idx " << e.petsc_gloabl_dof_idx 
    << " ( " << e.petsc_local_dof_idx <<  " ) "
    << *e.field_ptr;
  return os;
}

FEDofMoFEMEntity::FEDofMoFEMEntity(
  SideNumber *_side_number_ptr,
  const DofMoFEMEntity *_DofMoFEMEntity_ptr): 
  BaseFEDofMoFEMEntity(_side_number_ptr), interface_DofMoFEMEntity<DofMoFEMEntity>(_DofMoFEMEntity_ptr) {}
ostream& operator<<(ostream& os,const FEDofMoFEMEntity& e) {
  os << "local dof MoFEMFiniteElement idx " 
    << "side_number " << e.side_number_ptr->side_number << " "
    << "sense " << e.side_number_ptr->sense << " "
    << *e.field_ptr;
  return os;
}

FENumeredDofMoFEMEntity::FENumeredDofMoFEMEntity(
  SideNumber *_side_number_ptr,
  const NumeredDofMoFEMEntity *_NumeredDofMoFEMEntity_ptr): 
  BaseFEDofMoFEMEntity(_side_number_ptr), interface_NumeredDofMoFEMEntity<NumeredDofMoFEMEntity>(_NumeredDofMoFEMEntity_ptr) {}
ostream& operator<<(ostream& os,const FENumeredDofMoFEMEntity& e) {
  os << "local dof MoFEMFiniteElement idx " 
    << "side_number " << e.side_number_ptr->side_number << " "
    << "sense " << e.side_number_ptr->sense << " "
    << *e.field_ptr;
  return os;
}

//MoFEMFiniteElement
MoFEMFiniteElement::MoFEMFiniteElement(Interface &moab,const EntityHandle _meshset): meshset(_meshset) {
  ErrorCode rval;
  Tag th_FEId;
  rval = moab.tag_get_handle("_FEId",th_FEId); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FEId,&meshset,1,(const void **)&tag_id_data); CHKERR(rval);
  Tag th_FEName;
  rval = moab.tag_get_handle("_FEName",th_FEName); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FEName,&meshset,1,(const void **)&tag_name_data,&tag_name_size); CHKERR(rval);
  Tag th_FEIdCol,th_FEIdRow,th_FEIdData;
  rval = moab.tag_get_handle("_FEIdCol",th_FEIdCol); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FEIdCol,&meshset,1,(const void **)&tag_BitFieldId_col_data); CHKERR(rval);
  rval = moab.tag_get_handle("_FEIdRow",th_FEIdRow); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FEIdRow,&meshset,1,(const void **)&tag_BitFieldId_row_data); CHKERR(rval);
  rval = moab.tag_get_handle("_FEIdData",th_FEIdData); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_FEIdData,&meshset,1,(const void **)&tag_BitFieldId_data); CHKERR(rval);
  string Tag_mat_name = "_FE_MatData_"+get_name();
  rval = moab.tag_get_handle(Tag_mat_name.c_str(),th_FEMatData); CHKERR(rval);
  string Tag_vec_name = "_FE_VecData_"+get_name();
  rval = moab.tag_get_handle(Tag_vec_name.c_str(),th_FEVecData); CHKERR(rval);
  string Tag_DofUidRow_name = "_DofUidRow_"+get_name();
  rval = moab.tag_get_handle(Tag_DofUidRow_name.c_str(),th_DofUidRow); CHKERR(rval);
  string Tag_DofUidCol_name = "_DofUidCol_"+get_name();
  rval = moab.tag_get_handle(Tag_DofUidCol_name.c_str(),th_DofUidCol); CHKERR(rval);
  string Tag_DofUidData_name = "_DofUidData_"+get_name();
  rval = moab.tag_get_handle(Tag_DofUidData_name.c_str(),th_DofUidData); CHKERR(rval);
}
ostream& operator<<(ostream& os,const MoFEMFiniteElement& e) {
    os << "id " << e.get_id() << " name " << e.get_name() << " f_id_row " << e.get_BitFieldId_row() 
    << " f_id_col " << e.get_BitFieldId_col() << " BitFEId_data " << e.get_BitFieldId_data();
    return os;
}
void MoFEMFiniteElement_col_change_bit_add::operator()(MoFEMFiniteElement &MoFEMFiniteElement) {
  *((BitFieldId*)(MoFEMFiniteElement.tag_BitFieldId_col_data)) |= f_id_col;
}
void MoFEMFiniteElement_row_change_bit_add::operator()(MoFEMFiniteElement &MoFEMFiniteElement) {
  *((BitFieldId*)(MoFEMFiniteElement.tag_BitFieldId_row_data)) |= f_id_row;
}
void EntMoFEMFiniteElement_change_bit_add::operator()(MoFEMFiniteElement &MoFEMFiniteElement) {
  *((BitFieldId*)(MoFEMFiniteElement.tag_BitFieldId_data)) |= f_id_data;
}
static void EntMoFEMFiniteElement_dofs_change(
  Interface &moab,const DofMoFEMEntity_multiIndex_uid_view &uids_view,const EntityHandle ent,const Tag th_DofUid,
  const void** tag_uids_data,int *tag_uids_size) {
  if(uids_view.empty()) return;
  ErrorCode rval;
  vector<UId> data;
  data.resize(uids_view.size());
  DofMoFEMEntity_multiIndex_uid_view::iterator miit = uids_view.begin();
  vector<UId>::iterator vit = data.begin();
  for(;miit!=uids_view.end();miit++,vit++) *vit = (*miit)->get_unique_id();
  assert(vit==data.end());
  int tag_sizes[] = { data.size()*sizeof(UId) };
  void const* tag_data[] = { &data[0] };
  rval = moab.tag_set_by_ptr(th_DofUid,&ent,1,tag_data,tag_sizes); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_DofUid,&ent,1,tag_uids_data,tag_uids_size); CHKERR_THROW(rval);
}
void EntMoFEMFiniteElement_row_dofs_change::operator()(EntMoFEMFiniteElement &MoFEMFiniteElement) { 
  EntMoFEMFiniteElement_dofs_change(moab,uids_view,MoFEMFiniteElement.get_ent(),MoFEMFiniteElement.fe_ptr->th_DofUidRow,(const void **)&MoFEMFiniteElement.tag_row_uids_data,&MoFEMFiniteElement.tag_row_uids_size);
}
void EntMoFEMFiniteElement_col_dofs_change::operator()(EntMoFEMFiniteElement &MoFEMFiniteElement) { 
  EntMoFEMFiniteElement_dofs_change(moab,uids_view,MoFEMFiniteElement.get_ent(),MoFEMFiniteElement.fe_ptr->th_DofUidCol,(const void **)&MoFEMFiniteElement.tag_col_uids_data,&MoFEMFiniteElement.tag_col_uids_size);
}
void EntMoFEMFiniteElement_data_dofs_change::operator()(EntMoFEMFiniteElement &MoFEMFiniteElement) { 
  EntMoFEMFiniteElement_dofs_change(moab,uids_view,MoFEMFiniteElement.get_ent(),MoFEMFiniteElement.fe_ptr->th_DofUidData,(const void **)&MoFEMFiniteElement.tag_data_uids_data,&MoFEMFiniteElement.tag_data_uids_size);
}

//MoFEMFiniteElement data
EntMoFEMFiniteElement::EntMoFEMFiniteElement(Interface &moab,const RefMoFEMElement *_ref_MoFEMFiniteElement,const MoFEMFiniteElement *_MoFEMFiniteElement_ptr): 
  interface_MoFEMFiniteElement<MoFEMFiniteElement>(_MoFEMFiniteElement_ptr),interface_RefMoFEMElement<RefMoFEMElement>(_ref_MoFEMFiniteElement) {
  ErrorCode rval;
  EntityHandle ent = get_ent();
  rval = moab.tag_get_by_ptr(fe_ptr->th_DofUidRow,&ent,1,(const void **)&tag_row_uids_data,&tag_row_uids_size); 
  if(rval != MB_SUCCESS) tag_row_uids_size = 0;
  rval = moab.tag_get_by_ptr(fe_ptr->th_DofUidCol,&ent,1,(const void **)&tag_row_uids_data,&tag_col_uids_size); 
  if(rval != MB_SUCCESS) tag_col_uids_size = 0;
  rval = moab.tag_get_by_ptr(fe_ptr->th_DofUidData,&ent,1,(const void **)&tag_row_uids_data,&tag_data_uids_size);
  if(rval != MB_SUCCESS) tag_data_uids_size = 0;
}
ostream& operator<<(ostream& os,const EntMoFEMFiniteElement& e) {
  os << *e.fe_ptr << " ent " << e.get_ent() << endl; 
  DofMoFEMEntity_multiIndex_uid_view::iterator miit;
  unsigned int ii = 0;
  if(e.tag_row_uids_size/sizeof(UId)>0) os << "row dof_uids ";
  for(;ii<e.tag_row_uids_size/sizeof(UId);ii++) os << ((UId*)e.tag_row_uids_data)[ii] << " ";
  if(ii!=0) os << endl;
  if(e.tag_col_uids_size/sizeof(UId)>0) os << "col dof_uids ";
  for(ii = 0;ii<e.tag_col_uids_size/sizeof(UId);ii++) os << ((UId*)e.tag_col_uids_data)[ii] << " ";
  if(ii!=0) os << endl;
  if(e.tag_data_uids_size/sizeof(UId)>0) os << "data dof_uids ";
  for(ii = 0;ii<e.tag_data_uids_size/sizeof(UId);ii++) os << ((UId*)e.tag_data_uids_data)[ii] << " ";
  return os;
}
PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_row_dof_uid_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_MoFEMFiniteElement_dof_uid_view(dofs,dofs_view,operation_type,tag_row_uids_data,tag_row_uids_size); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_col_dof_uid_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_MoFEMFiniteElement_dof_uid_view(dofs,dofs_view,operation_type,tag_col_uids_data,tag_col_uids_size); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_row_dof_uid_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_MoFEMFiniteElement_dof_uid_view(dofs,dofs_view,operation_type,tag_row_uids_data,tag_row_uids_size); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_col_dof_uid_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_MoFEMFiniteElement_dof_uid_view(dofs,dofs_view,operation_type,tag_col_uids_data,tag_col_uids_size); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntMoFEMFiniteElement::get_uid_side_number(
  Interface &moab,const UId _ent_uid_,
  const DofMoFEMEntity_multiIndex &dofs_moabfield,
  int &side_number, int &sense, int &offset) const { 
  PetscFunctionBegin;
  EntityHandle child = dofs_moabfield.get<Unique_mi_tag>().find(_ent_uid_)->get_ent();
  PetscErrorCode ierr;
  ierr = moab.side_number(get_ent(),child,side_number,sense,offset); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

//MoFEMAdjacencies
MoFEMAdjacencies::MoFEMAdjacencies(const MoFEMEntity *_MoFEMEntity_ptr,const EntMoFEMFiniteElement *_EntMoFEMFiniteElement_ptr):
  by_other(0),MoFEMEntity_ptr(_MoFEMEntity_ptr),EntMoFEMFiniteElement_ptr(_EntMoFEMFiniteElement_ptr) {}
ostream& operator<<(ostream& os,const MoFEMAdjacencies& e) {
  os << "by_other " << bitset<3>(e.by_other) << " "
    << *e.MoFEMEntity_ptr << endl << *e.EntMoFEMFiniteElement_ptr->fe_ptr;
  return os;
}
PetscErrorCode MoFEMAdjacencies::get_ent_adj_dofs_bridge(
    const DofMoFEMEntity_multiIndex &dofs_moabfield,const by_what _by,
    DofMoFEMEntity_multiIndex_uid_view &uids_view,const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  switch (_by) {
    case by_row:  
      ierr = EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_row_dof_uid_view(dofs_moabfield,uids_view,operation_type); CHKERRQ(ierr);
      break;
    case by_col: 
      ierr = EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_col_dof_uid_view(dofs_moabfield,uids_view,operation_type); CHKERRQ(ierr);
      break;
    default:
      ostringstream ss;
      ss << "don't know that to do for elem " << EntMoFEMFiniteElement_ptr->get_name() << " and field " << MoFEMEntity_ptr->get_name();
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode MoFEMAdjacencies::get_ent_adj_dofs_bridge(
    const NumeredDofMoFEMEntity_multiIndex &dofs_moabproblem,const by_what _by,
    NumeredDofMoFEMEntity_multiIndex_uid_view &uids_view,const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  switch (_by) {
    case by_row:  
      ierr = EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_row_dof_uid_view(dofs_moabproblem,uids_view,operation_type); CHKERRQ(ierr);
      break;
    case by_col: 
      ierr = EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_col_dof_uid_view(dofs_moabproblem,uids_view,operation_type); CHKERRQ(ierr);
      break;
    default: 
      ostringstream ss;
      ss << *this << endl;
      ss << "don't know that to do for elem " << EntMoFEMFiniteElement_ptr->get_name() << " and field " << MoFEMEntity_ptr->get_name();
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

//....
PetscErrorCode test_moab(Interface &moab,const EntityHandle ent) {
  PetscFunctionBegin;
  //tets type
  EntityType type = (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH);
  if(type != moab.type_from_handle(ent)) SETERRQ(PETSC_COMM_SELF,1,"incositencies with type_from_handle"); 
  //tets id
  EntityID id = (EntityType)(ent&MB_ID_MASK);
  if(id != moab.id_from_handle(ent)) SETERRQ(PETSC_COMM_SELF,1,"incositencies with id_from_handle"); 
  PetscFunctionReturn(0);
}

}
