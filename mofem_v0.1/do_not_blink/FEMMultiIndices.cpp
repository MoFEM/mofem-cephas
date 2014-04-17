/** \file CoreDataStructures.cpp
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


#include <CoreDataStructures.hpp>

namespace MoFEM {

//ref moab MoFEMFiniteElement
RefMoFEMElement::RefMoFEMElement(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr):
  interface_RefMoFEMEntity<RefMoFEMEntity>(_RefMoFEMEntity_ptr) {
 
  ErrorCode rval;
  const EntityHandle root_meshset = moab.get_root_set();

  Tag th_RefFEMeshset;
  rval = moab.tag_get_handle("_RefFEMeshset",th_RefFEMeshset); CHKERR_THROW(rval);
  EntityHandle meshset;
  rval = moab.tag_get_data(th_RefFEMeshset,&root_meshset,1,&meshset); CHKERR_THROW(rval);
  EntityHandle ent = get_ref_ent();
  rval = moab.add_entities(meshset,&ent,1); CHKERR_THROW(rval);

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
  ErrorCode rval;
  Tag th_RefBitEdge;
  rval = moab.tag_get_handle("_RefBitEdge",th_RefBitEdge); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitEdge,&ref_ptr->ent,1,(const void **)&tag_BitRefEdges); CHKERR_THROW(rval);
  switch (ref_ptr->get_ent_type()) {
    case MBPRISM:
    break;
    default:
      THROW_AT_LINE("this work only for PRISMs");
  }
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
      int sense_m1_map[3][3] = { {0,2,1}, {1,0,2}, {2,1,0} };
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
	  } else {
	    cerr << conn_ent[0] << " " << conn_ent[1] << " " << conn_ent[2] << endl;
	    cerr << face3[0] << " " << face3[1] << " " << face3[2] << endl;
	    cerr << face4[0] << " " << face4[1] << " " << face4[2] << endl;
	    cerr << offset << endl;
	    THROW_AT_LINE("Huston we have problem");
	  }
      } 
      THROW_AT_LINE("Huston we have problem");
    }
    if(num_nodes == 2) {
      EntityHandle edges[6][2] = {
	{ conn[0], conn[1] } /*0*/, { conn[1], conn[2] } /*1*/, { conn[2], conn[0] } /*2*/,
	{ conn[3], conn[4] } /*3+3*/, { conn[4], conn[5] } /*3+4*/, { conn[5], conn[3] } /*3+5*/ };
      for(int ee = 0;ee<6;ee++) {
	if(
	  (( conn_ent[0] == edges[ee][0] )&&( conn_ent[1] == edges[ee][1] ))||
	  (( conn_ent[0] == edges[ee][1] )&&( conn_ent[1] == edges[ee][0] )) ) {
	  side_number = ee;
	  if(ee>=3) {
	    side_number += 3;
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
RefMoFEMElement_TET::RefMoFEMElement_TET(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr): 
  RefMoFEMElement(moab,_RefMoFEMEntity_ptr),tag_BitRefEdges(NULL) {
  ErrorCode rval;
  Tag th_RefBitEdge;
  rval = moab.tag_get_handle("_RefBitEdge",th_RefBitEdge); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitEdge,&ref_ptr->ent,1,(const void **)&tag_BitRefEdges); CHKERR_THROW(rval);
  Tag th_RefType;
  switch (ref_ptr->get_ent_type()) {
    case MBTET:
    break;
    default:
      PetscTraceBackErrorHandler(PETSC_COMM_WORLD,__LINE__,PETSC_FUNCTION_NAME,__FILE__,__SDIR__,1,PETSC_ERROR_INITIAL,
	"this work only for TETs",PETSC_NULL);
      THROW_AT_LINE("this work only for TETs");
      //PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,__LINE__,PETSC_FUNCTION_NAME,__FILE__,__SDIR__,1,PETSC_ERROR_INITIAL,
	//"this work only for TETs",PETSC_NULL);
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
RefMoFEMElement_TRI::RefMoFEMElement_TRI(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr): RefMoFEMElement(moab,_RefMoFEMEntity_ptr) {
  switch (ref_ptr->get_ent_type()) {
    case MBTRI:
    break;
    default:
      THROW_AT_LINE("this work only for TRIs");
  }
}
SideNumber* RefMoFEMElement_TRI::get_side_number_ptr(Interface &moab,EntityHandle ent) const {
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
ostream& operator<<(ostream& os,const RefMoFEMElement_TRI& e) {
  os << *e.ref_ptr;
  return os;
}
RefMoFEMElement_EDGE::RefMoFEMElement_EDGE(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr): RefMoFEMElement(moab,_RefMoFEMEntity_ptr) {
  switch (ref_ptr->get_ent_type()) {
    case MBEDGE:
    break;
    default:
      THROW_AT_LINE("this work only for TRIs");
  }
}
SideNumber* RefMoFEMElement_EDGE::get_side_number_ptr(Interface &moab,EntityHandle ent) const {
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
  if(side_number==-1) THROW_AT_LINE("this is not working");
  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,side_number,sense,offset)).first;
  //cerr << side_number << " " << sense << " " << offset << endl;
  return const_cast<SideNumber*>(&*miit);
}
ostream& operator<<(ostream& os,const RefMoFEMElement_EDGE& e) {
  os << *e.ref_ptr;
  return os;
}
RefMoFEMElement_VERTEX::RefMoFEMElement_VERTEX(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr): RefMoFEMElement(moab,_RefMoFEMEntity_ptr) {
  switch (ref_ptr->get_ent_type()) {
    case MBVERTEX:
    break;
    default:
      THROW_AT_LINE("this works only for TRIs");
  }
}
SideNumber* RefMoFEMElement_VERTEX::get_side_number_ptr(Interface &moab,EntityHandle ent) const {
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
  THROW_AT_LINE("no side entitiy for vertex if its is not an vertex itself");
  return NULL;
}
ostream& operator<<(ostream& os,const RefMoFEMElement_VERTEX& e) {
  os << *e.ref_ptr;
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
    os << "id " << e.get_id() << " name " << e.get_name_ref() << " f_id_row " << e.get_BitFieldId_row() 
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

void MoFEMFiniteElement_col_change_bit_off::operator()(MoFEMFiniteElement &MoFEMFiniteElement) {
  *((BitFieldId*)(MoFEMFiniteElement.tag_BitFieldId_col_data)) &= f_id_col.flip();
}

void MoFEMFiniteElement_row_change_bit_off::operator()(MoFEMFiniteElement &MoFEMFiniteElement) {
  *((BitFieldId*)(MoFEMFiniteElement.tag_BitFieldId_row_data)) &= f_id_row.flip();
}

void EntMoFEMFiniteElement_change_bit_off::operator()(MoFEMFiniteElement &MoFEMFiniteElement) {
  *((BitFieldId*)(MoFEMFiniteElement.tag_BitFieldId_data)) &= f_id_data.flip();
}

//MoFEMFiniteElement data
EntMoFEMFiniteElement::EntMoFEMFiniteElement(Interface &moab,const RefMoFEMElement *_ref_MoFEMFiniteElement,const MoFEMFiniteElement *_MoFEMFiniteElement_ptr): 
  interface_MoFEMFiniteElement<MoFEMFiniteElement>(_MoFEMFiniteElement_ptr),interface_RefMoFEMElement<RefMoFEMElement>(_ref_MoFEMFiniteElement) {
  ErrorCode rval;
  PetscErrorCode ierr;
  //get finite element entity
  EntityHandle ent = get_ent();
  uid = get_unique_id_calculate();
  //set tags
  rval = moab.tag_get_by_ptr(fe_ptr->th_DofUidRow,&ent,1,(const void **)&tag_row_uids_data,&tag_row_uids_size); 
  if(rval != MB_SUCCESS) tag_row_uids_size = 0;
  rval = moab.tag_get_by_ptr(fe_ptr->th_DofUidCol,&ent,1,(const void **)&tag_row_uids_data,&tag_col_uids_size); 
  if(rval != MB_SUCCESS) tag_col_uids_size = 0;
  rval = moab.tag_get_by_ptr(fe_ptr->th_DofUidData,&ent,1,(const void **)&tag_row_uids_data,&tag_data_uids_size);
  if(rval != MB_SUCCESS) tag_data_uids_size = 0;
  //add ents to meshset
  EntityHandle meshset = get_meshset();
  ierr = moab.add_entities(meshset,&ent,1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
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
  const DofMoFEMEntity_multiIndex &dofsMoabField,
  int &side_number, int &sense, int &offset) const { 
  PetscFunctionBegin;
  EntityHandle child = dofsMoabField.get<Unique_mi_tag>().find(_ent_uid_)->get_ent();
  PetscErrorCode ierr;
  ierr = moab.side_number(get_ent(),child,side_number,sense,offset); CHKERRQ(ierr);
  PetscFunctionReturn(0);
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
  int tag_sizes[1]; 
  tag_sizes[0] = data.size()*sizeof(UId);
  void const* tag_data[] = { &data[0] };
  rval = moab.tag_set_by_ptr(th_DofUid,&ent,1,tag_data,tag_sizes); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_DofUid,&ent,1,tag_uids_data,tag_uids_size); CHKERR_THROW(rval);
}

void EntMoFEMFiniteElement_row_dofs_change::operator()(EntMoFEMFiniteElement &MoFEMFiniteElement) { 
  EntMoFEMFiniteElement_dofs_change(
    moab,uids_view,MoFEMFiniteElement.get_ent(),MoFEMFiniteElement.fe_ptr->th_DofUidRow,(const void **)&MoFEMFiniteElement.tag_row_uids_data,&MoFEMFiniteElement.tag_row_uids_size);
}

void EntMoFEMFiniteElement_col_dofs_change::operator()(EntMoFEMFiniteElement &MoFEMFiniteElement) { 
  EntMoFEMFiniteElement_dofs_change(
    moab,uids_view,MoFEMFiniteElement.get_ent(),MoFEMFiniteElement.fe_ptr->th_DofUidCol,(const void **)&MoFEMFiniteElement.tag_col_uids_data,&MoFEMFiniteElement.tag_col_uids_size);
}

void EntMoFEMFiniteElement_data_dofs_change::operator()(EntMoFEMFiniteElement &MoFEMFiniteElement) { 
  EntMoFEMFiniteElement_dofs_change(
    moab,uids_view,MoFEMFiniteElement.get_ent(),MoFEMFiniteElement.fe_ptr->th_DofUidData,(const void **)&MoFEMFiniteElement.tag_data_uids_data,&MoFEMFiniteElement.tag_data_uids_size);
}



}
