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

#include <petscsys.h>
#include <cblas.h>

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>

#include <Common.hpp>
#include <CoreDataStructures.hpp>

namespace MoFEM {

//ref moab MoFEMFiniteElement
RefMoFEMElement::RefMoFEMElement(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr):
  interface_RefMoFEMEntity<RefMoFEMEntity>(_RefMoFEMEntity_ptr) {} 
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

    if(moab.type_from_handle(ent)==MBVERTEX) {
      THROW_AT_LINE("Huston we have problem, vertex (specified by ent) is not part of prism, that is impossible (top tip: check your prisms)");
    } 

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
      PetscTraceBackErrorHandler(
	PETSC_COMM_WORLD,
	__LINE__,PETSC_FUNCTION_NAME,__FILE__,
	MOFEM_DATA_INSONSISTENCY,PETSC_ERROR_INITIAL,
	"this work only for TETs",PETSC_NULL);
      THROW_AT_LINE("this work only for TETs");
  }
  rval = moab.tag_get_handle("_RefType",th_RefType); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefType,&ref_ptr->ent,1,(const void **)&tag_type_data); CHKERR_THROW(rval);
  const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ref_ptr->ent,0,0,0));
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
  pair<SideNumber_multiIndex::iterator,bool> p_miit;
  p_miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,side_number,sense,offset));
  miit = p_miit.first;
  if(miit->ent != ent) {
    PetscTraceBackErrorHandler(
      PETSC_COMM_WORLD,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INSONSISTENCY,PETSC_ERROR_INITIAL,"data inconstency",PETSC_NULL);
    PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INSONSISTENCY,PETSC_ERROR_INITIAL,"data inconstency",PETSC_NULL);
  }
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
  ErrorCode rval;
  int side_number,sense,offset;
  EntityHandle tri = get_ref_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = moab.get_connectivity(tri,conn,num_nodes,true); CHKERR_THROW(rval);
  for(int nn = 0;nn<3; nn++) {
    const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(conn[nn],nn,0,-1));
  }
  for(int ee = 0;ee<3; ee++) {
    EntityHandle edge;
    rval = moab.side_element(tri,1,ee,edge); CHKERR_THROW(rval);
    rval = moab.side_number(tri,edge,side_number,sense,offset); CHKERR_THROW(rval);
    if(ee != side_number) {
      PetscTraceBackErrorHandler(
	PETSC_COMM_WORLD,
	__LINE__,PETSC_FUNCTION_NAME,__FILE__,
	MOFEM_DATA_INSONSISTENCY,PETSC_ERROR_INITIAL,"data inconstency",PETSC_NULL);
      PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,
	__LINE__,PETSC_FUNCTION_NAME,__FILE__,
	MOFEM_DATA_INSONSISTENCY,PETSC_ERROR_INITIAL,"data inconstency",PETSC_NULL);
    }
    const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(edge,ee,sense,offset));
  }
  const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(tri,0,0,0));
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

PetscErrorCode DefaultElementAdjacency::defaultVertex(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency) {
  PetscFunctionBegin;
  switch (field_ptr->get_space()) {
    case H1: 
      adjacency.insert(fe_ptr->get_ent());
      break;
    case NOFIELD: {
      adjacency.insert(field_ptr->get_meshset());
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this field is not implemented for VERTEX finite element");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DefaultElementAdjacency::defaultEdge(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency) {
  PetscFunctionBegin;
  ErrorCode rval;
  EntityHandle fe_ent = fe_ptr->get_ent();
  Range nodes;
  switch (field_ptr->get_space()) {
    case H1: 
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
      adjacency.insert(nodes.begin(),nodes.end());
      adjacency.insert(fe_ent);
      break;
    case NOFIELD: 
      adjacency.insert(field_ptr->get_meshset());
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this field is not implemented for EDGE finite element");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DefaultElementAdjacency::defaultTri(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency) {
  PetscFunctionBegin;
  ErrorCode rval;
  Range nodes,edges;
  EntityHandle fe_ent = fe_ptr->get_ent();
  switch (field_ptr->get_space()) {
    case H1: 
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
      adjacency.insert(nodes.begin(),nodes.end());
      rval = moab.get_adjacencies(&fe_ent,1,1,false,edges); CHKERR_PETSC(rval);
      adjacency.insert(edges.begin(),edges.end());
      for(Range::iterator eeit = edges.begin();eeit!=edges.end();eeit++) {
	fe_ptr->get_side_number_ptr(moab,*eeit);
      }
      //add faces
      adjacency.insert(fe_ent);
      break;
    case HDIV:
      adjacency.insert(fe_ent);
      break;
    case NOFIELD:
      adjacency.insert(field_ptr->get_meshset());
      break;
    case L2:
      //FIXME this is matter of convention what should be done here
      //no ajacencies for L2 field
      //adjacency.insert(fe_ent); // add this just in case, if L2 is on skeleton
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this field is not implemented for TRI finite element");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DefaultElementAdjacency::defaultTet(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency) {
  PetscFunctionBegin;
  ErrorCode rval;
  Range nodes,edges,faces;
  EntityHandle fe_ent = fe_ptr->get_ent();
  switch (field_ptr->get_space()) {
    case H1: 
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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"TET has at least 4 adjacent nodes; it can has more if user add more adjacencies");
      }
      adjacency.insert(nodes.begin(),nodes.end());
      case HCURL: 
	rval = moab.get_adjacencies(&fe_ent,1,1,false,edges); CHKERR_PETSC(rval);
  	adjacency.insert(edges.begin(),edges.end());
	for(Range::iterator eeit = edges.begin();eeit!=edges.end();eeit++) {
	  fe_ptr->get_side_number_ptr(moab,*eeit);
	}
      case HDIV: 
	rval = moab.get_adjacencies(&fe_ent,1,2,false,faces); CHKERR_PETSC(rval);
  	adjacency.insert(faces.begin(),faces.end());
	for(Range::iterator fit = faces.begin();fit!=faces.end();fit++) {
	  fe_ptr->get_side_number_ptr(moab,*fit);
	}
      case L2:
  	adjacency.insert(fe_ent);
      break;
    case NOFIELD:
      adjacency.insert(field_ptr->get_meshset());
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this field is not implemented for TRI finite element");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DefaultElementAdjacency::defaultPrism(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency) {
  PetscFunctionBegin;
  ErrorCode rval;
  EntityHandle fe_ent = fe_ptr->get_ent();
  Range nodes;
  //initialize side sets
  try {
    EntityHandle prism = fe_ent;
    EntityHandle face_side3,face_side4;
    rval = moab.side_element(prism,2,3,face_side3); CHKERR_PETSC(rval);
    rval = moab.side_element(prism,2,4,face_side4); CHKERR_PETSC(rval);
    fe_ptr->get_RefMoFEMElement()->get_side_number_ptr(moab,face_side3);
    fe_ptr->get_RefMoFEMElement()->get_side_number_ptr(moab,face_side4);
    int ee = 0;
    for(;ee<3;ee++) {
      EntityHandle edge;
      rval = moab.side_element(prism,1,ee,edge); CHKERR_PETSC(rval);
      SideNumber *side_ptr = fe_ptr->get_RefMoFEMElement()->get_side_number_ptr(moab,edge);
      if(side_ptr->side_number!=ee) SETERRQ1(PETSC_COMM_SELF,1,"data insonsitency for edge %d",ee);
      rval = moab.side_element(prism,1,6+ee,edge); CHKERR_PETSC(rval);
      side_ptr = fe_ptr->get_RefMoFEMElement()->get_side_number_ptr(moab,edge);
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
      SideNumber *side_ptr = fe_ptr->get_RefMoFEMElement()->get_side_number_ptr(moab,node);
      if(side_ptr->side_number!=nn) SETERRQ1(PETSC_COMM_SELF,1,"data insonsitency for node %d",nn);
      rval = moab.side_element(prism,0,nn+3,node); CHKERR_PETSC(rval);
      side_ptr = fe_ptr->get_RefMoFEMElement()->get_side_number_ptr(moab,node);
      if(side_ptr->side_number!=nn+3) {
	if(side_ptr->side_number!=nn) {
	  SETERRQ1(PETSC_COMM_SELF,1,"data insonsitency for node %d",nn);
	} else {
	  side_ptr->brother_side_number = nn+3; 
	}
      }
    }
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW ,msg);
  }
  //get adjacencies
  SideNumber_multiIndex &side_table = fe_ptr->get_RefMoFEMElement()->get_side_number_table();
  switch(field_ptr->get_space()) {
    case H1: 	
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
      adjacency.insert(nodes.begin(),nodes.end());
    case HCURL: {
      SideNumber_multiIndex::nth_index<2>::type::iterator
      siit = side_table.get<2>().lower_bound(MBEDGE), hi_siit = side_table.get<2>().upper_bound(MBEDGE);
      for(;siit!=hi_siit;siit++) adjacency.insert(siit->ent);
    }
    case HDIV: {
      SideNumber_multiIndex::nth_index<2>::type::iterator
      siit = side_table.get<2>().lower_bound(MBTRI), hi_siit = side_table.get<2>().upper_bound(MBTRI);
      for(;siit!=hi_siit;siit++) adjacency.insert(siit->ent);
    }
    case L2:
      adjacency.insert(fe_ent); 
      break;
    case NOFIELD:
      adjacency.insert(field_ptr->get_meshset());
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this field is not implemented for TRI finite element");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DefaultElementAdjacency::defaultMeshset(Interface &moab,const MoFEMField *field_ptr,const EntMoFEMFiniteElement *fe_ptr,Range &adjacency) {
  PetscFunctionBegin;
  ErrorCode rval;
  Range ent_ents;
  EntityHandle fe_ent = fe_ptr->get_ent();
  //get all meshsets in finite element meshset 
  rval = moab.get_entities_by_type(fe_ent,MBENTITYSET,ent_ents,false); CHKERR_PETSC(rval);
  //resolve recusively all ents in the meshset
  rval = moab.get_entities_by_handle(fe_ent,ent_ents,true); CHKERR_PETSC(rval); 
  Range::iterator eit_eit = ent_ents.begin();
  for(;eit_eit!=ent_ents.end();eit_eit++) {
    switch (field_ptr->get_space()) {
      case NOFIELD:
	if(moab.type_from_handle(*eit_eit)==MBENTITYSET) {
	adjacency.insert(*eit_eit);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
  }
  PetscFunctionReturn(0);
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

  //custom adjacency map
  for(int tt = 0;tt<MBMAXTYPE;tt++) {
    element_adjacency_table[tt] = NULL;
  }

  element_adjacency_table[MBVERTEX] = DefaultElementAdjacency::defaultVertex;
  element_adjacency_table[MBEDGE] = DefaultElementAdjacency::defaultEdge;
  element_adjacency_table[MBTRI] = DefaultElementAdjacency::defaultTri;
  element_adjacency_table[MBTET] = DefaultElementAdjacency::defaultTet;
  element_adjacency_table[MBPRISM] = DefaultElementAdjacency::defaultPrism;
  element_adjacency_table[MBENTITYSET] = DefaultElementAdjacency::defaultMeshset;

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
  //get finite element entity
  global_uid =  get_global_unique_id_calculate();
  //add ents to meshset
  //EntityHandle meshset = get_meshset();
  //EntityHandle ent = get_ent();
  //ierr = moab.add_entities(meshset,&ent,1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

ostream& operator<<(ostream& os,const EntMoFEMFiniteElement& e) {
  os << *e.fe_ptr << endl;
  os << *e.ref_ptr << endl; 
  os << "row dof_uids ";
  DofMoFEMEntity_multiIndex_uid_view::iterator rit;
  rit = e.row_dof_view.begin();
  for(;rit!=e.row_dof_view.end();rit++) {
    os << (*rit)->get_global_unique_id() << " ";
  }
  os << "col dof_uids ";
  DofMoFEMEntity_multiIndex_uid_view::iterator cit;
  cit = e.col_dof_view.begin();
  for(;cit!=e.col_dof_view.end();cit++) {
    os << (*cit)->get_global_unique_id() << " ";
  }
  os << "data dof_uids ";
  DofMoFEMEntity_multiIndex_uid_view::iterator dit;
  dit = e.data_dof_view.begin();
  for(;dit!=e.data_dof_view.end();dit++) {
    os << (*dit)->get_global_unique_id() << " ";
  }
  return os;
}

template <typename MOFEM_DOFS,typename MOFEM_DOFS_VIEW>
static PetscErrorCode get_fe_MoFEMFiniteElement_dof_view(
    const DofMoFEMEntity_multiIndex_uid_view &fe_dofs_view,
    const MOFEM_DOFS &mofem_dofs,
    MOFEM_DOFS_VIEW &mofem_dofs_view,
    const int operation_type) {
  PetscFunctionBegin;
  GlobalUId global_uid;
  typename boost::multi_index::index<MOFEM_DOFS,Unique_mi_tag>::type::iterator mofem_it,mofem_it_end;
  DofMoFEMEntity_multiIndex_uid_view::iterator it,it_end;
  if(operation_type==Interface::UNION) {
    mofem_it = mofem_dofs.template get<Unique_mi_tag>().begin();
    mofem_it_end = mofem_dofs.template get<Unique_mi_tag>().end();
    it = fe_dofs_view.begin();
    it_end = fe_dofs_view.end();
    for(;it!=it_end;it++) {
      global_uid = (*it)->get_global_unique_id();
      if(mofem_it != mofem_it_end) {
	if(mofem_it->get_global_unique_id() != global_uid) {
	  mofem_it = mofem_dofs.template get<Unique_mi_tag>().find(global_uid);
	}
      } else {
	mofem_it = mofem_dofs.template get<Unique_mi_tag>().find(global_uid);
      }
      if(mofem_it != mofem_it_end) {
	mofem_dofs_view.insert(&*mofem_it);
	mofem_it++;
      }
    }
  } else {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_row_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_active_view &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_MoFEMFiniteElement_dof_view(row_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_col_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_active_view &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_MoFEMFiniteElement_dof_view(col_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_data_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_active_view &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_MoFEMFiniteElement_dof_view(data_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_row_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_MoFEMFiniteElement_dof_view(row_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_col_dof_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_MoFEMFiniteElement_dof_view(col_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_row_dof_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view_ordered &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_MoFEMFiniteElement_dof_view(row_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_col_dof_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view_ordered &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_MoFEMFiniteElement_dof_view(col_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_row_dof_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view_hashed &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_MoFEMFiniteElement_dof_view(row_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EntMoFEMFiniteElement::get_MoFEMFiniteElement_col_dof_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view_hashed &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_MoFEMFiniteElement_dof_view(col_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NumeredMoFEMFiniteElement::get_row_dofs_by_petsc_gloabl_dof_idx(DofIdx idx,const FENumeredDofMoFEMEntity **dof_ptr) const {
  PetscFunctionBegin;
  FENumeredDofMoFEMEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = rows_dofs.get<PetscGlobalIdx_mi_tag>().find(idx);
  if(dit == rows_dofs.get<PetscGlobalIdx_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"dof which index < %d > not found",idx);
  }
  *dof_ptr = &*dit;
  PetscFunctionReturn(0);
}

PetscErrorCode NumeredMoFEMFiniteElement::get_col_dofs_by_petsc_gloabl_dof_idx(DofIdx idx,const FENumeredDofMoFEMEntity **dof_ptr) const {
  PetscFunctionBegin;
  FENumeredDofMoFEMEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = rows_dofs.get<PetscGlobalIdx_mi_tag>().find(idx);
  if(dit == rows_dofs.get<PetscGlobalIdx_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"dof which index < %d > not found",idx);
  }
  *dof_ptr = &*dit;
  PetscFunctionReturn(0);
}


}
