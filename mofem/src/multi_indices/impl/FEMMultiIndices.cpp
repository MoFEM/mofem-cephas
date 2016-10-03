/** \file CoreDataStructures.cpp
 * \brief Mylti-index contains data structures and other low-level functions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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

#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>

namespace MoFEM {

//ref moab FiniteElement
RefElement::RefElement(Interface &moab,const boost::shared_ptr<RefEntity> ref_ent_ptr):
interface_RefEntity<RefEntity>(ref_ent_ptr) {}

std::ostream& operator<<(std::ostream& os,const RefElement& e) {
  os << " ref egdes " << e.getBitRefEdges();
  os << " " << *(e.sPtr);
  return os;
}

RefElement_MESHSET::RefElement_MESHSET(Interface &moab,const boost::shared_ptr<RefEntity> ref_ent_ptr):
RefElement(moab,ref_ent_ptr) {
  switch (ref_ent_ptr->getEntType()) {
    case MBENTITYSET:
    break;
    default:
      THROW_MESSAGE("this work only for MESHSETs");
  }
}
boost::shared_ptr<SideNumber> RefElement_MESHSET::getSideNumberPtr(Interface &moab,EntityHandle ent) const {
  NOT_USED(moab);
  NOT_USED(ent);
  SideNumber_multiIndex::iterator miit;
  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
    boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
  ).first;
  return *miit;
  THROW_MESSAGE("not implemented");
  return boost::shared_ptr<SideNumber>();
}
RefElement_PRISM::RefElement_PRISM(
  Interface &moab,const boost::shared_ptr<RefEntity> ref_ent_ptr
):
RefElement(moab,ref_ent_ptr) {
  ErrorCode rval;
  Tag th_RefBitEdge;
  rval = moab.tag_get_handle("_RefBitEdge",th_RefBitEdge); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitEdge,&ref_ent_ptr->ent,1,(const void **)&tag_BitRefEdges); MOAB_THROW(rval);
  switch (ref_ent_ptr->getEntType()) {
    case MBPRISM:
    break;
    default:
      THROW_MESSAGE("this work only for PRISMs");
  }
  EntityHandle prism = getRefEnt();
  int num_nodes;
  const EntityHandle* conn;
  rval = moab.get_connectivity(prism,conn,num_nodes,true); MOAB_THROW(rval);
  assert(num_nodes == 6);
  for(int nn = 0;nn<6; nn++) {
    const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(conn[nn],nn,0,-1))
    );
  }
}
boost::shared_ptr<SideNumber> RefElement_PRISM::getSideNumberPtr(Interface &moab,EntityHandle ent) const {

  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  // this int is in table then return pointer
  if(miit!=side_number_table.end()) return *miit;

  // if ent is a this prism
  if(sPtr->ent == ent) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }

  // if ent is meshset
  if(moab.type_from_handle(ent)==MBENTITYSET) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }

  // use moab to get sense, side and offset
  MoABErrorCode rval;
  int side_number,sense,offset;
  rval = moab.side_number(sPtr->ent,ent,side_number,sense,offset);

  // it has to be degenerated prism, get sense from nodes topology
  if(side_number==-1 || rval != MB_SUCCESS) {

    if(moab.type_from_handle(ent)==MBVERTEX) {
      THROW_MESSAGE(
        "Huston we have problem, vertex (specified by ent) is not part of prism, that is impossible (top tip: check your prisms)"
      );
    }

    //get prism connectivity
    int num_nodes;
    const EntityHandle* conn;
    rval = moab.get_connectivity(sPtr->ent,conn,num_nodes,true); MOAB_THROW(rval);
    assert(num_nodes==6);
    //get ent connectivity
    const EntityHandle* conn_ent;
    rval = moab.get_connectivity(ent,conn_ent,num_nodes,true); MOAB_THROW(rval);

    // for(int nn = 0; nn<6;nn++) {
    //   std::cerr << conn[nn] << " ";
    // };
    // std::cerr << std::endl;
    // for(int nn = 0; nn<num_nodes;nn++) {
    //   std::cerr << conn_ent[nn] << " ";
    // }
    // std::cerr << std::endl;

    //bottom face
    EntityHandle face3[3] = { conn[0], conn[1], conn[2] };
    //top face
    EntityHandle face4[3] = { conn[3], conn[4], conn[5] };
    if(num_nodes == 3) {
      int sense_p1_map[3][3] = { {0,1,2}, {1,2,0}, {2,0,1} };
      int sense_m1_map[3][3] = { {0,2,1}, {1,0,2}, {2,1,0} };
      EntityHandle* conn0_3_ptr = std::find( face3, &face3[3], conn_ent[0] );
      if( conn0_3_ptr != &face3[3] ) {
        offset = std::distance( face3, conn0_3_ptr );
        if(
          face3[ sense_p1_map[offset][0] ] == conn_ent[0] &&
          face3[ sense_p1_map[offset][1] ] == conn_ent[1] &&
          face3[ sense_p1_map[offset][2] ] == conn_ent[2]
        )
        {
          miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
            boost::shared_ptr<SideNumber>(new SideNumber(ent,3,1,offset))
          ).first;
          return *miit;
        } else if (
          face3[ sense_m1_map[offset][0] ] == conn_ent[0] &&
          face3[ sense_m1_map[offset][1] ] == conn_ent[1] &&
          face3[ sense_m1_map[offset][2] ] == conn_ent[2]
        ) {
          miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
            boost::shared_ptr<SideNumber>(new SideNumber(ent,3,-1,offset))
          ).first;
          return *miit;
        }
      }
      EntityHandle* conn0_4_ptr = std::find( face4, &face4[3], conn_ent[0] );
      if( conn0_4_ptr != &face4[3] ) {
        offset = std::distance( face4, conn0_4_ptr );
        if(
          face4[ sense_p1_map[offset][0] ] == conn_ent[0] &&
          face4[ sense_p1_map[offset][1] ] == conn_ent[1] &&
          face4[ sense_p1_map[offset][2] ] == conn_ent[2]
        ) {
          miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
            boost::shared_ptr<SideNumber>(new SideNumber(ent,4,1,3+offset))
          ).first;
          return *miit;
        } else if (
          face4[ sense_m1_map[offset][0] ] == conn_ent[0] &&
          face4[ sense_m1_map[offset][1] ] == conn_ent[1] &&
          face4[ sense_m1_map[offset][2] ] == conn_ent[2]
        ) {
          miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
            boost::shared_ptr<SideNumber>(new SideNumber(ent,4,-1,3+offset))
          ).first;
          return *miit;
        } else {
          std::cerr << conn_ent[0] << " " << conn_ent[1] << " " << conn_ent[2] << std::endl;
          std::cerr << face3[0] << " " << face3[1] << " " << face3[2] << std::endl;
          std::cerr << face4[0] << " " << face4[1] << " " << face4[2] << std::endl;
          std::cerr << offset << std::endl;
          THROW_MESSAGE("Huston we have problem");
        }
      }
      THROW_MESSAGE("Huston we have problem");
    }

    if(num_nodes == 2) {
      {
        // Triangle edges
        EntityHandle edges[6][2] = {
          { conn[0], conn[1] } /*0*/, { conn[1], conn[2] } /*1*/, { conn[2], conn[0] } /*2*/,
          { conn[3], conn[4] } /*3+3*/, { conn[4], conn[5] } /*3+4*/, { conn[5], conn[3] } /*3+5*/
        };
        for(int ee = 0;ee<6;ee++) {
          if(
            (( conn_ent[0] == edges[ee][0] )&&( conn_ent[1] == edges[ee][1] ))||
            (( conn_ent[0] == edges[ee][1] )&&( conn_ent[1] == edges[ee][0] ))
          ) {
            side_number = ee;
            if(ee>=3) {
              side_number += 3;
              EntityHandle* conn0_4_ptr = std::find( face4, &face4[3], conn_ent[0] );
              offset = std::distance( face4, conn0_4_ptr ) + 3;
            } else {
              EntityHandle* conn0_3_ptr = std::find( face3, &face3[3], conn_ent[0] );
              offset = std::distance( face3, conn0_3_ptr );
            }
            sense = 1;
            if(( conn_ent[0] == edges[ee][1] )&&( conn_ent[1] == edges[ee][0] ))  sense = -1;
            miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
              boost::shared_ptr<SideNumber>(new SideNumber(ent,side_number,sense,offset))
            ).first;
            return *miit;
          }
        }
      }
      // {
      //   // Edges through thickness
      //   EntityHandle edges[3][2] = {
      //     { conn[0], conn[3] }, { conn[1], conn[4] }, { conn[2], conn[5] }
      //   };
      //   for(int ee = 0;ee<3;ee++) {
      //     if(
      //       (( conn_ent[0] == edges[ee][0] )&&( conn_ent[1] == edges[ee][1] ))||
      //       (( conn_ent[0] == edges[ee][1] )&&( conn_ent[1] == edges[ee][0] ))
      //     ) {
      //       side_number = 3+ee;
      //       offset = distance(conn,find(conn,&conn[6],conn_ent[0]));
      //       sense = 1;
      //       if(( conn_ent[0] == edges[ee][1] )&&( conn_ent[1] == edges[ee][0] ))  sense = -1;
      //       miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,side_number,sense,offset)).first;
      //       return const_cast<SideNumber*>(&*miit);
      //     }
      //   }
      // }
      // for(int nn = 0; nn<6;nn++) {
      //   std::cerr << conn[nn] << " ";
      // };
      // std::cerr << std::endl;
      // std::cerr << conn_ent[0] << " " << conn_ent[1] << std::endl;
      THROW_MESSAGE("Huston we have problem");
    }
    std::ostringstream sss;
    sss << "this not working: " << ent << " type: " << moab.type_from_handle(ent) << " " << MBEDGE << " " << MBTRI << std::endl;
    THROW_MESSAGE(sss.str().c_str());
  }
  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
    boost::shared_ptr<SideNumber>(new SideNumber(ent,side_number,sense,offset))
  ).first;
  return *miit;
  THROW_MESSAGE("not implemented");
  return boost::shared_ptr<SideNumber>();
}
RefElement_TET::RefElement_TET(Interface &moab,const boost::shared_ptr<RefEntity> ref_ent_ptr):
RefElement(moab,ref_ent_ptr),tag_BitRefEdges(NULL) {
  ErrorCode rval;
  Tag th_RefBitEdge;
  rval = moab.tag_get_handle("_RefBitEdge",th_RefBitEdge); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitEdge,&ref_ent_ptr->ent,1,(const void **)&tag_BitRefEdges); MOAB_THROW(rval);
  Tag th_RefType;
  switch (ref_ent_ptr->getEntType()) {
    case MBTET:
    break;
    default:
    PetscTraceBackErrorHandler(
      PETSC_COMM_WORLD,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,
      "this work only for TETs",PETSC_NULL
    );
    THROW_MESSAGE("this work only for TETs");
  }
  rval = moab.tag_get_handle("_RefType",th_RefType); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefType,&sPtr->ent,1,(const void **)&tag_type_data); MOAB_THROW(rval);
  const_cast<SideNumber_multiIndex&>(side_number_table).insert(
    boost::shared_ptr<SideNumber>(new SideNumber(sPtr->ent,0,0,0))
  );
}
boost::shared_ptr<SideNumber> RefElement_TET::getSideNumberPtr(Interface &moab,EntityHandle ent) const {
  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  if(miit!=side_number_table.end()) return *miit;
  if(sPtr->ent == ent) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }
  if(moab.type_from_handle(ent)==MBENTITYSET) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }
  ErrorCode rval;
  int side_number,sense,offset;
  rval = moab.side_number(sPtr->ent,ent,side_number,sense,offset); MOAB_THROW(rval);
  std::pair<SideNumber_multiIndex::iterator,bool> p_miit;
  p_miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
    boost::shared_ptr<SideNumber>(new SideNumber(ent,side_number,sense,offset))
  );
  miit = p_miit.first;
  if(miit->get()->ent != ent) {
    THROW_MESSAGE("this not working");
  }
  //std::cerr << side_number << " " << sense << " " << offset << std::endl;
  return *miit;
}
std::ostream& operator<<(std::ostream& os,const RefElement_TET& e) {
  os << "ref type " << e.tag_type_data[0] << " ref sub type " << e.tag_type_data[1];
  os << " ref egdes " << e.getBitRefEdges();
  os << " " << *e.sPtr;
  return os;
}

RefElement_TRI::RefElement_TRI(Interface &moab,const boost::shared_ptr<RefEntity> ref_ent_ptr):
RefElement(moab,ref_ent_ptr) {
  switch (ref_ent_ptr->getEntType()) {
    case MBTRI:
    break;
    default:
    THROW_MESSAGE("this work only for TRIs");
  }
  ErrorCode rval;
  int side_number,sense,offset;
  EntityHandle tri = getRefEnt();
  int num_nodes;
  const EntityHandle* conn;
  rval = moab.get_connectivity(tri,conn,num_nodes,true); MOAB_THROW(rval);
  for(int nn = 0;nn<3; nn++) {
    const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(conn[nn],nn,0,0))
    );
  }
  for(int ee = 0;ee<3; ee++) {
    EntityHandle edge;
    rval = moab.side_element(tri,1,ee,edge); MOAB_THROW(rval);
    rval = moab.side_number(tri,edge,side_number,sense,offset); MOAB_THROW(rval);
    const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(edge,ee,sense,offset))
    );
  }
  const_cast<SideNumber_multiIndex&>(side_number_table).insert(
    boost::shared_ptr<SideNumber>(new SideNumber(tri,0,0,0))
  );
}
boost::shared_ptr<SideNumber> RefElement_TRI::getSideNumberPtr(Interface &moab,EntityHandle ent) const {
  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  if(miit!=side_number_table.end()) return *miit;
  if(sPtr->ent == ent) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }
  if(moab.type_from_handle(ent)==MBENTITYSET) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }
  ErrorCode rval;
  int side_number,sense,offset;
  rval = moab.side_number(sPtr->ent,ent,side_number,sense,offset); MOAB_THROW(rval);
  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
    boost::shared_ptr<SideNumber>(new SideNumber(ent,side_number,sense,offset))
  ).first;
  //std::cerr << side_number << " " << sense << " " << offset << std::endl;
  return *miit;
}
std::ostream& operator<<(std::ostream& os,const RefElement_TRI& e) {
  os << *e.sPtr;
  return os;
}
RefElement_EDGE::RefElement_EDGE(Interface &moab,const boost::shared_ptr<RefEntity> ref_ent_ptr):
RefElement(moab,ref_ent_ptr) {
  switch (ref_ent_ptr->getEntType()) {
    case MBEDGE:
    break;
    default:
      THROW_MESSAGE("this work only for TRIs");
  }
}
boost::shared_ptr<SideNumber> RefElement_EDGE::getSideNumberPtr(Interface &moab,EntityHandle ent) const {
  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  if(miit!=side_number_table.end()) return *miit;
  if(sPtr->ent == ent) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }
  if(moab.type_from_handle(ent)==MBENTITYSET) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }
  ErrorCode rval;
  int side_number,sense,offset;
  rval = moab.side_number(sPtr->ent,ent,side_number,sense,offset); MOAB_THROW(rval);
  miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
    boost::shared_ptr<SideNumber>(new SideNumber(ent,side_number,sense,offset))
  ).first;
  //std::cerr << side_number << " " << sense << " " << offset << std::endl;
  return *miit;
}
std::ostream& operator<<(std::ostream& os,const RefElement_EDGE& e) {
  os << *e.sPtr;
  return os;
}
RefElement_VERTEX::RefElement_VERTEX(Interface &moab,const boost::shared_ptr<RefEntity> ref_ent_ptr):
RefElement(moab,ref_ent_ptr) {
  switch (ref_ent_ptr->getEntType()) {
    case MBVERTEX:
    break;
    default:
      THROW_MESSAGE("this works only for TRIs");
  }
}
boost::shared_ptr<SideNumber> RefElement_VERTEX::getSideNumberPtr(
  Interface &moab,EntityHandle ent
) const {
  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  if(miit!=side_number_table.end()) return *miit;
  if(sPtr->ent == ent) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }
  if(moab.type_from_handle(ent)==MBENTITYSET) {
    miit = const_cast<SideNumber_multiIndex&>(side_number_table).insert(
      boost::shared_ptr<SideNumber>(new SideNumber(ent,0,0,0))
    ).first;
    return *miit;
  }
  THROW_MESSAGE("no side entity for vertex if its is not an vertex itself");
  return boost::shared_ptr<SideNumber>();
}
std::ostream& operator<<(std::ostream& os,const RefElement_VERTEX& e) {
  os << *e.sPtr;
  return os;
}

PetscErrorCode DefaultElementAdjacency::defaultVertex(
  Interface &moab,const Field &field_ptr,const EntFiniteElement &fe_ptr,Range &adjacency
) {
  PetscFunctionBegin;
  MoABErrorCode rval;
  switch (field_ptr.getSpace()) {
    case H1:
    adjacency.insert(fe_ptr.getEnt());
    break;
    case NOFIELD:
    {
      Range ents;
      rval = moab.get_entities_by_handle(field_ptr.getMeshset(),ents,false); CHKERRQ_MOAB(rval);
      adjacency.merge(ents);
      for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
        const_cast<SideNumber_multiIndex&>(fe_ptr.getSideNumberTable()).insert(
          boost::shared_ptr<SideNumber>(new SideNumber(*eit,-1,0,0))
        );
      }
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this field is not implemented for VERTEX finite element");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DefaultElementAdjacency::defaultEdge(
  Interface &moab,const Field &field_ptr,const EntFiniteElement &fe_ptr,Range &adjacency
) {
  PetscFunctionBegin;
  ErrorCode rval;
  EntityHandle fe_ent = fe_ptr.getEnt();
  Range nodes;
  switch (field_ptr.getSpace()) {
    case H1:
    //moab.get_connectivity(&fe_ent,1,nodes,true);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(&fe_ent,1,0,false,nodes,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(&fe_ent,1,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(&fe_ent,1,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    adjacency.insert(nodes.begin(),nodes.end());
    adjacency.insert(fe_ent);
    break;
    case HCURL:
    adjacency.insert(fe_ent);
    break;
    case NOFIELD:
    {
      Range ents;
      rval = moab.get_entities_by_handle(field_ptr.getMeshset(),ents,false); CHKERRQ_MOAB(rval);
      adjacency.merge(ents);
      for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
        const_cast<SideNumber_multiIndex&>(fe_ptr.getSideNumberTable()).insert(
          boost::shared_ptr<SideNumber>(new SideNumber(*eit,-1,0,0))
        );
      }
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this field is not implemented for EDGE finite element");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DefaultElementAdjacency::defaultTri(
  Interface &moab,const Field &field_ptr,const EntFiniteElement &fe_ptr,Range &adjacency
) {
  PetscFunctionBegin;
  ErrorCode rval;
  Range nodes,edges;
  EntityHandle fe_ent = fe_ptr.getEnt();
  switch (field_ptr.getSpace()) {
    case H1:
    //moab.get_connectivity(&fe_ent,1,nodes,true);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(&fe_ent,1,0,false,nodes,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(&fe_ent,1,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(&fe_ent,1,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    adjacency.insert(nodes.begin(),nodes.end());
    rval = moab.get_adjacencies(&fe_ent,1,1,false,edges); CHKERRQ_MOAB(rval);
    adjacency.insert(edges.begin(),edges.end());
    for(Range::iterator eeit = edges.begin();eeit!=edges.end();eeit++) {
      fe_ptr.getSideNumberPtr(moab,*eeit);
    }
    //add faces
    adjacency.insert(fe_ent);
    break;
    case HCURL:
    rval = moab.get_adjacencies(&fe_ent,1,1,false,edges); CHKERRQ_MOAB(rval);
    adjacency.insert(edges.begin(),edges.end());
    for(Range::iterator eeit = edges.begin();eeit!=edges.end();eeit++) {
      fe_ptr.getSideNumberPtr(moab,*eeit);
    }
    //add faces
    adjacency.insert(fe_ent);
    break;
    case HDIV:
    adjacency.insert(fe_ent);
    break;
    case NOFIELD:
    {
      Range ents;
      rval = moab.get_entities_by_handle(field_ptr.getMeshset(),ents,false); CHKERRQ_MOAB(rval);
      adjacency.merge(ents);
      for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
        const_cast<SideNumber_multiIndex&>(fe_ptr.getSideNumberTable()).insert(
          boost::shared_ptr<SideNumber>(new SideNumber(*eit,-1,0,0))
        );
      }
    }
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
PetscErrorCode DefaultElementAdjacency::defaultTet(
  Interface &moab,const Field &field_ptr,const EntFiniteElement &fe_ptr,Range &adjacency
) {
  PetscFunctionBegin;
  ErrorCode rval;
  Range nodes,edges,faces;
  EntityHandle fe_ent = fe_ptr.getEnt();
  switch (field_ptr.getSpace()) {
    case H1:
    //moab.get_connectivity(&fe_ent,1,nodes,true);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(&fe_ent,1,0,false,nodes,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(&fe_ent,1,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(&fe_ent,1,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    if(nodes.size()<4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"TET has at least 4 adjacent nodes; it can has more if user add more adjacencies");
    }
    adjacency.insert(nodes.begin(),nodes.end());
    case HCURL:
    rval = moab.get_adjacencies(&fe_ent,1,1,false,edges); CHKERRQ_MOAB(rval);
    adjacency.insert(edges.begin(),edges.end());
    for(Range::iterator eeit = edges.begin();eeit!=edges.end();eeit++) {
      fe_ptr.getSideNumberPtr(moab,*eeit);
    }
    case HDIV:
    rval = moab.get_adjacencies(&fe_ent,1,2,false,faces); CHKERRQ_MOAB(rval);
    adjacency.insert(faces.begin(),faces.end());
    for(Range::iterator fit = faces.begin();fit!=faces.end();fit++) {
      fe_ptr.getSideNumberPtr(moab,*fit);
    }
    case L2:
    adjacency.insert(fe_ent);
    break;
    case NOFIELD:
    {
      Range ents;
      rval = moab.get_entities_by_handle(field_ptr.getMeshset(),ents,false); CHKERRQ_MOAB(rval);
      adjacency.merge(ents);
      for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
        const_cast<SideNumber_multiIndex&>(fe_ptr.getSideNumberTable()).insert(
          boost::shared_ptr<SideNumber>(new SideNumber(*eit,-1,0,0))
        );
      }
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this field is not implemented for TRI finite element");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DefaultElementAdjacency::defaultPrism(
  Interface &moab,const Field &field_ptr,const EntFiniteElement &fe_ptr,Range &adjacency
) {
  PetscFunctionBegin;
  ErrorCode rval;
  EntityHandle fe_ent = fe_ptr.getEnt();
  Range nodes;
  //initialize side sets
  try {
    EntityHandle prism = fe_ent;
    EntityHandle face_side3,face_side4;
    rval = moab.side_element(prism,2,3,face_side3); CHKERRQ_MOAB(rval);
    rval = moab.side_element(prism,2,4,face_side4); CHKERRQ_MOAB(rval);
    fe_ptr.getRefElement()->getSideNumberPtr(moab,face_side3);
    fe_ptr.getRefElement()->getSideNumberPtr(moab,face_side4);
    for(int qq = 0;qq<3;qq++) {
      EntityHandle quad = 0;
      rval = moab.side_element(prism,2,qq,quad);
      if(rval != MB_SUCCESS || quad == 0) continue;
      int side_number,sense,offset;
      rval = moab.side_number(prism,quad,side_number,sense,offset);
      if(side_number==-1 || rval != MB_SUCCESS) continue;
      const_cast<SideNumber_multiIndex&>(fe_ptr.getSideNumberTable()).insert(
        boost::shared_ptr<SideNumber>(new SideNumber(quad,side_number,sense,offset))
      );
    }
    int ee = 0;
    for(;ee<3;ee++) {
      EntityHandle edge = 0;
      rval = moab.side_element(prism,1,ee,edge); CHKERRQ_MOAB(rval);
      boost::shared_ptr<SideNumber> side_ptr =
      fe_ptr.getRefElement()->getSideNumberPtr(moab,edge);
      if(side_ptr->side_number!=ee) {
        SETERRQ1(PETSC_COMM_SELF,1,"data insistency for edge %d",ee);
      }
      rval = moab.side_element(prism,1,6+ee,edge); CHKERRQ_MOAB(rval);
      side_ptr = fe_ptr.getRefElement()->getSideNumberPtr(moab,edge);
      if(side_ptr->side_number!=ee+6) {
        if(side_ptr->side_number!=ee) {
          SETERRQ1(PETSC_COMM_SELF,1,"data insistency for edge %d",ee);
        } else {
          side_ptr->brother_side_number = ee+6;
        }
      }
    }
    for(;ee<6;ee++) {
      EntityHandle edge = 0;
      rval = moab.side_element(prism,1,ee,edge);
      if(rval != MB_SUCCESS || edge == 0) continue;
      int side_number,sense,offset;
      rval = moab.side_number(prism,edge,side_number,sense,offset);
      if(side_number==-1 || rval != MB_SUCCESS) continue;
      const_cast<SideNumber_multiIndex&>(fe_ptr.getSideNumberTable()).insert(
        boost::shared_ptr<SideNumber>(new SideNumber(edge,side_number,sense,offset))
      );
    }
    int nn = 0;
    for(;nn<3;nn++) {
      EntityHandle node;
      rval = moab.side_element(prism,0,nn,node); CHKERRQ_MOAB(rval);
      boost::shared_ptr<SideNumber> side_ptr = fe_ptr.getRefElement()->getSideNumberPtr(moab,node);
      if(side_ptr->side_number!=nn) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data insistency for node %d",nn);
      }
      rval = moab.side_element(prism,0,nn+3,node); CHKERRQ_MOAB(rval);
      side_ptr = fe_ptr.getRefElement()->getSideNumberPtr(moab,node);
      if(side_ptr->side_number!=nn+3) {
        if(side_ptr->side_number!=nn) {
          SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data insistency for node %d",nn);
        } else {
          side_ptr->brother_side_number = nn+3;
        }
      }
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  //get adjacencies
  SideNumber_multiIndex &side_table = fe_ptr.getRefElement()->getSideNumberTable();
  switch(field_ptr.getSpace()) {
    case H1:
    //moab.get_connectivity(&fe_ent,1,nodes,true);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(&fe_ent,1,0,false,nodes,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(&fe_ent,1,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(&fe_ent,1,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    adjacency.insert(nodes.begin(),nodes.end());
    case HCURL: {
      SideNumber_multiIndex::nth_index<2>::type::iterator siit,hi_siit;
      siit = side_table.get<2>().lower_bound(MBEDGE);
      hi_siit = side_table.get<2>().upper_bound(MBEDGE);
      for(;siit!=hi_siit;siit++) adjacency.insert(siit->get()->ent);
    }
    case HDIV: {
      SideNumber_multiIndex::nth_index<2>::type::iterator siit,hi_siit;
      siit = side_table.get<2>().lower_bound(MBTRI);
      hi_siit = side_table.get<2>().upper_bound(MBTRI);
      for(;siit!=hi_siit;siit++) adjacency.insert(siit->get()->ent);
      siit = side_table.get<2>().lower_bound(MBQUAD);
      hi_siit = side_table.get<2>().upper_bound(MBQUAD);
      for(;siit!=hi_siit;siit++) adjacency.insert(siit->get()->ent);
    }
    case L2:
    adjacency.insert(fe_ent);
    break;
    case NOFIELD:
    {
      Range ents;
      rval = moab.get_entities_by_handle(field_ptr.getMeshset(),ents,false); CHKERRQ_MOAB(rval);
      adjacency.merge(ents);
      for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
        const_cast<SideNumber_multiIndex&>(fe_ptr.getSideNumberTable()).insert(
          boost::shared_ptr<SideNumber>(new SideNumber(*eit,-1,0,0))
        );
      }
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"this field is not implemented for TRI finite element");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode DefaultElementAdjacency::defaultMeshset(
  Interface &moab,const Field &field_ptr,const EntFiniteElement &fe_ptr,Range &adjacency
) {
  PetscFunctionBegin;
  ErrorCode rval;
  EntityHandle fe_ent = fe_ptr.getEnt();
  //get all meshsets in finite element meshset
  Range ent_ents_meshset;
  rval = moab.get_entities_by_type(fe_ent,MBENTITYSET,ent_ents_meshset,false); CHKERRQ_MOAB(rval);
  //resolve recursively all ents in the meshset
  Range ent_ents;
  rval = moab.get_entities_by_handle(fe_ent,ent_ents,true); CHKERRQ_MOAB(rval);
  switch (field_ptr.getSpace()) {
    case H1:
    adjacency.merge(ent_ents.subset_by_type(MBVERTEX));
    case HCURL:
    adjacency.merge(ent_ents.subset_by_type(MBEDGE));
    case HDIV:
    adjacency.merge(ent_ents.subset_by_type(MBTRI));
    case L2:
    adjacency.merge(ent_ents.subset_by_type(MBTET));
    break;
    case NOFIELD:
    {
      Range ents;
      rval = moab.get_entities_by_handle(field_ptr.getMeshset(),ents,false); CHKERRQ_MOAB(rval);
      adjacency.merge(ents);
      for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
        const_cast<SideNumber_multiIndex&>(fe_ptr.getSideNumberTable()).insert(
          boost::shared_ptr<SideNumber>(new SideNumber(*eit,-1,0,0))
        );
      }
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  PetscFunctionReturn(0);
}

//FiniteElement
FiniteElement::FiniteElement(Interface &moab,const EntityHandle _meshset): meshset(_meshset) {
  ErrorCode rval;
  Tag th_FEId;
  rval = moab.tag_get_handle("_FEId",th_FEId); CHKERR_MOAB(rval);
  rval = moab.tag_get_by_ptr(th_FEId,&meshset,1,(const void **)&tag_id_data); CHKERR_MOAB(rval);
  Tag th_FEName;
  rval = moab.tag_get_handle("_FEName",th_FEName); CHKERR_MOAB(rval);
  rval = moab.tag_get_by_ptr(th_FEName,&meshset,1,(const void **)&tag_name_data,&tag_name_size); CHKERR_MOAB(rval);
  Tag th_FEIdCol,th_FEIdRow,th_FEIdData;
  rval = moab.tag_get_handle("_FEIdCol",th_FEIdCol); CHKERR_MOAB(rval);
  rval = moab.tag_get_by_ptr(th_FEIdCol,&meshset,1,(const void **)&tag_BitFieldId_col_data); CHKERR_MOAB(rval);
  rval = moab.tag_get_handle("_FEIdRow",th_FEIdRow); CHKERR_MOAB(rval);
  rval = moab.tag_get_by_ptr(th_FEIdRow,&meshset,1,(const void **)&tag_BitFieldId_row_data); CHKERR_MOAB(rval);
  rval = moab.tag_get_handle("_FEIdData",th_FEIdData); CHKERR_MOAB(rval);
  rval = moab.tag_get_by_ptr(th_FEIdData,&meshset,1,(const void **)&tag_BitFieldId_data); CHKERR_MOAB(rval);
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

std::ostream& operator<<(std::ostream& os,const FiniteElement& e) {
    os << "id " << e.getId() << " name " << e.getNameRef() << " f_id_row " << e.getBitFieldIdRow()
    << " f_id_col " << e.getBitFieldIdCol() << " BitFEId_data " << e.getBitFieldIdData();
    return os;
}

void MoFEMFiniteElement_col_change_bit_add::operator()(boost::shared_ptr<FiniteElement> &fe) {
  *((BitFieldId*)(fe->tag_BitFieldId_col_data)) |= fIdCol;
}

void MoFEMFiniteElement_row_change_bit_add::operator()(boost::shared_ptr<FiniteElement> &fe) {
  *((BitFieldId*)(fe->tag_BitFieldId_row_data)) |= fIdRow;
}

void MoFEMFiniteElement_change_bit_add::operator()(boost::shared_ptr<FiniteElement> &fe) {
  *((BitFieldId*)(fe->tag_BitFieldId_data)) |= fIdData;
}

void MoFEMFiniteElement_col_change_bit_off::operator()(boost::shared_ptr<FiniteElement> &fe) {
  *((BitFieldId*)(fe->tag_BitFieldId_col_data)) &= fIdCol.flip();
}

void MoFEMFiniteElement_row_change_bit_off::operator()(boost::shared_ptr<FiniteElement> &fe) {
  *((BitFieldId*)(fe->tag_BitFieldId_row_data)) &= fIdRow.flip();
}

void MoFEMFiniteElement_change_bit_off::operator()(boost::shared_ptr<FiniteElement> &fe) {
  *((BitFieldId*)(fe->tag_BitFieldId_data)) &= fIdData.flip();
}

//FiniteElement data
EntFiniteElement::EntFiniteElement(
  Interface &moab,
  const boost::shared_ptr<RefElement> ref_finite_element,
  const boost::shared_ptr<FiniteElement> fe_ptr
):
interface_FiniteElement<FiniteElement>(fe_ptr),
interface_RefElement<RefElement>(ref_finite_element),
row_dof_view(boost::shared_ptr<DofEntity_multiIndex_uid_view>(new DofEntity_multiIndex_uid_view)),
col_dof_view(boost::shared_ptr<DofEntity_multiIndex_uid_view>(new DofEntity_multiIndex_uid_view)),
data_dof_view(boost::shared_ptr<DofEntity_multiIndex_uid_view>(new DofEntity_multiIndex_uid_view)) {

  //get finite element entity
  global_uid =  getGlobalUniqueIdCalculate();
  //add ents to meshset
  //EntityHandle meshset = getMeshset();
  //EntityHandle ent = getEnt();
  //ierr = moab.add_entities(meshset,&ent,1); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

std::ostream& operator<<(std::ostream& os,const EntFiniteElement& e) {
  os << *e.sFePtr << std::endl;
  os << *e.sPtr << std::endl;
  os << "row dof_uids ";
  DofEntity_multiIndex_uid_view::iterator rit;
  rit = e.row_dof_view->begin();
  for(;rit!=e.row_dof_view->end();rit++) {
    os << (*rit)->getGlobalUniqueId() << " ";
  }
  os << "col dof_uids ";
  DofEntity_multiIndex_uid_view::iterator cit;
  cit = e.col_dof_view->begin();
  for(;cit!=e.col_dof_view->end();cit++) {
    os << (*cit)->getGlobalUniqueId() << " ";
  }
  os << "data dof_uids ";
  DofEntity_multiIndex_uid_view::iterator dit;
  dit = e.data_dof_view->begin();
  for(;dit!=e.data_dof_view->end();dit++) {
    os << (*dit)->getGlobalUniqueId() << " ";
  }
  return os;
}

template <typename MOFEM_DOFS,typename MOFEM_DOFS_VIEW>
static PetscErrorCode get_fe_dof_view(
  const DofEntity_multiIndex_uid_view &fe_dofs_view,
  const MOFEM_DOFS &mofem_dofs,
  MOFEM_DOFS_VIEW &mofem_dofs_view,
  const int operation_type
) {
  PetscFunctionBegin;
  typename boost::multi_index::index<MOFEM_DOFS,Unique_mi_tag>::type::iterator mofem_it,mofem_it_end;
  DofEntity_multiIndex_uid_view::iterator it,it_end;
  if(operation_type==moab::Interface::UNION) {
    mofem_it = mofem_dofs.template get<Unique_mi_tag>().begin();
    mofem_it_end = mofem_dofs.template get<Unique_mi_tag>().end();
    it = fe_dofs_view.begin();
    it_end = fe_dofs_view.end();
    for(;it!=it_end;it++) {
      const GlobalUId &global_uid = (*it)->getGlobalUniqueId();
      if(mofem_it != mofem_it_end) {
        if((*mofem_it)->getGlobalUniqueId() != global_uid) {
          mofem_it = mofem_dofs.template get<Unique_mi_tag>().find(global_uid);
        }  // else lucky guess
      } else {
        mofem_it = mofem_dofs.template get<Unique_mi_tag>().find(global_uid);
      }
      if(mofem_it != mofem_it_end) {
        mofem_dofs_view.insert(mofem_dofs_view.end(),*mofem_it);
        mofem_it++;
      }
    }
  } else {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EntFiniteElement::getRowDofView(
  const DofEntity_multiIndex &dofs,DofEntity_multiIndex_active_view &dofs_view,
  const int operation_type
) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_dof_view(*row_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntFiniteElement::getColDofView(
  const DofEntity_multiIndex &dofs,DofEntity_multiIndex_active_view &dofs_view,
  const int operation_type
) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_dof_view(*col_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntFiniteElement::getDataDofView(
  const DofEntity_multiIndex &dofs,DofEntity_multiIndex_active_view &dofs_view,
  const int operation_type
) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_dof_view(*data_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EntFiniteElement::getRowDofView(
  const DofEntity_multiIndex &dofs,DofEntity_multiIndex_uid_view &dofs_view,
  const int operation_type
) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_dof_view(*row_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntFiniteElement::getColDofView(
  const DofEntity_multiIndex &dofs,DofEntity_multiIndex_uid_view &dofs_view,
  const int operation_type
) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_dof_view(*col_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EntFiniteElement::getRowDofView(
  const NumeredDofEntity_multiIndex &dofs,NumeredDofEntity_multiIndex_uid_view_ordered &dofs_view,
  const int operation_type
) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_dof_view(*row_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EntFiniteElement::getColDofView(
  const NumeredDofEntity_multiIndex &dofs,NumeredDofEntity_multiIndex_uid_view_ordered &dofs_view,
  const int operation_type
) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_dof_view(*col_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode EntFiniteElement::getRowDofView(
  const NumeredDofEntity_multiIndex &dofs,NumeredDofEntity_multiIndex_uid_view_hashed &dofs_view,
  const int operation_type
) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_dof_view(*row_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EntFiniteElement::getColDofView(
    const NumeredDofEntity_multiIndex &dofs,NumeredDofEntity_multiIndex_uid_view_hashed &dofs_view,
    const int operation_type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = get_fe_dof_view(*col_dof_view,dofs,dofs_view,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NumeredEntFiniteElement::getRowDofsByPetscGlobalDofIdx(DofIdx idx,const FENumeredDofEntity **dof_ptr) const {
  PetscFunctionBegin;
  FENumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = rows_dofs->get<PetscGlobalIdx_mi_tag>().find(idx);
  if(dit == rows_dofs->get<PetscGlobalIdx_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"dof which index < %d > not found",idx);
  }
  *dof_ptr = &*dit->get();
  PetscFunctionReturn(0);
}

PetscErrorCode NumeredEntFiniteElement::getColDofsByPetscGlobalDofIdx(DofIdx idx,const FENumeredDofEntity **dof_ptr) const {
  PetscFunctionBegin;
  FENumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit;
  dit = rows_dofs->get<PetscGlobalIdx_mi_tag>().find(idx);
  if(dit == rows_dofs->get<PetscGlobalIdx_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"dof which index < %d > not found",idx);
  }
  *dof_ptr = &*dit->get();
  PetscFunctionReturn(0);
}


}
