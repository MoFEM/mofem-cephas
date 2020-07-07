/** \file FEMultiIndices.cpp
 * \brief Multi-index containers for finite elements
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

namespace MoFEM {

// ref moab ent
BitRefEdges MoFEM::RefElement::DummyBitRefEdges = BitRefEdges(0);

const boost::shared_ptr<SideNumber> RefElement::nullSideNumber =
    boost::shared_ptr<SideNumber>();

// ref moab FiniteElement
RefElement::RefElement(const boost::shared_ptr<RefEntity> &ref_ents_ptr)
    : interface_RefEntity<RefEntity>(ref_ents_ptr) {}

std::ostream &operator<<(std::ostream &os, const RefElement &e) {
  os << " ref egdes " << e.getBitRefEdges();
  os << " " << *(e.sPtr);
  return os;
}

RefElement_MESHSET::RefElement_MESHSET(
    const boost::shared_ptr<RefEntity> &ref_ents_ptr)
    : RefElement(ref_ents_ptr) {
  switch (ref_ents_ptr->getEntType()) {
  case MBENTITYSET:
    break;
  default:
    THROW_MESSAGE("this work only for MESHSETs");
  }
}
const boost::shared_ptr<SideNumber> &
RefElement_MESHSET::getSideNumberPtr(const EntityHandle ent) const {
  NOT_USED(ent);
  SideNumber_multiIndex::iterator miit;
  miit =
      const_cast<SideNumber_multiIndex &>(side_number_table)
          .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, 0, 0, 0)))
          .first;
  return *miit;
}
RefElement_PRISM::RefElement_PRISM(
    const boost::shared_ptr<RefEntity> &ref_ents_ptr)
    : RefElement(ref_ents_ptr) {
  Tag th_RefBitEdge;
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  rval = moab.tag_get_handle("_RefBitEdge", th_RefBitEdge);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitEdge, &ref_ents_ptr->ent, 1,
                             (const void **)&tag_BitRefEdges);
  MOAB_THROW(rval);
  switch (ref_ents_ptr->getEntType()) {
  case MBPRISM:
    break;
  default:
    THROW_MESSAGE("this work only for PRISMs");
  }
  EntityHandle prism = getRefEnt();
  int num_nodes;
  const EntityHandle *conn;
  rval = moab.get_connectivity(prism, conn, num_nodes, true);
  MOAB_THROW(rval);
  assert(num_nodes == 6);
  for (int nn = 0; nn != 6; ++nn) {
    const_cast<SideNumber_multiIndex &>(side_number_table)
        .insert(
            boost::shared_ptr<SideNumber>(new SideNumber(conn[nn], nn, 0, -1)));
  }
  // Range face_side3, face_side4;
  // CHKERR moab.get_adjacencies(conn, 3, 2, true, face_side3);
  // CHKERR moab.get_adjacencies(&conn[3], 3, 2, true, face_side4);
  // if (face_side3.size() != 1)
  //   THROW_MESSAGE("prism don't have side face 3");
  // if (face_side4.size() != 1)
  //   THROW_MESSAGE("prims don't have side face 4");
  // getSideNumberPtr(*face_side3.begin());
  // getSideNumberPtr(*face_side4.begin());
}
const boost::shared_ptr<SideNumber> &
RefElement_PRISM::getSideNumberPtr(const EntityHandle ent) const {

  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;

  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  // this int is in table then return pointer
  if (miit != side_number_table.end())
    return *miit;

  // if ent is a this prism
  if (sPtr->ent == ent) {
    miit =
        const_cast<SideNumber_multiIndex &>(side_number_table)
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, 0, 0, 0)))
            .first;
    return *miit;
  }

  // if ent is meshset
  if (moab.type_from_handle(ent) == MBENTITYSET) {
    miit =
        const_cast<SideNumber_multiIndex &>(side_number_table)
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, 0, 0, 0)))
            .first;
    return *miit;
  }

  // use moab to get sense, side and offset

  int side_number, sense, offset;
  rval = moab.side_number(sPtr->ent, ent, side_number, sense, offset);

  // it has to be degenerated prism, get sense from nodes topology
  if (side_number == -1 || rval != MB_SUCCESS) {

    if (moab.type_from_handle(ent) == MBVERTEX) {
      THROW_MESSAGE("Huston we have problem, vertex (specified by ent) is not "
                    "part of prism, that is impossible (top tip: check your "
                    "prisms)");
    }

    // get prism connectivity
    int num_nodes;
    const EntityHandle *conn;
    rval = moab.get_connectivity(sPtr->ent, conn, num_nodes, true);
    MOAB_THROW(rval);
    assert(num_nodes == 6);
    // get ent connectivity
    const EntityHandle *conn_ent;
    rval = moab.get_connectivity(ent, conn_ent, num_nodes, true);
    MOAB_THROW(rval);

    // for(int nn = 0; nn<6;nn++) {
    //   std::cerr << conn[nn] << " ";
    // };
    // std::cerr << std::endl;
    // for(int nn = 0; nn<num_nodes;nn++) {
    //   std::cerr << conn_ent[nn] << " ";
    // }
    // std::cerr << std::endl;

    // bottom face
    EntityHandle face3[3] = {conn[0], conn[1], conn[2]};
    // top face
    EntityHandle face4[3] = {conn[3], conn[4], conn[5]};
    if (num_nodes == 3) {
      int sense_p1_map[3][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
      int sense_m1_map[3][3] = {{0, 2, 1}, {1, 0, 2}, {2, 1, 0}};
      EntityHandle *conn0_3_ptr = std::find(face3, &face3[3], conn_ent[0]);
      if (conn0_3_ptr != &face3[3]) {
        offset = std::distance(face3, conn0_3_ptr);
        if (face3[sense_p1_map[offset][0]] == conn_ent[0] &&
            face3[sense_p1_map[offset][1]] == conn_ent[1] &&
            face3[sense_p1_map[offset][2]] == conn_ent[2]) {
          miit = const_cast<SideNumber_multiIndex &>(side_number_table)
                     .insert(boost::shared_ptr<SideNumber>(
                         new SideNumber(ent, 3, 1, offset)))
                     .first;
          return *miit;
        } else if (face3[sense_m1_map[offset][0]] == conn_ent[0] &&
                   face3[sense_m1_map[offset][1]] == conn_ent[1] &&
                   face3[sense_m1_map[offset][2]] == conn_ent[2]) {
          miit = const_cast<SideNumber_multiIndex &>(side_number_table)
                     .insert(boost::shared_ptr<SideNumber>(
                         new SideNumber(ent, 3, -1, offset)))
                     .first;
          return *miit;
        }
      }
      EntityHandle *conn0_4_ptr = std::find(face4, &face4[3], conn_ent[0]);
      if (conn0_4_ptr != &face4[3]) {
        offset = std::distance(face4, conn0_4_ptr);
        if (face4[sense_p1_map[offset][0]] == conn_ent[0] &&
            face4[sense_p1_map[offset][1]] == conn_ent[1] &&
            face4[sense_p1_map[offset][2]] == conn_ent[2]) {
          miit = const_cast<SideNumber_multiIndex &>(side_number_table)
                     .insert(boost::shared_ptr<SideNumber>(
                         new SideNumber(ent, 4, 1, 3 + offset)))
                     .first;
          return *miit;
        } else if (face4[sense_m1_map[offset][0]] == conn_ent[0] &&
                   face4[sense_m1_map[offset][1]] == conn_ent[1] &&
                   face4[sense_m1_map[offset][2]] == conn_ent[2]) {
          miit = const_cast<SideNumber_multiIndex &>(side_number_table)
                     .insert(boost::shared_ptr<SideNumber>(
                         new SideNumber(ent, 4, -1, 3 + offset)))
                     .first;
          return *miit;
        } else {
          std::cerr << conn_ent[0] << " " << conn_ent[1] << " " << conn_ent[2]
                    << std::endl;
          std::cerr << face3[0] << " " << face3[1] << " " << face3[2]
                    << std::endl;
          std::cerr << face4[0] << " " << face4[1] << " " << face4[2]
                    << std::endl;
          std::cerr << offset << std::endl;
          THROW_MESSAGE("Huston we have problem");
        }
      }
      THROW_MESSAGE("Huston we have problem");
    }

    if (num_nodes == 2) {
      {
        // Triangle edges
        EntityHandle edges[6][2] = {
            {conn[0], conn[1]} /*0*/,   {conn[1], conn[2]} /*1*/,
            {conn[2], conn[0]} /*2*/,   {conn[3], conn[4]} /*3+3*/,
            {conn[4], conn[5]} /*3+4*/, {conn[5], conn[3]} /*3+5*/
        };
        for (int ee = 0; ee < 6; ee++) {
          if (((conn_ent[0] == edges[ee][0]) &&
               (conn_ent[1] == edges[ee][1])) ||
              ((conn_ent[0] == edges[ee][1]) &&
               (conn_ent[1] == edges[ee][0]))) {
            side_number = ee;
            if (ee >= 3) {
              side_number += 3;
              EntityHandle *conn0_4_ptr =
                  std::find(face4, &face4[3], conn_ent[0]);
              offset = std::distance(face4, conn0_4_ptr) + 3;
            } else {
              EntityHandle *conn0_3_ptr =
                  std::find(face3, &face3[3], conn_ent[0]);
              offset = std::distance(face3, conn0_3_ptr);
            }
            sense = 1;
            if ((conn_ent[0] == edges[ee][1]) && (conn_ent[1] == edges[ee][0]))
              sense = -1;
            miit = const_cast<SideNumber_multiIndex &>(side_number_table)
                       .insert(boost::shared_ptr<SideNumber>(
                           new SideNumber(ent, side_number, sense, offset)))
                       .first;
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
      //       (( conn_ent[0] == edges[ee][0] )&&( conn_ent[1] == edges[ee][1]
      //       ))||
      //       (( conn_ent[0] == edges[ee][1] )&&( conn_ent[1] == edges[ee][0]
      //       ))
      //     ) {
      //       side_number = 3+ee;
      //       offset = std::distance(conn,find(conn,&conn[6],conn_ent[0]));
      //       sense = 1;
      //       if(( conn_ent[0] == edges[ee][1] )&&( conn_ent[1] == edges[ee][0]
      //       ))  sense = -1; miit =
      //       const_cast<SideNumber_multiIndex&>(side_number_table).insert(SideNumber(ent,side_number,sense,offset)).first;
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
    sss << "this not working: " << ent
        << " type: " << moab.type_from_handle(ent) << " " << MBEDGE << " "
        << MBTRI << std::endl;
    THROW_MESSAGE(sss.str().c_str());
  }
  miit = const_cast<SideNumber_multiIndex &>(side_number_table)
             .insert(boost::shared_ptr<SideNumber>(
                 new SideNumber(ent, side_number, sense, offset)))
             .first;
  return *miit;
  THROW_MESSAGE("not implemented");
  return nullSideNumber;
}

RefElement_TET::RefElement_TET(const boost::shared_ptr<RefEntity> &ref_ents_ptr)
    : RefElement(ref_ents_ptr), tag_BitRefEdges(NULL) {
  Tag th_RefBitEdge;
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  rval = moab.tag_get_handle("_RefBitEdge", th_RefBitEdge);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_RefBitEdge, &ref_ents_ptr->ent, 1,
                             (const void **)&tag_BitRefEdges);
  MOAB_THROW(rval);
  switch (ref_ents_ptr->getEntType()) {
  case MBTET:
    break;
  default:
    THROW_MESSAGE("this work only for TETs");
  }
  const_cast<SideNumber_multiIndex &>(side_number_table)
      .insert(boost::make_shared<SideNumber>(sPtr->ent, 0, 0, 0));
}

const boost::shared_ptr<SideNumber> &
RefElement_TET::getSideNumberPtr(const EntityHandle ent) const {
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  auto miit = side_number_table.find(ent);
  if (miit != side_number_table.end())
    return *miit;
  if (sPtr->ent == ent) {
    miit = const_cast<SideNumber_multiIndex &>(side_number_table)
               .insert(boost::make_shared<SideNumber>(ent, 0, 0, 0))
               .first;
    return *miit;
  }
  if (moab.type_from_handle(ent) == MBENTITYSET) {
    miit = const_cast<SideNumber_multiIndex &>(side_number_table)
               .insert(boost::make_shared<SideNumber>(ent, 0, 0, 0))
               .first;
    return *miit;
  }
  int side_number, sense, offset;
  rval = moab.side_number(sPtr->ent, ent, side_number, sense, offset);
  MOAB_THROW(rval);
  auto p_miit = const_cast<SideNumber_multiIndex &>(side_number_table)
                    .insert(boost::make_shared<SideNumber>(ent, side_number,
                                                           sense, offset));
  miit = p_miit.first;
  if (miit->get()->ent != ent)
    THROW_MESSAGE("this not working");

  return *miit;
}
std::ostream &operator<<(std::ostream &os, const RefElement_TET &e) {
  os << "ref type " << e.tag_type_data[0] << " ref sub type "
     << e.tag_type_data[1];
  os << " ref egdes " << e.getBitRefEdges();
  os << " " << *e.sPtr;
  return os;
}

RefElementFace::RefElementFace(const boost::shared_ptr<RefEntity> &ref_ents_ptr)
    : RefElement(ref_ents_ptr) {

  int nb_nodes = 0;
  int nb_edges = 0;
  switch (ref_ents_ptr->getEntType()) {
  case MBTRI:
    nb_nodes = nb_edges = 3;
    break;
  case MBQUAD:
    nb_nodes = nb_edges = 4;
    break;
  default:
    THROW_MESSAGE("this works only for TRIs and QUADs");
  }
  int side_number, sense, offset;
  EntityHandle tri = getRefEnt();
  int num_nodes;
  const EntityHandle *conn;
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  rval = moab.get_connectivity(tri, conn, num_nodes, true);
  MOAB_THROW(rval);
  for (int nn = 0; nn < nb_nodes; nn++) {
    const_cast<SideNumber_multiIndex &>(side_number_table)
        .insert(
            boost::shared_ptr<SideNumber>(new SideNumber(conn[nn], nn, 0, 0)));
  }
  for (int ee = 0; ee < nb_edges; ee++) {
    EntityHandle edge;
    rval = moab.side_element(tri, 1, ee, edge);
    MOAB_THROW(rval);
    rval = moab.side_number(tri, edge, side_number, sense, offset);
    MOAB_THROW(rval);
    const_cast<SideNumber_multiIndex &>(side_number_table)
        .insert(boost::shared_ptr<SideNumber>(
            new SideNumber(edge, ee, sense, offset)));
  }
  const_cast<SideNumber_multiIndex &>(side_number_table)
      .insert(boost::shared_ptr<SideNumber>(new SideNumber(tri, 0, 0, 0)));
}
const boost::shared_ptr<SideNumber> &
RefElementFace::getSideNumberPtr(const EntityHandle ent) const {
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  if (miit != side_number_table.end())
    return *miit;
  if (sPtr->ent == ent) {
    miit =
        const_cast<SideNumber_multiIndex &>(side_number_table)
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, 0, 0, 0)))
            .first;
    return *miit;
  }
  if (moab.type_from_handle(ent) == MBENTITYSET) {
    miit =
        const_cast<SideNumber_multiIndex &>(side_number_table)
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, 0, 0, 0)))
            .first;
    return *miit;
  }
  int side_number, sense, offset;
  rval = moab.side_number(sPtr->ent, ent, side_number, sense, offset);
  MOAB_THROW(rval);
  miit = const_cast<SideNumber_multiIndex &>(side_number_table)
             .insert(boost::shared_ptr<SideNumber>(
                 new SideNumber(ent, side_number, sense, offset)))
             .first;
  // std::cerr << side_number << " " << sense << " " << offset << std::endl;
  return *miit;
}
std::ostream &operator<<(std::ostream &os, const RefElementFace &e) {
  os << *e.sPtr;
  return os;
}
RefElement_EDGE::RefElement_EDGE(
    const boost::shared_ptr<RefEntity> &ref_ents_ptr)
    : RefElement(ref_ents_ptr) {
  switch (ref_ents_ptr->getEntType()) {
  case MBEDGE:
    break;
  default:
    THROW_MESSAGE("this work only for TRIs");
  }
}
const boost::shared_ptr<SideNumber> &
RefElement_EDGE::getSideNumberPtr(const EntityHandle ent) const {
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  if (miit != side_number_table.end())
    return *miit;
  if (sPtr->ent == ent) {
    miit =
        const_cast<SideNumber_multiIndex &>(side_number_table)
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, 0, 0, 0)))
            .first;
    return *miit;
  }
  if (moab.type_from_handle(ent) == MBENTITYSET) {
    miit =
        const_cast<SideNumber_multiIndex &>(side_number_table)
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, 0, 0, 0)))
            .first;
    return *miit;
  }
  int side_number, sense, offset;
  rval = moab.side_number(sPtr->ent, ent, side_number, sense, offset);
  MOAB_THROW(rval);
  miit = const_cast<SideNumber_multiIndex &>(side_number_table)
             .insert(boost::shared_ptr<SideNumber>(
                 new SideNumber(ent, side_number, sense, offset)))
             .first;
  // std::cerr << side_number << " " << sense << " " << offset << std::endl;
  return *miit;
}
std::ostream &operator<<(std::ostream &os, const RefElement_EDGE &e) {
  os << *e.sPtr;
  return os;
}
RefElement_VERTEX::RefElement_VERTEX(
    const boost::shared_ptr<RefEntity> &ref_ents_ptr)
    : RefElement(ref_ents_ptr) {
  switch (ref_ents_ptr->getEntType()) {
  case MBVERTEX:
    break;
  default:
    THROW_MESSAGE("this works only for TRIs");
  }
}
const boost::shared_ptr<SideNumber> &
RefElement_VERTEX::getSideNumberPtr(const EntityHandle ent) const {
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  SideNumber_multiIndex::iterator miit = side_number_table.find(ent);
  if (miit != side_number_table.end())
    return *miit;
  if (sPtr->ent == ent) {
    miit =
        const_cast<SideNumber_multiIndex &>(side_number_table)
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, 0, 0, 0)))
            .first;
    return *miit;
  }
  if (moab.type_from_handle(ent) == MBENTITYSET) {
    miit =
        const_cast<SideNumber_multiIndex &>(side_number_table)
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, 0, 0, 0)))
            .first;
    return *miit;
  }
  THROW_MESSAGE("no side entity for vertex if its is not an vertex itself");
  return nullSideNumber;
}
std::ostream &operator<<(std::ostream &os, const RefElement_VERTEX &e) {
  os << *e.sPtr;
  return os;
}

MoFEMErrorCode DefaultElementAdjacency::defaultVertex(
    moab::Interface &moab, const Field &field, const EntFiniteElement &fe,
    Range &adjacency) {
  MoFEMFunctionBegin;
  switch (field.getSpace()) {
  case H1:
    adjacency.insert(fe.getEnt());
    break;
  case NOFIELD: {
    Range ents;
    CHKERR moab.get_entities_by_handle(field.getMeshset(), ents, false);
    adjacency.merge(ents);
    for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(
              boost::shared_ptr<SideNumber>(new SideNumber(*eit, -1, 0, 0)));
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "this field is not implemented for VERTEX finite element");
  }
  // build side table
  for (Range::iterator eit = adjacency.begin(); eit != adjacency.end(); eit++) {
    fe.getSideNumberPtr(*eit);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DefaultElementAdjacency::defaultEdge(moab::Interface &moab,
                                                    const Field &field,
                                                    const EntFiniteElement &fe,
                                                    Range &adjacency) {
  MoFEMFunctionBegin;
  EntityHandle fe_ent = fe.getEnt();
  // Range nodes;
  switch (field.getSpace()) {
  case H1:
    CHKERR moab.get_connectivity(&fe_ent, 1, adjacency, true);
  case L2:
  case HCURL:
    adjacency.insert(fe_ent);
    break;
  case NOFIELD: {
    Range ents;
    CHKERR moab.get_entities_by_handle(field.getMeshset(), ents, false);
    adjacency.merge(ents);
    for (auto e : ents) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "this field is not implemented for EDGE finite element");
  }
  // build side table
  for (auto e : adjacency)
    fe.getSideNumberPtr(e);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DefaultElementAdjacency::defaultFace(moab::Interface &moab,
                                                    const Field &field,
                                                    const EntFiniteElement &fe,
                                                    Range &adjacency) {
  MoFEMFunctionBegin;
  // Range nodes,edges;
  const EntityHandle fe_ent = fe.getEnt();
  switch (field.getSpace()) {
  case H1:
    CHKERR moab.get_connectivity(&fe_ent, 1, adjacency, true);
  case HCURL:
    CHKERR moab.get_adjacencies(&fe_ent, 1, 1, false, adjacency,
                                moab::Interface::UNION);
  case HDIV:
    adjacency.insert(fe_ent);
    break;
  case NOFIELD: {
    Range ents;
    CHKERR moab.get_entities_by_handle(field.getMeshset(), ents, false);
    adjacency.merge(ents);
    for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(
              boost::shared_ptr<SideNumber>(new SideNumber(*eit, -1, 0, 0)));
    }
  } break;
  case L2:
    adjacency.insert(fe_ent); // add this just in case, if L2 is on skeleton
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "this field is not implemented for TRI finite element");
  }
  // build side table
  for (Range::iterator eit = adjacency.begin(); eit != adjacency.end(); eit++)
    fe.getSideNumberPtr(*eit);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DefaultElementAdjacency::defaultTet(moab::Interface &moab,
                                                   const Field &field,
                                                   const EntFiniteElement &fe,
                                                   Range &adjacency) {
  MoFEMFunctionBegin;
  EntityHandle fe_ent = fe.getEnt();
  switch (field.getSpace()) {
  case H1:
    CHKERR moab.get_connectivity(&fe_ent, 1, adjacency, true);
  case HCURL:
    CHKERR moab.get_adjacencies(&fe_ent, 1, 1, false, adjacency,
                                moab::Interface::UNION);
  case HDIV:
    CHKERR moab.get_adjacencies(&fe_ent, 1, 2, false, adjacency,
                                moab::Interface::UNION);
  case L2:
    adjacency.insert(fe_ent);
    break;
  case NOFIELD: {
    Range ents;
    CHKERR moab.get_entities_by_handle(field.getMeshset(), ents, false);
    adjacency.merge(ents);
    for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(
              boost::shared_ptr<SideNumber>(new SideNumber(*eit, -1, 0, 0)));
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "this field is not implemented for TRI finite element");
  }
  // build side table
  for (Range::iterator eit = adjacency.begin(); eit != adjacency.end(); eit++)
    fe.getSideNumberPtr(*eit);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DefaultElementAdjacency::defaultPrism(moab::Interface &moab,
                                                     const Field &field,
                                                     const EntFiniteElement &fe,
                                                     Range &adjacency) {
  MoFEMFunctionBegin;
  const EntityHandle prism = fe.getEnt();
  Range nodes;
  // initialize side sets
  fe.getRefElement()->getSideNumberPtr(prism);
  EntityHandle face_side3, face_side4;
  CHKERR moab.side_element(prism, 2, 3, face_side3);
  CHKERR moab.side_element(prism, 2, 4, face_side4);
  fe.getRefElement()->getSideNumberPtr(face_side3);
  fe.getRefElement()->getSideNumberPtr(face_side4);
  for (int qq = 0; qq < 3; qq++) {
    EntityHandle quad = 0;
    rval = moab.side_element(prism, 2, qq, quad);
    if (rval != MB_SUCCESS || quad == 0)
      continue;
    int side_number, sense, offset;
    rval = moab.side_number(prism, quad, side_number, sense, offset);
    if (side_number == -1 || rval != MB_SUCCESS)
      continue;
    const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
        .insert(boost::shared_ptr<SideNumber>(
            new SideNumber(quad, side_number, sense, offset)));
  }
  int ee = 0;
  for (; ee < 3; ee++) {
    EntityHandle edge = 0;
    CHKERR moab.side_element(prism, 1, ee, edge);
    boost::shared_ptr<SideNumber> side_ptr =
        fe.getRefElement()->getSideNumberPtr(edge);
    if (side_ptr->side_number != ee) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "data inconsistency for edge %d while in FE datastructure is "
               "numbered %d.",
               ee, side_ptr->side_number);
    }
    CHKERR moab.side_element(prism, 1, 6 + ee, edge);
    side_ptr = fe.getRefElement()->getSideNumberPtr(edge);
    if (side_ptr->side_number != ee + 6) {
      if (side_ptr->side_number != ee) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "data inconsistency for edge %d while in FE datastructure "
                 "is numbered %d.",
                 ee, side_ptr->side_number);
      } else {
        side_ptr->brother_side_number = ee + 6;
      }
    }
  }
  for (; ee < 6; ee++) {
    EntityHandle edge = 0;
    rval = moab.side_element(prism, 1, ee, edge);
    if (rval != MB_SUCCESS || edge == 0)
      continue;
    int side_number, sense, offset;
    rval = moab.side_number(prism, edge, side_number, sense, offset);
    if (side_number == -1 || rval != MB_SUCCESS)
      continue;
    const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
        .insert(boost::shared_ptr<SideNumber>(
            new SideNumber(edge, side_number, sense, offset)));
  }
  int nn = 0;
  for (; nn < 3; nn++) {
    EntityHandle node;
    CHKERR moab.side_element(prism, 0, nn, node);
    boost::shared_ptr<SideNumber> side_ptr =
        fe.getRefElement()->getSideNumberPtr(node);
    if (side_ptr->side_number != nn) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "data inconsistency for node %d while in FE datastructure is "
               "numbered %d.",
               nn, side_ptr->side_number);
    }
    CHKERR moab.side_element(prism, 0, nn + 3, node);
    side_ptr = fe.getRefElement()->getSideNumberPtr(node);
    if (side_ptr->side_number != nn + 3) {
      if (side_ptr->side_number != nn) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "data inconsistency for node %d while in FE datastructure is "
                 "numbered %d.",
                 nn, side_ptr->side_number);
      } else {
        side_ptr->brother_side_number = nn + 3;
      }
    }
  }

  // get adjacencies
  SideNumber_multiIndex &side_table = fe.getRefElement()->getSideNumberTable();
  switch (field.getSpace()) {
  case H1:
    // moab.get_connectivity(&prism,1,nodes,true);
    // use get adjacencies, this will allow take in account adjacencies set user
    CHKERR moab.get_adjacencies(&prism, 1, 0, false, nodes,
                                moab::Interface::UNION);
    {
      Range topo_nodes;
      CHKERR moab.get_connectivity(&prism, 1, topo_nodes, true);
      Range mid_nodes;
      CHKERR moab.get_connectivity(&prism, 1, mid_nodes, false);
      mid_nodes = subtract(mid_nodes, topo_nodes);
      nodes = subtract(nodes, mid_nodes);
    }
    adjacency.insert(nodes.begin(), nodes.end());
  case HCURL: {
    SideNumber_multiIndex::nth_index<2>::type::iterator siit, hi_siit;
    siit = side_table.get<2>().lower_bound(MBEDGE);
    hi_siit = side_table.get<2>().upper_bound(MBEDGE);
    for (; siit != hi_siit; siit++)
      adjacency.insert(siit->get()->ent);
  }
  case HDIV: {
    SideNumber_multiIndex::nth_index<2>::type::iterator siit, hi_siit;
    siit = side_table.get<2>().lower_bound(MBTRI);
    hi_siit = side_table.get<2>().upper_bound(MBTRI);
    for (; siit != hi_siit; siit++)
      adjacency.insert(siit->get()->ent);
    siit = side_table.get<2>().lower_bound(MBQUAD);
    hi_siit = side_table.get<2>().upper_bound(MBQUAD);
    for (; siit != hi_siit; siit++)
      adjacency.insert(siit->get()->ent);
  }
  case L2:
    adjacency.insert(prism);
    break;
  case NOFIELD: {
    Range ents;
    CHKERR moab.get_entities_by_handle(field.getMeshset(), ents, false);
    adjacency.merge(ents);
    for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(
              boost::shared_ptr<SideNumber>(new SideNumber(*eit, -1, 0, 0)));
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "this field is not implemented for TRI finite element");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DefaultElementAdjacency::defaultMeshset(
    moab::Interface &moab, const Field &field, const EntFiniteElement &fe,
    Range &adjacency) {
  MoFEMFunctionBegin;
  EntityHandle fe_ent = fe.getEnt();
  // get all meshsets in finite element meshset
  Range ent_ents_meshset;
  CHKERR moab.get_entities_by_type(fe_ent, MBENTITYSET, ent_ents_meshset,
                                   false);
  // resolve recursively all ents in the meshset
  Range ent_ents;
  CHKERR moab.get_entities_by_handle(fe_ent, ent_ents, true);
  switch (field.getSpace()) {
  case H1:
    adjacency.merge(ent_ents.subset_by_type(MBVERTEX));
  case HCURL:
    adjacency.merge(ent_ents.subset_by_type(MBEDGE));
  case HDIV:
    adjacency.merge(ent_ents.subset_by_type(MBTRI));
  case L2:
    adjacency.merge(ent_ents.subset_by_type(MBTET));
    break;
  case NOFIELD: {
    Range ents;
    CHKERR moab.get_entities_by_handle(field.getMeshset(), ents, false);
    adjacency.merge(ents);
    for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(
              boost::shared_ptr<SideNumber>(new SideNumber(*eit, -1, 0, 0)));
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  MoFEMFunctionReturn(0);
}

// FiniteElement
FiniteElement::FiniteElement(moab::Interface &moab, const EntityHandle _meshset)
    : meshset(_meshset) {
  Tag th_FEId;
  rval = moab.tag_get_handle("_FEId", th_FEId);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_FEId, &meshset, 1, (const void **)&tagId);
  MOAB_THROW(rval);
  Tag th_FEName;
  rval = moab.tag_get_handle("_FEName", th_FEName);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_FEName, &meshset, 1, (const void **)&tagName,
                             &tagNameSize);
  MOAB_THROW(rval);
  Tag th_FEIdCol, th_FEIdRow, th_FEIdData;
  rval = moab.tag_get_handle("_FEIdCol", th_FEIdCol);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_FEIdCol, &meshset, 1,
                             (const void **)&tag_BitFieldId_col_data);
  MOAB_THROW(rval);
  rval = moab.tag_get_handle("_FEIdRow", th_FEIdRow);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_FEIdRow, &meshset, 1,
                             (const void **)&tag_BitFieldId_row_data);
  MOAB_THROW(rval);
  rval = moab.tag_get_handle("_FEIdData", th_FEIdData);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_FEIdData, &meshset, 1,
                             (const void **)&tag_BitFieldId_data);
  MOAB_THROW(rval);

  elementAdjacencyTable[MBVERTEX] = DefaultElementAdjacency::defaultVertex;
  elementAdjacencyTable[MBEDGE] = DefaultElementAdjacency::defaultEdge;
  elementAdjacencyTable[MBTRI] = DefaultElementAdjacency::defaultFace;
  elementAdjacencyTable[MBQUAD] = DefaultElementAdjacency::defaultFace;
  elementAdjacencyTable[MBTET] = DefaultElementAdjacency::defaultTet;
  elementAdjacencyTable[MBPRISM] = DefaultElementAdjacency::defaultPrism;
  elementAdjacencyTable[MBENTITYSET] = DefaultElementAdjacency::defaultMeshset;

  feUId = static_cast<UId>(getBitNumber()) << 8 * sizeof(EntityHandle);
}

std::ostream &operator<<(std::ostream &os, const FiniteElement &e) {
  os << e.getNameRef() << " fe_id " << e.getId().to_ulong() << " f_id_row "
     << e.getBitFieldIdRow() << " f_id_col " << e.getBitFieldIdCol()
     << " BitFEId_data " << e.getBitFieldIdData();
  return os;
}

void FiniteElement_col_change_bit_add::operator()(
    boost::shared_ptr<FiniteElement> &fe) {
  *static_cast<BitFieldId *>(fe->tag_BitFieldId_col_data) |= fIdCol;
}

void FiniteElement_row_change_bit_add::operator()(
    boost::shared_ptr<FiniteElement> &fe) {
  *static_cast<BitFieldId *>(fe->tag_BitFieldId_row_data) |= fIdRow;
}

void FiniteElement_change_bit_add::operator()(
    boost::shared_ptr<FiniteElement> &fe) {
  *static_cast<BitFieldId *>(fe->tag_BitFieldId_data) |= fIdData;
}

void FiniteElement_col_change_bit_off::operator()(
    boost::shared_ptr<FiniteElement> &fe) {
  *static_cast<BitFieldId *>(fe->tag_BitFieldId_col_data) &= fIdCol.flip();
}

void FiniteElement_row_change_bit_off::operator()(
    boost::shared_ptr<FiniteElement> &fe) {
  *static_cast<BitFieldId *>(fe->tag_BitFieldId_row_data) &= fIdRow.flip();
}

void FiniteElement_change_bit_off::operator()(
    boost::shared_ptr<FiniteElement> &fe) {
  *static_cast<BitFieldId *>(fe->tag_BitFieldId_data) &= fIdData.flip();
}

void FiniteElement_col_change_bit_reset::operator()(
    boost::shared_ptr<FiniteElement> &fe) {
  static_cast<BitFieldId *>(fe->tag_BitFieldId_col_data)->reset();
}

void FiniteElement_row_change_bit_reset::operator()(
    boost::shared_ptr<FiniteElement> &fe) {
  static_cast<BitFieldId *>(fe->tag_BitFieldId_row_data)->reset();
}

void FiniteElement_change_bit_reset::operator()(
    boost::shared_ptr<FiniteElement> &fe) {
  static_cast<BitFieldId *>(fe->tag_BitFieldId_data)->reset();
}

// FiniteElement data
EntFiniteElement::EntFiniteElement(
    const boost::shared_ptr<RefElement> &ref_finite_element,
    const boost::shared_ptr<FiniteElement> &fe_ptr)
    : interface_FiniteElement<FiniteElement>(fe_ptr),
      interface_RefElement<RefElement>(ref_finite_element),
      dataDofs(new FEDofEntity_multiIndex()),
      dataFieldEnts(new FieldEntity_multiIndex_spaceType_view()),
      rowFieldEnts(new FieldEntity_vector_view()),
      colFieldEnts(new FieldEntity_vector_view()) {}

std::ostream &operator<<(std::ostream &os, const EntFiniteElement &e) {
  os << *e.sFePtr << std::endl;
  os << *e.sPtr << std::endl;
  os << "data dof_uids ";
  for (auto &dit : *e.dataDofs) {
    if (!dit) {
      os << "null ptr";
    } else {
      if (!dit->getDofEntityPtr()) {
        os << "( null ptr to dof ) ";
      } else {
        if (!dit->getFieldEntityPtr()) {
          os << "(( null ptr to field entity )) ";
        } else {
          os << dit->getGlobalUniqueId() << " ";
        }
      }
    }
  }
  return os;
}

MoFEMErrorCode
EntFiniteElement::getElementAdjacency(const boost::shared_ptr<Field> field_ptr,
                                      Range &adjacency) {
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  MoFEMFunctionBegin;
  const EntFiniteElement *this_fe_ptr = this;
  if (get_MoFEMFiniteElementPtr()->elementAdjacencyTable[getEntType()] ==
      NULL) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  CHKERR get_MoFEMFiniteElementPtr()->elementAdjacencyTable[getEntType()](
      moab, *field_ptr, *this_fe_ptr, adjacency);
  MoFEMFunctionReturn(0);
}

/**
 * \Construct indexed finite element
 */
NumeredEntFiniteElement::NumeredEntFiniteElement(
    const boost::shared_ptr<EntFiniteElement> &sptr)
    : interface_EntFiniteElement<EntFiniteElement>(sptr), part(-1),
      rowDofs(boost::shared_ptr<FENumeredDofEntity_multiIndex>(
          new FENumeredDofEntity_multiIndex())),
      colDofs(boost::shared_ptr<FENumeredDofEntity_multiIndex>(
          new FENumeredDofEntity_multiIndex())){};

boost::weak_ptr<FENumeredDofEntity>
NumeredEntFiniteElement::getRowDofsByPetscGlobalDofIdx(const int idx) const {
  auto comp = [idx](const auto &a) { return a->getPetscGlobalDofIdx() == idx; };
  auto dit = std::find_if(rowDofs->begin(), rowDofs->end(), comp);
  if (dit != rowDofs->end())
    return *dit;
  else
    return boost::weak_ptr<FENumeredDofEntity>();
}

boost::weak_ptr<FENumeredDofEntity>
NumeredEntFiniteElement::getColDofsByPetscGlobalDofIdx(const int idx) const {
  auto comp = [idx](const auto &a) { return a->getPetscGlobalDofIdx() == idx; };
  auto dit = std::find_if(colDofs->begin(), colDofs->end(), comp);
  if (dit != colDofs->end())
    return *dit;
  else
    return boost::weak_ptr<FENumeredDofEntity>();
}

} // namespace MoFEM
