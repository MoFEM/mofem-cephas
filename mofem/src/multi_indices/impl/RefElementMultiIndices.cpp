/** \file RefElementMultiIndices.cpp
 * \brief Multi-index containers for finite elements
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
  EntityHandle prism = getEnt();
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
  if (type_from_handle(ent) == MBENTITYSET) {
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

    if (type_from_handle(ent) == MBVERTEX) {
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
    sss << "this not working: " << ent << " type: " << type_from_handle(ent)
        << " " << MBEDGE << " " << MBTRI << std::endl;
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

RefElementVolume::RefElementVolume(const boost::shared_ptr<RefEntity> &ref_ents_ptr)
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
  case MBHEX:
    break;
  default:
    THROW_MESSAGE("this work only for TETs or HEXs");
  }
  const_cast<SideNumber_multiIndex &>(side_number_table)
      .insert(boost::make_shared<SideNumber>(sPtr->ent, 0, 0, 0));
}

const boost::shared_ptr<SideNumber> &
RefElementVolume::getSideNumberPtr(const EntityHandle ent) const {
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
  if (type_from_handle(ent) == MBENTITYSET) {
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
std::ostream &operator<<(std::ostream &os, const RefElementVolume &e) {
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
  EntityHandle tri = getEnt();
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
  if (type_from_handle(ent) == MBENTITYSET) {
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
  if (type_from_handle(ent) == MBENTITYSET) {
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
  if (type_from_handle(ent) == MBENTITYSET) {
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

} // namespace MoFEM
