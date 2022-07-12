/** \file FEMultiIndices.cpp
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

constexpr DefaultElementAdjacency::DefEntTypeMap
    DefaultElementAdjacency::defVertexTypeMap;
constexpr DefaultElementAdjacency::DefEntTypeMap
    DefaultElementAdjacency::defEdgeTypeMap;
constexpr DefaultElementAdjacency::DefEntTypeMap
    DefaultElementAdjacency::defTriTypeMap;
constexpr DefaultElementAdjacency::DefEntTypeMap
    DefaultElementAdjacency::defQuadTypeMap;
constexpr DefaultElementAdjacency::DefEntTypeMap
    DefaultElementAdjacency::defTetTypeMap;
constexpr DefaultElementAdjacency::DefEntTypeMap
    DefaultElementAdjacency::defHexTypeMap;
constexpr DefaultElementAdjacency::DefEntTypeMap
    DefaultElementAdjacency::defPrismTypeMap;
constexpr DefaultElementAdjacency::DefEntTypeMap
    DefaultElementAdjacency::defMeshsetTypeMap;
constexpr std::array<const DefaultElementAdjacency::DefEntTypeMap *, MBMAXTYPE>
    DefaultElementAdjacency::defTypeMap;

MoFEMErrorCode DefaultElementAdjacency::defaultVertex(
    moab::Interface &moab, const Field &field, const EntFiniteElement &fe,
    std::vector<EntityHandle> &adjacency) {
  MoFEMFunctionBegin;
  switch (field.getSpace()) {
  case H1:
    adjacency.push_back(fe.getEnt());
    // build side table
    for (auto ent : adjacency)
      fe.getSideNumberPtr(ent);
    break;
  case NOFIELD: {
    CHKERR moab.get_entities_by_handle(field.getMeshset(), adjacency, false);
    for (auto ent : adjacency) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, -1, 0, 0)));
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "this field is not implemented for VERTEX finite element");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
DefaultElementAdjacency::defaultEdge(moab::Interface &moab, const Field &field,
                                     const EntFiniteElement &fe,
                                     std::vector<EntityHandle> &adjacency) {
  MoFEMFunctionBegin;
  EntityHandle fe_ent = fe.getEnt();
  // Range nodes;
  switch (field.getSpace()) {
  case H1:
    CHKERR moab.get_connectivity(&fe_ent, 1, adjacency, true);
  case L2:
  case HCURL:
    adjacency.push_back(fe_ent);
    // build side table
    for (auto e : adjacency)
      fe.getSideNumberPtr(e);
    break;
  case NOFIELD: {
    CHKERR moab.get_entities_by_handle(field.getMeshset(), adjacency, false);
    for (auto e : adjacency) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "this field is not implemented for EDGE finite element");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
DefaultElementAdjacency::defaultFace(moab::Interface &moab, const Field &field,
                                     const EntFiniteElement &fe,
                                     std::vector<EntityHandle> &adjacency) {
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
  case L2:
    adjacency.push_back(fe_ent);
    // build side table
    for (auto ent : adjacency)
      fe.getSideNumberPtr(ent);
    break;
  case NOFIELD: {
    CHKERR moab.get_entities_by_handle(field.getMeshset(), adjacency, false);
    for (auto ent : adjacency) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, -1, 0, 0)));
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "this field is not implemented for TRI finite element");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DefaultElementAdjacency::defaultVolume(
    moab::Interface &moab, const Field &field, const EntFiniteElement &fe,
    std::vector<EntityHandle> &adjacency) {
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
    adjacency.push_back(fe_ent);
    // build side table
    for (auto ent : adjacency)
      fe.getSideNumberPtr(ent);
    break;
  case NOFIELD: {
    CHKERR moab.get_entities_by_handle(field.getMeshset(), adjacency, false);
    for (auto ent : adjacency) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, -1, 0, 0)));
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "this field is not implemented for TRI finite element");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
DefaultElementAdjacency::defaultPrism(moab::Interface &moab, const Field &field,
                                      const EntFiniteElement &fe,
                                      std::vector<EntityHandle> &adjacency) {
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
    adjacency.insert(adjacency.end(), nodes.begin(), nodes.end());
  case HCURL: {
    auto siit = side_table.get<0>().lower_bound(get_id_for_min_type<MBEDGE>());
    auto hi_siit =
        side_table.get<0>().upper_bound(get_id_for_max_type<MBEDGE>());
    for (; siit != hi_siit; siit++)
      adjacency.push_back(siit->get()->ent);
  }
  case HDIV: {
    auto siit = side_table.get<0>().lower_bound(get_id_for_min_type<MBTRI>());
    auto hi_siit =
        side_table.get<0>().upper_bound(get_id_for_max_type<MBQUAD>());
    for (; siit != hi_siit; siit++)
      adjacency.push_back(siit->get()->ent);
  }
  case L2:
    adjacency.push_back(prism);
    break;
  case NOFIELD: {
    CHKERR moab.get_entities_by_handle(field.getMeshset(), adjacency, false);
    for (auto ent : adjacency) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, -1, 0, 0)));
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
    std::vector<EntityHandle> &adjacency) {
  MoFEMFunctionBegin;
  EntityHandle fe_ent = fe.getEnt();
  // resolve recursively all ents in the meshset
  switch (field.getSpace()) {
  case H1:
    CHKERR moab.get_entities_by_type(fe_ent, MBVERTEX, adjacency, true);
  case HCURL:
    CHKERR moab.get_entities_by_type(fe_ent, MBEDGE, adjacency, true);
  case HDIV:
    CHKERR moab.get_entities_by_dimension(fe_ent, 2, adjacency, true);
  case L2:
    CHKERR moab.get_entities_by_dimension(fe_ent, 3, adjacency, true);
    break;
  case NOFIELD: {
    CHKERR moab.get_entities_by_handle(field.getMeshset(), adjacency, false);
    for (auto ent : adjacency) {
      const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
          .insert(boost::shared_ptr<SideNumber>(new SideNumber(ent, -1, 0, 0)));
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
  elementAdjacencyTable[MBTET] = DefaultElementAdjacency::defaultVolume;
  elementAdjacencyTable[MBHEX] = DefaultElementAdjacency::defaultVolume;
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
    : interface_FiniteElement<FiniteElement, RefElement>(fe_ptr,
                                                         ref_finite_element),
      dataFieldEnts(new FieldEntity_vector_view()),
      rowFieldEnts(new FieldEntity_vector_view()),
      colFieldEnts(new FieldEntity_vector_view()) {}

std::ostream &operator<<(std::ostream &os, const EntFiniteElement &e) {
  os << *e.getFiniteElementPtr() << std::endl;
  os << *e.sPtr;
  return os;
}

MoFEMErrorCode
EntFiniteElement::getElementAdjacency(const boost::shared_ptr<Field> field_ptr,
                                      std::vector<EntityHandle> &adjacency) {
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  MoFEMFunctionBegin;
  const EntFiniteElement *this_fe_ptr = this;
  if (getFiniteElementPtr()->elementAdjacencyTable[getEntType()] == NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  CHKERR getFiniteElementPtr()->elementAdjacencyTable[getEntType()](
      moab, *field_ptr, *this_fe_ptr, adjacency);
  MoFEMFunctionReturn(0);
}

/**
 * \Construct indexed finite element
 */
NumeredEntFiniteElement::NumeredEntFiniteElement(
    const boost::shared_ptr<EntFiniteElement> &sptr)
    : interface_EntFiniteElement<EntFiniteElement>(sptr), part(-1){};

boost::weak_ptr<FENumeredDofEntity>
NumeredEntFiniteElement::getRowDofsByPetscGlobalDofIdx(const int idx) const {
  auto comp = [idx](const auto &a) { return a->getPetscGlobalDofIdx() == idx; };

  for (auto &it : getRowFieldEnts()) {
    if (auto e = it.lock()) {
      if (auto cache = e->entityCacheColDofs.lock()) {
        auto dit = std::find_if(cache->loHi[0], cache->loHi[1], comp);
        if (dit != cache->loHi[1])
          return boost::reinterpret_pointer_cast<FENumeredDofEntity>(*dit);
      } else
        THROW_MESSAGE("Cache not set");
    }
  }

  return boost::weak_ptr<FENumeredDofEntity>();
}

boost::weak_ptr<FENumeredDofEntity>
NumeredEntFiniteElement::getColDofsByPetscGlobalDofIdx(const int idx) const {

  auto comp = [idx](const auto &a) { return a->getPetscGlobalDofIdx() == idx; };

  for (auto &it : getColFieldEnts()) {
    if (auto e = it.lock()) {
      if (auto cache = e->entityCacheColDofs.lock()) {
        auto dit = std::find_if(cache->loHi[0], cache->loHi[1], comp);
        if (dit != cache->loHi[1])
          return boost::reinterpret_pointer_cast<FENumeredDofEntity>(*dit);
      } else
        THROW_MESSAGE("Cache not set");
    }
  }

  return boost::weak_ptr<FENumeredDofEntity>();
}

std::ostream &operator<<(std::ostream &os, const NumeredEntFiniteElement &e) {
  os << "part " << e.part << " " << *(e.getEntFiniteElement());
  return os;
}

template <typename ENTSVIEW, typename DOFSVIEW, typename EXTRACTOR,
          typename INSERTER>
inline static MoFEMErrorCode
get_cache_data_dofs_view(ENTSVIEW &ents_view, DOFSVIEW &dofs_view,
                         EXTRACTOR &&extractor, INSERTER &&inserter) {
  MoFEMFunctionBeginHot;

  auto hint = dofs_view->end();
  using ValType = typename std::remove_reference<decltype(**hint)>::type;

  for (auto &it : *ents_view) {
    if (auto e = it.lock()) {

      if (auto cache = extractor(e).lock())
        for (auto dit = cache->loHi[0]; dit != cache->loHi[1]; ++dit)
          hint = inserter(dofs_view, hint,
                          boost::reinterpret_pointer_cast<ValType>(*dit));
      else
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Cache not set");
    }
  }

  MoFEMFunctionReturnHot(0);
}

boost::shared_ptr<FEDofEntity_multiIndex>
EntFiniteElement::getDataDofsPtr() const {
  RefEntityTmp<0>::refElementPtr = this->getRefElement();
  struct Extractor {
    boost::weak_ptr<EntityCacheDofs>
    operator()(boost::shared_ptr<FieldEntity> &e) {
      return e->entityCacheDataDofs;
    }
  };

  struct Inserter {
    FEDofEntity_multiIndex::iterator
    operator()(boost::shared_ptr<FEDofEntity_multiIndex> &dofs_view,
               FEDofEntity_multiIndex::iterator &hint,
               boost::shared_ptr<FEDofEntity> &&dof) {
      return dofs_view->emplace_hint(hint, dof);
    }
  };

  auto data_dofs = boost::make_shared<FEDofEntity_multiIndex>();
  if (get_cache_data_dofs_view(dataFieldEnts, data_dofs, Extractor(),
                               Inserter()))
    THROW_MESSAGE("data_dofs can not be created");
  return data_dofs;
};

boost::shared_ptr<std::vector<boost::shared_ptr<FEDofEntity>>>
EntFiniteElement::getDataVectorDofsPtr() const {
  RefEntityTmp<0>::refElementPtr = this->getRefElement();

  struct Extractor {
    boost::weak_ptr<EntityCacheDofs>
    operator()(boost::shared_ptr<FieldEntity> &e) {
      return e->entityCacheDataDofs;
    }
  };

  struct Inserter {
    using Vec = std::vector<boost::shared_ptr<FEDofEntity>>;
    using It = Vec::iterator;
    It operator()(boost::shared_ptr<Vec> &dofs_view, It &hint,
                  boost::shared_ptr<FEDofEntity> &&dof) {
      dofs_view->emplace_back(dof);
      return dofs_view->end();
    }
  };

  auto data_vector_dofs =
      boost::make_shared<std::vector<boost::shared_ptr<FEDofEntity>>>();
  if (get_cache_data_dofs_view(dataFieldEnts, data_vector_dofs, Extractor(),
                               Inserter()))
    THROW_MESSAGE("dataDofs can not be created");

  return data_vector_dofs;
};

boost::shared_ptr<FENumeredDofEntity_multiIndex>
NumeredEntFiniteElement::getRowDofsPtr() const {
  RefEntityTmp<0>::refElementPtr = this->getRefElement();

  struct Extractor {
    boost::weak_ptr<EntityCacheNumeredDofs>
    operator()(boost::shared_ptr<FieldEntity> &e) {
      return e->entityCacheRowDofs;
    }
  };

  struct Inserter {
    using Idx = FENumeredDofEntity_multiIndex;
    using It = Idx::iterator;
    It operator()(boost::shared_ptr<Idx> &dofs_view, It &hint,
                  boost::shared_ptr<FENumeredDofEntity> &&dof) {
      return dofs_view->emplace_hint(hint, dof);
    }
  };

  auto row_dofs = boost::make_shared<FENumeredDofEntity_multiIndex>();
  if (get_cache_data_dofs_view(getRowFieldEntsPtr(), row_dofs, Extractor(),
                               Inserter()))
    THROW_MESSAGE("row_dofs can not be created");

  return row_dofs;
}

boost::shared_ptr<FENumeredDofEntity_multiIndex>
NumeredEntFiniteElement::getColDofsPtr() const {
  RefEntityTmp<0>::refElementPtr = this->getRefElement();

  struct Extractor {
    boost::weak_ptr<EntityCacheNumeredDofs>
    operator()(boost::shared_ptr<FieldEntity> &e) {
      return e->entityCacheColDofs;
    }
  };

  struct Inserter {
    using Idx = FENumeredDofEntity_multiIndex;
    using It = Idx::iterator;
    It operator()(boost::shared_ptr<Idx> &dofs_view, It &hint,
                  boost::shared_ptr<FENumeredDofEntity> &&dof) {
      return dofs_view->emplace_hint(hint, dof);
    }
  };

  auto col_dofs = boost::make_shared<FENumeredDofEntity_multiIndex>();
  if (get_cache_data_dofs_view(getColFieldEntsPtr(), col_dofs, Extractor(),
                               Inserter()))
    THROW_MESSAGE("col_dofs can not be created");

  return col_dofs;
}

} // namespace MoFEM
