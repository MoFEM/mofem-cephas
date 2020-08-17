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
    auto siit = side_table.get<0>().lower_bound(get_id_for_min_type<MBEDGE>());
    auto hi_siit =
        side_table.get<0>().upper_bound(get_id_for_max_type<MBEDGE>());
    for (; siit != hi_siit; siit++)
      adjacency.insert(siit->get()->ent);
  }
  case HDIV: {
    auto siit = side_table.get<0>().lower_bound(get_id_for_min_type<MBTRI>());
    auto hi_siit =
        side_table.get<0>().upper_bound(get_id_for_max_type<MBQUAD>());
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

const EntFiniteElement *EntFiniteElement::lastSeenDataEntFiniteElement =
    nullptr;
boost::shared_ptr<FEDofEntity_multiIndex> EntFiniteElement::dataDofs;

const EntFiniteElement *EntFiniteElement::lastSeenDataVectorEntFiniteElement =
    nullptr;
boost::shared_ptr<std::vector<boost::shared_ptr<FEDofEntity>>>
    EntFiniteElement::dataVectorDofs;

// FiniteElement data
EntFiniteElement::EntFiniteElement(
    const boost::shared_ptr<RefElement> &ref_finite_element,
    const boost::shared_ptr<FiniteElement> &fe_ptr)
    : interface_FiniteElementImpl<FiniteElement, RefElement>(
          fe_ptr, ref_finite_element),
      finiteElementPtr(fe_ptr), dataFieldEnts(new FieldEntity_vector_view()),
      rowFieldEnts(new FieldEntity_vector_view()),
      colFieldEnts(new FieldEntity_vector_view()) {}

EntFiniteElement::~EntFiniteElement() {
  dataVectorDofs.reset();
  lastSeenDataVectorEntFiniteElement = nullptr;
  dataVectorDofs.reset();
  lastSeenDataEntFiniteElement = nullptr;
}

std::ostream &operator<<(std::ostream &os, const EntFiniteElement &e) {
  os << *e.getFiniteElementPtr() << std::endl;
  os << *e.sPtr << std::endl;
  return os;
}

MoFEMErrorCode
EntFiniteElement::getElementAdjacency(const boost::shared_ptr<Field> field_ptr,
                                      Range &adjacency) {
  moab::Interface &moab = getRefEntityPtr()->getBasicDataPtr()->moab;
  MoFEMFunctionBegin;
  const EntFiniteElement *this_fe_ptr = this;
  if (getFiniteElementPtr()->elementAdjacencyTable[getEntType()] == NULL) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }
  CHKERR getFiniteElementPtr()->elementAdjacencyTable[getEntType()](
      moab, *field_ptr, *this_fe_ptr, adjacency);
  MoFEMFunctionReturn(0);
}

const NumeredEntFiniteElement
    *NumeredEntFiniteElement::lastSeenRowFiniteElement = nullptr;
boost::shared_ptr<FENumeredDofEntity_multiIndex>
    NumeredEntFiniteElement::rowDofs; ///< indexed dofs on rows

const NumeredEntFiniteElement
    *NumeredEntFiniteElement::lastSeenColFiniteElement = nullptr;
boost::shared_ptr<FENumeredDofEntity_multiIndex>
    NumeredEntFiniteElement::colDofs; ///< indexed dofs on columns

/**
 * \Construct indexed finite element
 */
NumeredEntFiniteElement::NumeredEntFiniteElement(
    const boost::shared_ptr<EntFiniteElement> &sptr)
    : interface_EntFiniteElement<EntFiniteElement>(sptr), part(-1){};

NumeredEntFiniteElement::~NumeredEntFiniteElement() {
  lastSeenRowFiniteElement = nullptr;
  rowDofs.reset();
  lastSeenColFiniteElement = nullptr;
  colDofs.reset();
}

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

boost::shared_ptr<FEDofEntity_multiIndex> &
EntFiniteElement::getDataDofsPtr() const {
  RefEntityTmp<0>::refElementPtr = this->getRefElement();
  if (lastSeenDataEntFiniteElement != this) {
    if (dataDofs)
      dataDofs->clear();
    else
      dataDofs = boost::make_shared<FEDofEntity_multiIndex>();

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

    if (get_cache_data_dofs_view(dataFieldEnts, dataDofs, Extractor(),
                                 Inserter()))
      THROW_MESSAGE("dataDofs can not be created");
    lastSeenDataEntFiniteElement = this;
  }
  return dataDofs;
};

boost::shared_ptr<std::vector<boost::shared_ptr<FEDofEntity>>> &
EntFiniteElement::getDataVectorDofsPtr() const {
  RefEntityTmp<0>::refElementPtr = this->getRefElement();
  if (lastSeenDataVectorEntFiniteElement != this) {
    if (dataVectorDofs)
      dataVectorDofs->clear();
    else
      dataVectorDofs =
          boost::make_shared<std::vector<boost::shared_ptr<FEDofEntity>>>();

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

    if (get_cache_data_dofs_view(dataFieldEnts, dataVectorDofs, Extractor(),
                                 Inserter()))
      THROW_MESSAGE("dataDofs can not be created");
    lastSeenDataVectorEntFiniteElement = this;
  }
  return dataVectorDofs;
};

boost::shared_ptr<FENumeredDofEntity_multiIndex> &
NumeredEntFiniteElement::getRowDofsPtr() const {
  RefEntityTmp<0>::refElementPtr = this->getRefElement();
  if (lastSeenRowFiniteElement != this) {
    if (rowDofs)
      rowDofs->clear();
    else
      rowDofs = boost::make_shared<FENumeredDofEntity_multiIndex>();

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

    if (get_cache_data_dofs_view(getRowFieldEntsPtr(), rowDofs, Extractor(),
                                 Inserter()))
      THROW_MESSAGE("rowDofs can not be created");

    lastSeenRowFiniteElement = this;
  }
  return rowDofs;
}

boost::shared_ptr<FENumeredDofEntity_multiIndex> &
NumeredEntFiniteElement::getColDofsPtr() const {
  RefEntityTmp<0>::refElementPtr = this->getRefElement();
  if (lastSeenColFiniteElement != this) {

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

    if (colDofs)
      colDofs->clear();
    else
      colDofs = boost::make_shared<FENumeredDofEntity_multiIndex>();

    if (get_cache_data_dofs_view(getColFieldEntsPtr(), colDofs, Extractor(),
                                 Inserter()))
      THROW_MESSAGE("colDofs can not be created");

    lastSeenColFiniteElement = this;
  }
  return colDofs;
}

} // namespace MoFEM
