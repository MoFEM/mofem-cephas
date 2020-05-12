/** \file FECore.cpp
 * \brief Core interface methods for managing deletions and insertion dofs
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

MoFEMErrorCode
Core::get_finite_elements(const FiniteElement_multiIndex **fe_ptr) const {
  MoFEMFunctionBeginHot;
  *fe_ptr = &finiteElements;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::get_ents_finite_elements(
    const EntFiniteElement_multiIndex **fe_ent_ptr) const {
  MoFEMFunctionBeginHot;
  *fe_ent_ptr = &entsFiniteElements;
  MoFEMFunctionReturnHot(0);
}

bool Core::check_finite_element(const std::string &name) const {
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FeSetByName;
  const FeSetByName &set = finiteElements.get<FiniteElement_name_mi_tag>();
  FeSetByName::iterator miit = set.find(name);
  if (miit == set.end())
    return false;
  return true;
}

MoFEMErrorCode Core::add_finite_element(const std::string &fe_name,
                                        enum MoFEMTypes bh, int verb) {
  MoFEMFunctionBegin;
  *buildMoFEM &= 1 << 0;
  if (verb == -1) {
    verb = verbose;
  }
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe =
      finite_element_name_set.find(fe_name);
  if (bh == MF_EXCL) {
    if (it_fe != finite_element_name_set.end()) {
      SETERRQ1(cOmm, MOFEM_NOT_FOUND, "this < %s > is there", fe_name.c_str());
    }
  } else {
    if (it_fe != finite_element_name_set.end())
      MoFEMFunctionReturnHot(0);
  }
  EntityHandle meshset;
  CHKERR get_moab().create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, meshset);

  // id
  BitFEId id = getFEShift();
  CHKERR get_moab().tag_set_data(th_FEId, &meshset, 1, &id);

  // Add finite element meshset to partion meshset. In case of no elements
  // on processor part, when mesh file is read, finite element meshset is
  // prevented from deletion by moab reader.
  auto add_meshset_to_partition = [&](auto meshset) {
    MoFEMFunctionBegin;
    const void *tag_vals[] = {&rAnk};
    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &get_moab(), get_basic_entity_data_ptr()->pcommID);
    Tag part_tag = pcomm->part_tag();
    Range tagged_sets;
    CHKERR get_moab().get_entities_by_type_and_tag(0, MBENTITYSET, &part_tag,
                                                   tag_vals, 1, tagged_sets,
                                                   moab::Interface::UNION);
    for (auto s : tagged_sets)
      CHKERR get_moab().add_entities(s, &meshset, 1);
    MoFEMFunctionReturn(0);
  };
  CHKERR add_meshset_to_partition(meshset);

  // id name
  void const *tag_data[] = {fe_name.c_str()};
  int tag_sizes[1];
  tag_sizes[0] = fe_name.size();
  CHKERR get_moab().tag_set_by_ptr(th_FEName, &meshset, 1, tag_data, tag_sizes);

  // add FiniteElement
  std::pair<FiniteElement_multiIndex::iterator, bool> p = finiteElements.insert(
      boost::shared_ptr<FiniteElement>(new FiniteElement(moab, meshset)));
  if (!p.second)
    SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "FiniteElement not inserted");
  if (verb > 0) {
    std::ostringstream ss;
    ss << "add finite element: " << fe_name << std::endl;
    PetscPrintf(cOmm, ss.str().c_str());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_finite_element_adjacency_table(const std::string &fe_name,
                                            const EntityType type,
                                            ElementAdjacencyFunct function) {
  MoFEMFunctionBeginHot;
  *buildMoFEM &= 1 << 0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe =
      finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ(cOmm, MOFEM_NOT_FOUND,
            "This finite element is not defined (advise: check spelling)");
  boost::shared_ptr<FiniteElement> fe;
  fe = *it_fe;
  fe->elementAdjacencyTable[type] = function;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_finite_element_add_field_data(const std::string &fe_name,
                                           const std::string &name_data) {
  MoFEMFunctionBegin;
  *buildMoFEM &= 1 << 0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe =
      finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ(cOmm, MOFEM_NOT_FOUND,
            "This finite element is not defined (advise: check spelling)");
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_change_bit_add(getBitFieldId(name_data)));
  if (!success)
    SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_finite_element_add_field_row(const std::string &fe_name,
                                          const std::string &name_row) {
  MoFEMFunctionBegin;
  *buildMoFEM &= 1 << 0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe =
      finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ1(cOmm, MOFEM_NOT_FOUND, "this < %s > is not there",
             fe_name.c_str());
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_row_change_bit_add(getBitFieldId(name_row)));
  if (!success)
    SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_finite_element_add_field_col(const std::string &fe_name,
                                          const std::string &name_col) {
  MoFEMFunctionBegin;
  *buildMoFEM &= 1 << 0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe =
      finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "this FiniteElement is there");
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_col_change_bit_add(getBitFieldId(name_col)));
  if (!success)
    SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_finite_element_off_field_data(const std::string &fe_name,
                                           const std::string &name_data) {
  MoFEMFunctionBegin;
  *buildMoFEM &= 1 << 0;
  auto &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  auto it_fe = finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ(cOmm, MOFEM_NOT_FOUND, "this FiniteElement is there");
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_change_bit_off(getBitFieldId(name_data)));
  if (!success)
    SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_finite_element_off_field_row(const std::string &fe_name,
                                          const std::string &name_row) {
  MoFEMFunctionBegin;
  *buildMoFEM &= 1 << 0;
  auto &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  auto it_fe = finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ1(cOmm, MOFEM_NOT_FOUND, "this < %s > is not there",
             fe_name.c_str());
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_row_change_bit_off(getBitFieldId(name_row)));
  if (!success)
    SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_finite_element_off_field_col(const std::string &fe_name,
                                          const std::string &name_col) {
  MoFEMFunctionBegin;
  *buildMoFEM &= 1 << 0;
  auto &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  auto it_fe = finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ(cOmm, MOFEM_NOT_FOUND, "this FiniteElement is there");
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_col_change_bit_off(getBitFieldId(name_col)));
  if (!success)
    SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

BitFEId Core::getBitFEId(const std::string &name) const {
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  const FiniteElements_by_name &set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator miit = set.find(name);
  if (miit == set.end())
    THROW_MESSAGE(
        ("finite element < " + name + " > not found (top tip: check spelling)")
            .c_str());
  return (*miit)->getId();
}

std::string Core::getBitFEIdName(const BitFEId id) const {
  typedef FiniteElement_multiIndex::index<BitFEId_mi_tag>::type
      finiteElements_by_id;
  const finiteElements_by_id &set = finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = set.find(id);
  assert(miit != set.end());
  return (*miit)->getName();
}

EntityHandle Core::get_finite_element_meshset(const BitFEId id) const {
  typedef FiniteElement_multiIndex::index<BitFEId_mi_tag>::type
      finiteElements_by_id;
  const finiteElements_by_id &set = finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = set.find(id);
  if (miit == set.end())
    THROW_MESSAGE("finite element not found");
  return (*miit)->meshset;
}

EntityHandle Core::get_finite_element_meshset(const std::string &name) const {
  return get_finite_element_meshset(getBitFEId(name));
}

MoFEMErrorCode
Core::get_finite_element_entities_by_dimension(const std::string name, int dim,
                                               Range &ents) const {

  MoFEMFunctionBegin;

  EntityHandle meshset = get_finite_element_meshset(name);
  CHKERR get_moab().get_entities_by_dimension(meshset, dim, ents, true);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::get_finite_element_entities_by_type(const std::string name,
                                                         EntityType type,
                                                         Range &ents) const {

  MoFEMFunctionBegin;

  EntityHandle meshset = get_finite_element_meshset(name);
  CHKERR get_moab().get_entities_by_type(meshset, type, ents, true);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::get_finite_element_entities_by_handle(const std::string name,
                                            Range &ents) const {

  MoFEMFunctionBegin;

  EntityHandle meshset = get_finite_element_meshset(name);
  CHKERR get_moab().get_entities_by_handle(meshset, ents, true);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::list_finite_elements() const {
  MoFEMFunctionBeginHot;
  typedef FiniteElement_multiIndex::index<BitFEId_mi_tag>::type
      finiteElements_by_id;
  const finiteElements_by_id &BitFEId_set =
      finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = BitFEId_set.begin();
  for (; miit != BitFEId_set.end(); miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscSynchronizedPrintf(cOmm, ss.str().c_str());
  }
  PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::add_ents_to_finite_element_by_type(
    const EntityHandle meshset, const EntityType type, const std::string &name,
    const bool recursive) {
  *buildMoFEM &= 1 << 0;
  EntityHandle idm = no_handle;
  MoFEMFunctionBegin;

  idm = get_finite_element_meshset(getBitFEId(name));
  Range ents;
  CHKERR get_moab().get_entities_by_type(meshset, type, ents, recursive);
  CHKERR getInterface<BitRefManager>()->setElementsBitRefLevel(ents);
  CHKERR get_moab().add_entities(idm, ents);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::add_ents_to_finite_element_by_dim(const EntityHandle meshset,
                                        const int dim, const std::string &name,
                                        const bool recursive) {
  EntityHandle idm = no_handle;
  *buildMoFEM &= 1 << 0;
  MoFEMFunctionBegin;
  idm = get_finite_element_meshset(getBitFEId(name));
  Range ents;
  CHKERR get_moab().get_entities_by_dimension(meshset, dim, ents, recursive);
  CHKERR getInterface<BitRefManager>()->setElementsBitRefLevel(ents);
  CHKERR get_moab().add_entities(idm, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_ents_to_finite_element_by_type(
    const Range &ents, const EntityType type, const std::string &name) {
  EntityHandle idm = no_handle;
  *buildMoFEM &= 1 << 0;
  MoFEMFunctionBegin;
  idm = get_finite_element_meshset(getBitFEId(name));
  CHKERR getInterface<BitRefManager>()->setElementsBitRefLevel(
      ents.subset_by_type(type));
  CHKERR get_moab().add_entities(idm, ents.subset_by_type(type));
  MoFEMFunctionReturn(0);
} // namespace MoFEM

MoFEMErrorCode
Core::add_ents_to_finite_element_by_dim(const Range &ents, const int dim,
                                        const std::string &name) {
  EntityHandle idm = no_handle;
  *buildMoFEM &= 1 << 0;
  MoFEMFunctionBegin;
  idm = get_finite_element_meshset(getBitFEId(name));
  CHKERR getInterface<BitRefManager>()->setElementsBitRefLevel(
      ents.subset_by_dimension(dim));
  CHKERR get_moab().add_entities(idm, ents.subset_by_dimension(dim));
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_ents_to_finite_element_by_EDGEs(
    const EntityHandle meshset, const std::string &name, const bool recursive) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_type(meshset, MBEDGE, name, recursive);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode
Core::add_ents_to_finite_element_by_EDGEs(const Range &edges,
                                          const std::string &name) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_type(edges, MBEDGE, name);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode
Core::add_ents_to_finite_element_by_VERTICEs(const Range &vert,
                                             const std::string &name) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_type(vert, MBVERTEX, name);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode
Core::add_ents_to_finite_element_by_TRIs(const Range &tris,
                                         const std::string &name) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_type(tris, MBTRI, name);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode Core::add_ents_to_finite_element_by_TRIs(
    const EntityHandle meshset, const std::string &name, const bool recursive) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_type(meshset, MBTRI, name, recursive);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode
Core::add_ents_to_finite_element_by_TETs(const Range &tets,
                                         const std::string &name) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_type(tets, MBTET, name);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode Core::add_ents_to_finite_element_by_TETs(
    const EntityHandle meshset, const std::string &name, const bool recursive) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_type(meshset, MBTET, name, recursive);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode
Core::add_ents_to_finite_element_by_PRISMs(const Range &prims,
                                           const std::string &name) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_type(prims, MBPRISM, name);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode Core::add_ents_to_finite_element_by_PRISMs(
    const EntityHandle meshset, const std::string &name, const bool recursive) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_type(meshset, MBPRISM, name, recursive);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,
                                                    const std::string &name,
                                                    EntityType type, int verb) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_bit_ref(bit, BitRefLevel().set(), name,
                                               type, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_ents_to_finite_element_EntType_by_bit_ref(
    const BitRefLevel &bit, const BitRefLevel &mask, const std::string &name,
    EntityType type, int verb) {
  MoFEMFunctionBegin;
  CHKERR add_ents_to_finite_element_by_bit_ref(bit, mask, name, type, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_ents_to_finite_element_by_bit_ref(
    const BitRefLevel &bit, const BitRefLevel &mask, const std::string &name,
    EntityType type, int verb) {
  MoFEMFunctionBegin;
  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_TAG("SYNC", PETSC_FUNCTION_NAME);
  
  if (verb == -1)
    verb = verbose;
  *buildMoFEM &= 1 << 0;
  const BitFEId id = getBitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);

  auto &ref_MoFEMFiniteElement = refinedFiniteElements.get<EntType_mi_tag>();
  auto miit = ref_MoFEMFiniteElement.lower_bound(type);
  auto hi_miit = ref_MoFEMFiniteElement.upper_bound(type);

  int nb_add_fes = 0;
  for (; miit != hi_miit; miit++) {
    BitRefLevel bit2 = miit->get()->getBitRefLevel();
    if ((bit2 & mask) != bit2)
      continue;
    if ((bit2 & bit).any()) {
      EntityHandle ent = miit->get()->getRefEnt();
      CHKERR get_moab().add_entities(idm, &ent, 1);
      nb_add_fes++;
    }
  }

  MOFEM_LOG("SYNC", LogManager::SeverityLevel::inform)
      << "Finite element " << name << " added. Nb. of elements added "
      << nb_add_fes << " out of " << std::distance(miit, hi_miit);

  MOFEM_LOG_SYNCHORMISE(cOmm)

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_ents_to_finite_element_by_MESHSET(
    const EntityHandle meshset, const std::string &name, const bool recursive) {
  MoFEMFunctionBegin;
  *buildMoFEM &= 1 << 0;
  const BitFEId id = getBitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);
  if (recursive == false) {
    CHKERR get_moab().add_entities(idm, &meshset, 1);
  } else {
    Range meshsets;
    CHKERR get_moab().get_entities_by_type(meshset, MBENTITYSET, meshsets,
                                           false);
    CHKERR get_moab().add_entities(idm, meshsets);
  }
  MoFEMFunctionReturn(0);
}

template <int I> struct BuildFiniteElements {

  template <typename T1, typename T2>
  static inline void addToData(T1 &range_dit, T2 &fe_vec) {
    static_assert(I == DATA, "t should be set to DATA");

    for (auto dit = range_dit.first; dit != range_dit.second; ++dit) {
      const EntityHandle dof_ent = dit->get()->getEnt();
      // Fill array
      for (auto fe_it : fe_vec) {
        if (auto fe_ptr = fe_it.lock()) {
          // Add FEDofEntity, first create dofs, one by one, note that memory
          // is already reserved. Then create shared pointers and finally add
          // th_FEName to element multi-index
          // There are data dofs on this element
          auto &side_number_ptr = fe_ptr->getSideNumberPtr(dof_ent);
          fe_ptr->getDofsSequence().lock()->emplace_back(side_number_ptr, *dit);
        }
      }
    }
  }

  template <typename T> static inline void emplaceHint(T &fe_vec) {
    static_assert(I == DATA, "t should be set to DATA");

    // Add to data in FE
    for (auto fe_it : fe_vec) {
      if (auto fe_ptr = fe_it.lock()) {
        // It is a but unsafe, since if one mess up something in
        // buildFiniteElements, weak_ptr will not give pointer
        auto data_dofs_array_vec = fe_ptr->getDofsSequence().lock();
        // Create shared pointers vector
        auto hint = fe_ptr->data_dofs->end();
        for (auto &vit : *data_dofs_array_vec)
          hint =
              fe_ptr->data_dofs->emplace_hint(hint, data_dofs_array_vec, &vit);
      }
    }
  }
};

MoFEMErrorCode
Core::buildFiniteElements(const boost::shared_ptr<FiniteElement> &fe,
                          const Range *ents_ptr, int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  auto &fields_by_id = fIelds.get<BitFieldId_mi_tag>();

  // Get id of mofem fields for row, col and data
  enum IntLoop { ROW = 0, COL, DATA, LAST };
  std::array<BitFieldId, LAST> fe_fields = {fe.get()->getBitFieldIdRow(),
                                            fe.get()->getBitFieldIdCol(),
                                            fe.get()->getBitFieldIdData()};

  // Get finite element meshset
  EntityHandle meshset = get_finite_element_meshset(fe.get()->getId());

  // Get entities from finite element meshset // if meshset
  Range fe_ents;
  CHKERR get_moab().get_entities_by_handle(meshset, fe_ents, false);

  if (ents_ptr)
    fe_ents = intersect(fe_ents, *ents_ptr);

  // Map entity uid to pointers
  typedef std::vector<boost::weak_ptr<EntFiniteElement>> VecOfWeakFEPtrs;
  typedef std::map<const UId *, VecOfWeakFEPtrs> MapEntUIdAndVecOfWeakFEPtrs;
  MapEntUIdAndVecOfWeakFEPtrs ent_uid_and_fe_vec;
  std::map<EntityHandle, boost::shared_ptr<std::vector<FEDofEntity>>>
      data_dofs_array;
  VecOfWeakFEPtrs processed_fes;
  processed_fes.reserve(fe_ents.size());

  int last_row_field_ents_view_size = 0;
  int last_col_field_ents_view_size = 0;

  // View of field entities on element
  FieldEntity_vector_view data_field_ents_view;

  // Loop meshset finite element ents and add finite elements
  for (Range::const_pair_iterator peit = fe_ents.const_pair_begin();
       peit != fe_ents.const_pair_end(); peit++) {

    EntityHandle first = peit->first;
    EntityHandle second = peit->second;

    // Find range of ref entities that is sequence
    // note: iterator is a wrapper
    // check if is in refinedFiniteElements database
    auto ref_fe_miit =
        refinedFiniteElements.get<Ent_mi_tag>().lower_bound(first);
    if (ref_fe_miit == refinedFiniteElements.get<Ent_mi_tag>().end()) {
      std::ostringstream ss;
      ss << "refinedFiniteElements not in database ent = " << first;
      ss << " type " << get_moab().type_from_handle(first);
      ss << " " << *fe;
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, ss.str().c_str());
    }
    auto hi_ref_fe_miit =
        refinedFiniteElements.get<Ent_mi_tag>().upper_bound(second);

    EntFiniteElement_multiIndex::iterator hint_p = entsFiniteElements.end();
    for (; ref_fe_miit != hi_ref_fe_miit; ref_fe_miit++) {

      // Add finite element to database
      hint_p = entsFiniteElements.emplace_hint(
          hint_p, boost::make_shared<EntFiniteElement>(*ref_fe_miit, fe));
      processed_fes.emplace_back(*hint_p);
      auto fe_raw_ptr = hint_p->get();

      // Allocate space for etities view
      data_field_ents_view.clear();
      fe_raw_ptr->row_field_ents_view->reserve(last_row_field_ents_view_size);
      // Create shared pointer for entities view
      if (fe_fields[ROW] == fe_fields[COL]) {
        fe_raw_ptr->col_field_ents_view = fe_raw_ptr->row_field_ents_view;
      } else {
        // row and columns are diffent
        if (fe_raw_ptr->col_field_ents_view == fe_raw_ptr->row_field_ents_view)
          fe_raw_ptr->col_field_ents_view =
              boost::make_shared<FieldEntity_vector_view>();
        fe_raw_ptr->col_field_ents_view->reserve(last_col_field_ents_view_size);
      }

      int nb_dofs_on_data = 0;

      // Iterate over all field and check which one is on the element
      for (unsigned int ii = 0; ii != BitFieldId().size(); ++ii) {

        // Common field id for ROW, COL and DATA
        BitFieldId id_common = 0;
        // Check if the field (ii) is added to finite element
        for (int ss = 0; ss < LAST; ss++) {
          id_common |= fe_fields[ss] & BitFieldId().set(ii);
        }
        if (id_common.none())
          continue;

        // Find in database data associated with the field (ii)
        const BitFieldId field_id = BitFieldId().set(ii);
        auto miit = fields_by_id.find(field_id);
        if (miit == fields_by_id.end()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Data inconsistency");
        }

        // Loop over adjacencies of element and find field entities on those
        // adjacencies, that create hash map map_uid_fe which is used later
        const std::string field_name = miit->get()->getName();
        const bool add_to_data = (field_id & fe_fields[DATA]).any();
        const bool add_to_row = (field_id & fe_fields[ROW]).any();
        const bool add_to_col = (field_id & fe_fields[COL]).any();

        // Entities adjacent to entities
        Range adj_ents;

        // Resolve entities on element, those entities are used to build tag
        // with dof uids on finite element tag
        CHKERR fe_raw_ptr->getElementAdjacency(*miit, adj_ents);

        for (Range::pair_iterator p_eit = adj_ents.pair_begin();
             p_eit != adj_ents.pair_end(); ++p_eit) {

          const EntityHandle first = p_eit->first;
          const EntityHandle second = p_eit->second;

          typedef FieldEntity_multiIndex::index<
              Composite_Name_And_Ent_mi_tag>::type FieldEntityByComposite;
          auto &field_ents_by_name_and_ent =
              entsFields.get<Composite_Name_And_Ent_mi_tag>();
          FieldEntityByComposite::iterator meit;

          // If one entity in the pair search for one, otherwise search for
          // range
          if (first == second)
            meit = field_ents_by_name_and_ent.find(
                boost::make_tuple(field_name, first));
          else
            meit = field_ents_by_name_and_ent.lower_bound(
                boost::make_tuple(field_name, first));

          if (meit != field_ents_by_name_and_ent.end()) {

            decltype(meit) hi_meit;

            if (first == second) {
              hi_meit = meit;
              ++hi_meit;
            } else
              hi_meit = field_ents_by_name_and_ent.upper_bound(
                  boost::make_tuple(field_name, second));

            // Add to view and create list of finite elements with this dof UId
            for (; meit != hi_meit; ++meit) {
              // Add entity to map with key entity uids pointers  and data
              // finite elements weak ptrs. I using pointers to uids instead
              // uids because this is faster.
              const UId *uid_ptr = &(meit->get()->getGlobalUniqueId());
              auto &fe_vec = ent_uid_and_fe_vec[uid_ptr];
              // get number of dofs on entities to pre-allocate memory for
              // element
              const int nb_dofs_on_ent = meit->get()->getNbDofsOnEnt();
              if (add_to_data) {
                nb_dofs_on_data += nb_dofs_on_ent;
                data_field_ents_view.emplace_back(*meit);
              }
              if (add_to_row) {
                fe_raw_ptr->row_field_ents_view->emplace_back(*meit);
              }
              if (add_to_col) {
                if (fe_raw_ptr->col_field_ents_view !=
                    fe_raw_ptr->row_field_ents_view)
                  fe_raw_ptr->col_field_ents_view->emplace_back(*meit);
              }
              // add finite element to processed list
              fe_vec.emplace_back(*hint_p);
            }
          }
        }
      }

      // Sort field ents by uid
      auto uid_comp = [](const auto &a, const auto &b) {
        return a.lock()->getGlobalUniqueId() < b.lock()->getGlobalUniqueId();
      };

      // Sort all views

      // Data
      sort(data_field_ents_view.begin(), data_field_ents_view.end(), uid_comp);
      for (auto e : data_field_ents_view)
        fe_raw_ptr->data_field_ents_view->emplace_back(e);

      // Row
      sort(fe_raw_ptr->row_field_ents_view->begin(),
           fe_raw_ptr->row_field_ents_view->end(), uid_comp);
      last_row_field_ents_view_size = fe_raw_ptr->row_field_ents_view->size();

      // Column
      if (fe_raw_ptr->col_field_ents_view != fe_raw_ptr->row_field_ents_view) {
        sort(fe_raw_ptr->col_field_ents_view->begin(),
             fe_raw_ptr->col_field_ents_view->end(), uid_comp);
        last_col_field_ents_view_size = fe_raw_ptr->col_field_ents_view->size();
      }

      // Clear finite element data structures
      fe_raw_ptr->data_dofs->clear();

      // Reserve memory for data FE Dofs
      auto data_dofs_array_vec = boost::make_shared<std::vector<FEDofEntity>>();
      data_dofs_array[fe_raw_ptr->getEnt()] = data_dofs_array_vec;
      data_dofs_array_vec->reserve(nb_dofs_on_data);

      fe_raw_ptr->getDofsSequence() = data_dofs_array_vec;
    }
  }

  auto &dofs_by_ent_uid = dofsField.get<Unique_Ent_mi_tag>();

  // Loop over hash map, which has all entities on given elemnts
  boost::shared_ptr<SideNumber> side_number_ptr;
  for (auto &mit : ent_uid_and_fe_vec) {
    auto range_dit = dofs_by_ent_uid.equal_range(*mit.first);
    if (range_dit.first != range_dit.second) {
      const BitFieldId field_id = range_dit.first->get()->getId();
      if ((field_id & fe_fields[DATA]).any())
        BuildFiniteElements<DATA>::addToData(range_dit, mit.second);
    }
  }

  BuildFiniteElements<DATA>::emplaceHint(processed_fes);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_finite_elements(int verb) {
  MoFEMFunctionBeginHot;
  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_TAG("SYNC", PETSC_FUNCTION_NAME);
  MOFEM_LOG_CHANNEL("WORD");
  MOFEM_LOG_TAG("WORD", PETSC_FUNCTION_NAME);

  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  // loop Finite Elements
  for (auto &fe : finiteElements) {
    MOFEM_LOG("SYNC", LogManager::SeverityLevel::verbose)
        << "Build Finite Elements " << fe->getName();
    CHKERR buildFiniteElements(fe, NULL, verb);
  }

  if (verb >= VERBOSE) {

    auto &finite_elements_by_id = entsFiniteElements.get<BitFEId_mi_tag>();
    for (auto &fe : finiteElements) {
      auto miit = finite_elements_by_id.lower_bound(fe->getId());
      auto hi_miit = finite_elements_by_id.upper_bound(fe->getId());
      const auto count = std::distance(miit, hi_miit);
      MOFEM_LOG("SYNC", LogManager::SeverityLevel::inform)
          << "Finite element " << fe->getName()
          << " added. Nb. of elements added " << count;
    }

    MOFEM_LOG_SYNCHORMISE(cOmm)
  }

  *buildMoFEM |= 1 << 1;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::build_finite_elements(const BitRefLevel &bit, int verb) {
  MoFEMFunctionBeginHot;
  SETERRQ(cOmm, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::build_finite_elements(const string fe_name,
                                           const Range *const ents_ptr,
                                           int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;

  auto fe_miit = finiteElements.get<FiniteElement_name_mi_tag>().find(fe_name);
  if (fe_miit == finiteElements.get<FiniteElement_name_mi_tag>().end())
    SETERRQ1(cOmm, MOFEM_NOT_FOUND, "Finite element <%s> not found",
             fe_name.c_str());

  if (verb >= VERBOSE)
    PetscPrintf(cOmm, "Build Finite Elements %s\n", fe_name.c_str());
  CHKERR buildFiniteElements(*fe_miit, ents_ptr, verb);

  if (verb >= VERBOSE) {
    auto &finite_elements_by_id = entsFiniteElements.get<BitFEId_mi_tag>();
    auto miit = finite_elements_by_id.lower_bound((*fe_miit)->getId());
    auto hi_miit = finite_elements_by_id.upper_bound((*fe_miit)->getId());
    int count = std::distance(miit, hi_miit);
    std::ostringstream ss;
    ss << *(*fe_miit) << " Nb. FEs " << count << std::endl;
    PetscSynchronizedPrintf(cOmm, ss.str().c_str());
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  }

  *buildMoFEM |= 1 << 1;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::build_adjacencies(const Range &ents, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;
  MOFEM_LOG_CHANNEL("WORLD");
  MOFEM_LOG_TAG("WORLD", PETSC_FUNCTION_NAME)

  if (!((*buildMoFEM) & BUILD_FIELD))
    SETERRQ(cOmm, MOFEM_NOT_FOUND, "field not build");
  if (!((*buildMoFEM) & BUILD_FE))
    SETERRQ(cOmm, MOFEM_NOT_FOUND, "fe not build");
  for (Range::const_pair_iterator peit = ents.pair_begin();
       peit != ents.pair_end(); ++peit) {
    EntFiniteElement_multiIndex::index<Ent_mi_tag>::type::iterator fit, hi_fit;
    fit = entsFiniteElements.get<Ent_mi_tag>().lower_bound(peit->first);
    hi_fit = entsFiniteElements.get<Ent_mi_tag>().upper_bound(peit->second);
    for (; fit != hi_fit; ++fit) {
      if ((*fit)->getBitFieldIdRow().none() &&
          (*fit)->getBitFieldIdCol().none() &&
          (*fit)->getBitFieldIdData().none())
        continue;
      int by = BYROW;
      if ((*fit)->getBitFieldIdRow() != (*fit)->getBitFieldIdCol())
        by |= BYCOL;
      if ((*fit)->getBitFieldIdRow() != (*fit)->getBitFieldIdData())
        by |= BYDATA;
      FieldEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_row(by);
      auto hint = entFEAdjacencies.end();
      for (auto e : *(*fit)->row_field_ents_view) {
        hint = entFEAdjacencies.emplace_hint(hint, e.lock(), *fit);
        bool success = entFEAdjacencies.modify(hint, modify_row);
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
      }
      if ((*fit)->getBitFieldIdRow() != (*fit)->getBitFieldIdCol()) {
        int by = BYCOL;
        if ((*fit)->getBitFieldIdCol() != (*fit)->getBitFieldIdData())
          by |= BYDATA;
        FieldEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_col(by);
        auto hint = entFEAdjacencies.end();
        for (auto e : *(*fit)->col_field_ents_view) {
          hint = entFEAdjacencies.emplace_hint(hint, e.lock(), *fit);
          bool success = entFEAdjacencies.modify(hint, modify_col);
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
        }
      }
      if ((*fit)->getBitFieldIdRow() != (*fit)->getBitFieldIdData() ||
          (*fit)->getBitFieldIdCol() != (*fit)->getBitFieldIdData()) {
        FieldEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_data(
            BYDATA);
        auto hint = entFEAdjacencies.end();
        for (auto &e : *(*fit)->data_field_ents_view) {
          hint = entFEAdjacencies.emplace_hint(hint, e, *fit);
          bool success = entFEAdjacencies.modify(hint, modify_data);
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
        }
      }
    }
  }
  
  if (verb >= VERBOSE) {
    MOFEM_LOG("WORLD", LogManager::SeverityLevel::inform) <<
      "Number of adjacencies " << entFEAdjacencies.size();
    MOFEM_LOG_SYNCHORMISE(cOmm)
  }

  *buildMoFEM |= 1 << 2;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::build_adjacencies(const BitRefLevel &bit,
                                       const BitRefLevel &mask, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents);

  CHKERR build_adjacencies(ents, verb);

  MoFEMFunctionReturn(0);
}
MoFEMErrorCode Core::build_adjacencies(const BitRefLevel &bit, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR build_adjacencies(bit, BitRefLevel().set(), verb);

  MoFEMFunctionReturn(0);
}

EntFiniteElementByName::iterator
Core::get_fe_by_name_begin(const std::string &fe_name) const {
  return entsFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(
      fe_name);
}
EntFiniteElementByName::iterator
Core::get_fe_by_name_end(const std::string &fe_name) const {
  return entsFiniteElements.get<FiniteElement_name_mi_tag>().upper_bound(
      fe_name);
}

MoFEMErrorCode Core::check_number_of_ents_in_ents_finite_element(
    const std::string &name) const {
  MoFEMFunctionBegin;
  FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator it;
  it = finiteElements.get<FiniteElement_name_mi_tag>().find(name);
  if (it == finiteElements.get<FiniteElement_name_mi_tag>().end()) {
    SETERRQ1(cOmm, 1, "finite element not found < %s >", name.c_str());
  }
  EntityHandle meshset = (*it)->getMeshset();

  int num_entities;
  CHKERR get_moab().get_number_entities_by_handle(meshset, num_entities);

  if (entsFiniteElements.get<FiniteElement_name_mi_tag>().count(
          (*it)->getName().c_str()) != (unsigned int)num_entities) {
    SETERRQ1(cOmm, 1,
             "not equal number of entities in meshset and finite elements "
             "multiindex < %s >",
             (*it)->getName().c_str());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::check_number_of_ents_in_ents_finite_element() const {
  MoFEMFunctionBegin;
  FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator it;
  it = finiteElements.get<FiniteElement_name_mi_tag>().begin();
  for (; it != finiteElements.get<FiniteElement_name_mi_tag>().end(); it++) {
    EntityHandle meshset = (*it)->getMeshset();

    int num_entities;
    CHKERR get_moab().get_number_entities_by_handle(meshset, num_entities);

    if (entsFiniteElements.get<FiniteElement_name_mi_tag>().count(
            (*it)->getName().c_str()) != (unsigned int)num_entities) {
      SETERRQ1(cOmm, 1,
               "not equal number of entities in meshset and finite elements "
               "multiindex < %s >",
               (*it)->getName().c_str());
    }
  }
  MoFEMFunctionReturn(0);
}
} // namespace MoFEM
