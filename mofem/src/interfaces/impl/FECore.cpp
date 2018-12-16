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
    // list_finiteElements();
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
  if (it_fe == finite_element_name_set.end()) {
    SETERRQ(cOmm, MOFEM_NOT_FOUND, "this FiniteElement is there");
  }
  boost::shared_ptr<FiniteElement> fe;
  fe = *it_fe;
  fe->elementAdjacencyTable[type] = function;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_finite_element_add_field_data(const std::string &fe_name,
                                           const std::string &name_data) {
  MoFEMFunctionBeginHot;
  *buildMoFEM &= 1 << 0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe =
      finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ(cOmm, MOFEM_NOT_FOUND, "this FiniteElement is there");
  try {
    bool success = finite_element_name_set.modify(
        it_fe, FiniteElement_change_bit_add(getBitFieldId(name_data)));
    if (!success)
      SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm, e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_finite_element_add_field_row(const std::string &fe_name,
                                          const std::string &name_row) {
  MoFEMFunctionBeginHot;
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
  try {
    bool success = finite_element_name_set.modify(
        it_fe, FiniteElement_row_change_bit_add(getBitFieldId(name_row)));
    if (!success)
      SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm, e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_finite_element_add_field_col(const std::string &fe_name,
                                          const std::string &name_col) {
  MoFEMFunctionBeginHot;
  *buildMoFEM &= 1 << 0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe =
      finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "this FiniteElement is there");
  try {
    bool success = finite_element_name_set.modify(
        it_fe, FiniteElement_col_change_bit_add(getBitFieldId(name_col)));
    if (!success)
      SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm, e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_finite_element_off_field_data(const std::string &fe_name,
                                           const std::string &name_data) {
  MoFEMFunctionBeginHot;
  *buildMoFEM &= 1 << 0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe =
      finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ(cOmm, MOFEM_NOT_FOUND, "this FiniteElement is there");
  try {
    bool success = finite_element_name_set.modify(
        it_fe, FiniteElement_change_bit_off(getBitFieldId(name_data)));
    if (!success)
      SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm, e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_finite_element_off_field_row(const std::string &fe_name,
                                          const std::string &name_row) {
  MoFEMFunctionBeginHot;
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
  try {
    bool success = finite_element_name_set.modify(
        it_fe, FiniteElement_row_change_bit_off(getBitFieldId(name_row)));
    if (!success)
      SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm, e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_finite_element_off_field_col(const std::string &fe_name,
                                          const std::string &name_col) {
  MoFEMFunctionBeginHot;
  *buildMoFEM &= 1 << 0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type
      FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe =
      finite_element_name_set.find(fe_name);
  if (it_fe == finite_element_name_set.end())
    SETERRQ(cOmm, MOFEM_NOT_FOUND, "this FiniteElement is there");
  try {
    bool success = finite_element_name_set.modify(
        it_fe, FiniteElement_col_change_bit_off(getBitFieldId(name_col)));
    if (!success)
      SETERRQ(cOmm, MOFEM_OPERATION_UNSUCCESSFUL, "modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm, e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
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

  if (verb == -1)
    verb = verbose;
  *buildMoFEM &= 1 << 0;
  const BitFEId id = getBitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);
  typedef RefElement_multiIndex::index<EntType_mi_tag>::type refMoabFE_by_type;
  refMoabFE_by_type &ref_MoFEMFiniteElement =
      refinedFiniteElements.get<EntType_mi_tag>();
  refMoabFE_by_type::iterator miit = ref_MoFEMFiniteElement.lower_bound(type);
  refMoabFE_by_type::iterator hi_miit =
      ref_MoFEMFiniteElement.upper_bound(type);
  if (verb > 1) {
    PetscSynchronizedPrintf(cOmm, "nb. ref elements in database %d\n",
                            std::distance(miit, hi_miit));
  }
  int nb_add_FEs = 0;
  for (; miit != hi_miit; miit++) {
    BitRefLevel bit2 = miit->get()->getBitRefLevel();
    if ((bit2 & mask) != bit2)
      continue;
    if ((bit2 & bit).any()) {
      EntityHandle ent = miit->get()->getRefEnt();
      CHKERR get_moab().add_entities(idm, &ent, 1);
      nb_add_FEs++;
    }
  }
  if (verb > 0) {
    std::ostringstream ss;
    ss << "Add Nb. FEs " << nb_add_FEs << " form BitRef " << bit << std::endl;
    PetscSynchronizedPrintf(cOmm, "%s", ss.str().c_str());
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  }

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

template <int I>
struct AddToFE {

  template <typename T1, typename T2, typename T3>
  static inline void addToView(T1 &range_dit, T2 &fe_vec, T3 &dofs_size) {

    static_assert(I == ROW || I == COL, "t should be set to ROW or COL");

    auto uid_comp = [](auto a, auto b) {
      return a.lock()->getGlobalUniqueId() < b.lock()->getGlobalUniqueId();
    };

    // Reserve size
    for (auto dit = range_dit.first; dit != range_dit.second; ++dit)
      for (auto fe_it : fe_vec) {
        auto fe_raw_ptr = fe_it.lock().get();
        if (I == ROW)
          fe_raw_ptr->row_dof_view->reserve(dofs_size.at(fe_raw_ptr->getEnt()));
        else
          fe_raw_ptr->col_dof_view->reserve(dofs_size.at(fe_raw_ptr->getEnt()));
      }

    // Add DOFs
    for (auto dit = range_dit.first; dit != range_dit.second; ++dit)
      for (auto fe_it : fe_vec) {
        auto fe_raw_ptr = fe_it.lock().get();
        if (I == ROW)
          fe_raw_ptr->row_dof_view->emplace_back(*dit);
        else
          fe_raw_ptr->col_dof_view->emplace_back(*dit);
      }

    // Sort
    for (auto dit = range_dit.first; dit != range_dit.second; ++dit)
      for (auto fe_it : fe_vec) {
        auto fe_raw_ptr = fe_it.lock().get();
        if (I == ROW)
          sort(fe_raw_ptr->row_dof_view->begin(),
               fe_raw_ptr->row_dof_view->end(), uid_comp);
        else
          sort(fe_raw_ptr->row_dof_view->begin(),
               fe_raw_ptr->row_dof_view->end(), uid_comp);
      }
  }
};

MoFEMErrorCode
Core::buildFiniteElements(const boost::shared_ptr<FiniteElement> &fe,
                          const Range *ents_ptr, int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  typedef RefElement_multiIndex::index<Ent_mi_tag>::type RefFiniteElementByEnt;
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldById;
  FieldById &fields_by_id = fIelds.get<BitFieldId_mi_tag>();

  // get id of mofem fields for row, col and data
  enum IntLoop { ROW = 0, COL, DATA, LAST };
  BitFieldId fe_fields[LAST] = {fe.get()->getBitFieldIdRow(),
                                fe.get()->getBitFieldIdCol(),
                                fe.get()->getBitFieldIdData()};

  // get finite element meshset
  EntityHandle meshset = get_finite_element_meshset(fe.get()->getId());
  // get entities from finite element meshset // if meshset
  Range fe_ents;
  CHKERR get_moab().get_entities_by_handle(meshset, fe_ents, false);

  if (ents_ptr)
    fe_ents = intersect(fe_ents, *ents_ptr);

  // map entity uid to pointers
  typedef std::vector<boost::weak_ptr<EntFiniteElement>> VecOfWeakFEPtrs;
  typedef std::map<const UId *, VecOfWeakFEPtrs> MapEntUIdAndVecOfWeakFEPtrs;
  MapEntUIdAndVecOfWeakFEPtrs ent_uid_and_fe_vec;
  std::map<EntityHandle, int> data_dofs_size, row_dofs_size, col_dofs_size;

  // loop meshset Ents and add finite elements
  for (Range::const_pair_iterator peit = fe_ents.const_pair_begin();
       peit != fe_ents.const_pair_end(); peit++) {

    EntityHandle first = peit->first;
    EntityHandle second = peit->second;

    // note: iterator is a wrapper
    // check if is in refinedFiniteElements database
    RefFiniteElementByEnt::iterator ref_fe_miit, hi_ref_fe_miit;
    ref_fe_miit = refinedFiniteElements.get<Ent_mi_tag>().lower_bound(first);
    if (ref_fe_miit == refinedFiniteElements.get<Ent_mi_tag>().end()) {
      std::ostringstream ss;
      ss << "refinedFiniteElements not in database ent = " << first;
      ss << " type " << get_moab().type_from_handle(first);
      ss << " " << *fe;
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, ss.str().c_str());
    }
    hi_ref_fe_miit =
        refinedFiniteElements.get<Ent_mi_tag>().upper_bound(second);

    EntFiniteElement_multiIndex::iterator hint_p = entsFiniteElements.end();

    for (; ref_fe_miit != hi_ref_fe_miit; ref_fe_miit++) {

      hint_p = entsFiniteElements.emplace_hint(
          hint_p, boost::make_shared<EntFiniteElement>(*ref_fe_miit, fe));

      auto fe_raw_ptr = hint_p->get();

      if (fe_fields[ROW] == fe_fields[COL]) {
        fe_raw_ptr->col_dof_view = fe_raw_ptr->row_dof_view;
      } else if (fe_raw_ptr->col_dof_view == fe_raw_ptr->row_dof_view) {
        fe_raw_ptr->col_dof_view =
            boost::make_shared<DofEntity_vector_view>();
      }

      fe_raw_ptr->row_dof_view->clear();
      fe_raw_ptr->col_dof_view->clear();
      fe_raw_ptr->data_dofs->clear();

      for (unsigned int ii = 0; ii < BitFieldId().size(); ii++) {

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
        FieldById::iterator miit = fields_by_id.find(field_id);
        if (miit == fields_by_id.end()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Data inconsistency");
        }

        // Entities adjacent to entities
        Range adj_ents;

        // Resolve entities on element, those entities are used to build tag
        // with dof uids on finite element tag
        CHKERR fe_raw_ptr->getElementAdjacency(*miit, adj_ents);

        // Loop over adjacencies of element and find field entities on those
        // adjacencies, that create hash map map_uid_fe which is used later
        const std::string field_name = miit->get()->getName();
        const bool add_to_data = (field_id & fe_fields[DATA]).any();
        const bool add_to_row = (field_id & fe_fields[ROW]).any();
        const bool add_to_col = (field_id & fe_fields[COL]).any();

        for (Range::pair_iterator p_eit = adj_ents.pair_begin();
             p_eit != adj_ents.pair_end(); ++p_eit) {
          const EntityHandle first = p_eit->first;
          const EntityHandle second = p_eit->second;
          auto meit =
              entsFields.get<Composite_Name_And_Ent_mi_tag>().lower_bound(
                  boost::make_tuple(field_name, first));
          if (meit == entsFields.get<Composite_Name_And_Ent_mi_tag>().end())
            continue;
          auto hi_meit =
              entsFields.get<Composite_Name_And_Ent_mi_tag>().upper_bound(
                  boost::make_tuple(field_name, second));
          // create list of finite elements with this dof UId
          for (; meit != hi_meit; ++meit) {
            const UId *uid_ptr = &(meit->get()->getGlobalUniqueId());
            auto &fe_vec = ent_uid_and_fe_vec[uid_ptr];
            fe_vec.emplace_back(*hint_p);
            if (add_to_data)
              data_dofs_size[fe_raw_ptr->getEnt()] +=
                  meit->get()->getNbDofsOnEnt();
            if (add_to_row)
              row_dofs_size[fe_raw_ptr->getEnt()] +=
                  meit->get()->getNbDofsOnEnt();
            if (add_to_col)
              col_dofs_size[fe_raw_ptr->getEnt()] +=
                  meit->get()->getNbDofsOnEnt();
          }
        }
      }
    }
  }

  // Reserve memory
  std::map<EntityHandle, boost::shared_ptr<std::vector<FEDofEntity>>>
      data_dofs_array;
  for (auto &mit : data_dofs_size) {
    data_dofs_array[mit.first] = boost::make_shared<std::vector<FEDofEntity>>();
    data_dofs_array[mit.first]->reserve(mit.second);
  }

  typedef DofEntity_multiIndex::index<Unique_Ent_mi_tag>::type DofsByEntUId;
  DofsByEntUId &dofs_by_ent_uid = dofsField.get<Unique_Ent_mi_tag>();

  // Loop over hash map, which has all entities on given elemnts
  boost::shared_ptr<SideNumber> side_number_ptr;
  for (auto &mit : ent_uid_and_fe_vec) {

    auto range_dit = dofs_by_ent_uid.equal_range(*mit.first);

    if (range_dit.first != range_dit.second) {

      const BitFieldId field_id = range_dit.first->get()->getId();

      if ((field_id & fe_fields[ROW]).any())
        AddToFE<ROW>::addToView(range_dit, mit.second, row_dofs_size);

      if (fe_fields[ROW] != fe_fields[COL] && (field_id & fe_fields[COL]).any())
        AddToFE<COL>::addToView(range_dit, mit.second, col_dofs_size);

      for (auto dit = range_dit.first; dit != range_dit.second; ++dit) {
        const BitFieldId field_id = dit->get()->getId();
        const EntityHandle dof_ent = dit->get()->getEnt();

        for (auto fe_it : mit.second) {

          auto fe_raw_ptr = fe_it.lock().get();

          // Add FEDofEntity, first create dofs, one by one, note that memory
          // is already reserved. Then create shared pointers and finally add
          // th_FEName to element multi-index
          const EntityHandle fe_ent = fe_raw_ptr->getEnt();
          boost::shared_ptr<std::vector<FEDofEntity>> &data_dofs_array_vec =
              data_dofs_array[fe_ent];
          if (data_dofs_size[fe_ent] != 0 &&
              (field_id & fe_raw_ptr->getBitFieldIdData()).any()) {

            // There are data dofs on this element
            side_number_ptr = fe_raw_ptr->getSideNumberPtr(dof_ent);
            data_dofs_array_vec->emplace_back(side_number_ptr, *dit);
            if (data_dofs_array_vec->size() ==
                (unsigned int)data_dofs_size[fe_ent]) {
              // That means that FEDofEntity vector is full, and can be added to
              // multi-index

              // Create shared pointers vector
              auto hint = fe_raw_ptr->data_dofs->end();
              for (std::vector<FEDofEntity>::iterator vit =
                       data_dofs_array_vec->begin();
                   vit != data_dofs_array_vec->end(); vit++) {
                hint = fe_raw_ptr->data_dofs->emplace_hint(
                    hint, data_dofs_array_vec, &(*vit));
              }
              fe_raw_ptr->getDofsSequence() = data_dofs_array_vec;
            }
          }
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_finite_elements(int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;

  FiniteElement_multiIndex::iterator fe_miit = finiteElements.begin();

  // loop Finite Elements
  for (; fe_miit != finiteElements.end(); fe_miit++) {
    if (verb > 0)
      PetscPrintf(cOmm, "Build Finite Elements %s\n",
                  (*fe_miit)->getName().c_str());
    CHKERR buildFiniteElements(*fe_miit, NULL, verb);
  }

  if (verb > 0) {
    PetscSynchronizedPrintf(cOmm, "Nb. FEs %u\n", entsFiniteElements.size());
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
    typedef EntFiniteElement_multiIndex::index<BitFEId_mi_tag>::type
        FiniteElementById;
    FiniteElementById &finite_elements_by_id =
        entsFiniteElements.get<BitFEId_mi_tag>();
    FiniteElement_multiIndex::iterator fe_id_it = finiteElements.begin();
    for (; fe_id_it != finiteElements.end(); fe_id_it++) {
      FiniteElementById::iterator miit =
          finite_elements_by_id.lower_bound((*fe_id_it)->getId());
      FiniteElementById::iterator hi_miit =
          finite_elements_by_id.upper_bound((*fe_id_it)->getId());
      int count = std::distance(miit, hi_miit);
      std::ostringstream ss;
      ss << *(*fe_id_it) << " Nb. FEs " << count << std::endl;
      PetscSynchronizedPrintf(cOmm, ss.str().c_str());
      PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
    }
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
                                           const Range *ents_ptr, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;

  FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator
      fe_miit = finiteElements.get<FiniteElement_name_mi_tag>().find(fe_name);
  if (fe_miit == finiteElements.get<FiniteElement_name_mi_tag>().end()) {
    SETERRQ1(cOmm, MOFEM_NOT_FOUND, "Finite element <%s> not found",
             fe_name.c_str());
  }

  if (verb >= VERBOSE)
    PetscPrintf(cOmm, "Build Finite Elements %s\n", fe_name.c_str());
  CHKERR buildFiniteElements(*fe_miit, ents_ptr, verb);

  if (verb >= VERBOSE) {
    typedef EntFiniteElement_multiIndex::index<BitFEId_mi_tag>::type
        FiniteElementById;
    FiniteElementById &finite_elements_by_id =
        entsFiniteElements.get<BitFEId_mi_tag>();
    FiniteElementById::iterator miit =
        finite_elements_by_id.lower_bound((*fe_miit)->getId());
    FiniteElementById::iterator hi_miit =
        finite_elements_by_id.upper_bound((*fe_miit)->getId());
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
  if (verb == -1)
    verb = verbose;
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
      UId ent_uid = UId(0);
      for (auto rvit = (*fit)->row_dof_view->begin();
           rvit != (*fit)->row_dof_view->end(); ++rvit) {
        if (ent_uid == (*rvit).lock()->getFieldEntityPtr()->getGlobalUniqueId())
          continue;
        ent_uid = (*rvit).lock()->getFieldEntityPtr()->getGlobalUniqueId();
        std::pair<FieldEntityEntFiniteElementAdjacencyMap_multiIndex::iterator,
                  bool>
            p;
        p = entFEAdjacencies.insert(FieldEntityEntFiniteElementAdjacencyMap(
            (*rvit).lock()->getFieldEntityPtr(), *fit));
        bool success = entFEAdjacencies.modify(p.first, modify_row);
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
      }
      if ((*fit)->getBitFieldIdRow() != (*fit)->getBitFieldIdCol()) {
        int by = BYCOL;
        if ((*fit)->getBitFieldIdCol() != (*fit)->getBitFieldIdData())
          by |= BYDATA;
        FieldEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_col(by);
        ent_uid = UId(0);
        for (auto cvit = (*fit)->col_dof_view->begin();
             cvit != (*fit)->col_dof_view->end(); cvit++) {
          if (ent_uid ==
              (*cvit).lock()->getFieldEntityPtr()->getGlobalUniqueId())
            continue;
          ent_uid = (*cvit).lock()->getFieldEntityPtr()->getGlobalUniqueId();
          std::pair<
              FieldEntityEntFiniteElementAdjacencyMap_multiIndex::iterator,
              bool>
              p;
          p = entFEAdjacencies.insert(FieldEntityEntFiniteElementAdjacencyMap(
              (*cvit).lock()->getFieldEntityPtr(), *fit));
          bool success = entFEAdjacencies.modify(p.first, modify_col);
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
        }
      }
      if ((*fit)->getBitFieldIdRow() != (*fit)->getBitFieldIdData() ||
          (*fit)->getBitFieldIdCol() != (*fit)->getBitFieldIdData()) {
        FieldEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_data(
            BYDATA);
        ent_uid = UId(0);
        for (FEDofEntity_multiIndex::iterator dvit = (*fit)->data_dofs->begin();
             dvit != (*fit)->data_dofs->end(); dvit++) {
          if (ent_uid == (*dvit)->getFieldEntityPtr()->getGlobalUniqueId())
            continue;
          ent_uid = (*dvit)->getFieldEntityPtr()->getGlobalUniqueId();
          std::pair<
              FieldEntityEntFiniteElementAdjacencyMap_multiIndex::iterator,
              bool>
              p;
          p = entFEAdjacencies.insert(FieldEntityEntFiniteElementAdjacencyMap(
              (*dvit)->getFieldEntityPtr(), *fit));
          bool success = entFEAdjacencies.modify(p.first, modify_data);
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
        }
      }
    }
  }
  if (verb >= VERY_NOISY) {
    list_adjacencies();
  }
  if (verb >= VERBOSE) {
    PetscSynchronizedPrintf(cOmm, "Nb. entFEAdjacencies %u\n",
                            entFEAdjacencies.size());
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
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
