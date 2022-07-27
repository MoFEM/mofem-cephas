/** \file FECore.cpp
 * \brief Core interface methods for managing deletions and insertion dofs
 */


#include <MoFEM.hpp>

#define FECoreFunctionBegin                                                    \
  MoFEMFunctionBegin;                                                          \
  MOFEM_LOG_CHANNEL("WORLD");                                                  \
  MOFEM_LOG_CHANNEL("SYNC");                                                   \
  MOFEM_LOG_FUNCTION();                                                        \
  MOFEM_LOG_TAG("WORLD", "FECore");                                            \
  MOFEM_LOG_TAG("SYNC", "FECore");

namespace MoFEM {

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
  FECoreFunctionBegin;
  *buildMoFEM &= 1 << 0;
  if (verb == -1) {
    verb = verbose;
  }

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

  auto &finite_element_name_set =
      finiteElements.get<FiniteElement_name_mi_tag>();
  auto it_fe = finite_element_name_set.find(fe_name);

  if (bh == MF_EXCL) {
    if (it_fe != finite_element_name_set.end()) {
      SETERRQ1(mofemComm, MOFEM_NOT_FOUND, "this < %s > is there",
               fe_name.c_str());
    }

  } else {
    if (it_fe != finite_element_name_set.end())
      MoFEMFunctionReturnHot(0);
  }
  EntityHandle meshset;
  CHKERR get_moab().create_meshset(MESHSET_SET, meshset);
  CHKERR add_meshset_to_partition(meshset);

  // id
  int fe_shift = 0;
  for (; finiteElements.get<BitFEId_mi_tag>().find(BitFEId().set(fe_shift)) !=
         finiteElements.get<BitFEId_mi_tag>().end();
       ++fe_shift) {
  }

  auto id = BitFEId().set(fe_shift);
  CHKERR get_moab().tag_set_data(th_FEId, &meshset, 1, &id);

  // id name
  void const *tag_data[] = {fe_name.c_str()};
  int tag_sizes[1];
  tag_sizes[0] = fe_name.size();
  CHKERR get_moab().tag_set_by_ptr(th_FEName, &meshset, 1, tag_data, tag_sizes);

  // add FiniteElement
  auto p = finiteElements.insert(
      boost::shared_ptr<FiniteElement>(new FiniteElement(moab, meshset)));
  if (!p.second)
    SETERRQ(mofemComm, MOFEM_OPERATION_UNSUCCESSFUL,
            "FiniteElement not inserted");

  if (verb > QUIET)
    MOFEM_LOG("WORLD", Sev::inform) << "Add finite element " << fe_name;

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
    SETERRQ(mofemComm, MOFEM_NOT_FOUND,
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
    SETERRQ(mofemComm, MOFEM_NOT_FOUND,
            "This finite element is not defined (advise: check spelling)");
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_change_bit_add(get_field_id(name_data)));
  if (!success)
    SETERRQ(mofemComm, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
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
    SETERRQ1(mofemComm, MOFEM_NOT_FOUND, "this < %s > is not there",
             fe_name.c_str());
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_row_change_bit_add(get_field_id(name_row)));
  if (!success)
    SETERRQ(mofemComm, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
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
    SETERRQ(mofemComm, MOFEM_OPERATION_UNSUCCESSFUL,
            "this FiniteElement is there");
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_col_change_bit_add(get_field_id(name_col)));
  if (!success)
    SETERRQ(mofemComm, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
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
    SETERRQ(mofemComm, MOFEM_NOT_FOUND, "this FiniteElement is there");
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_change_bit_off(get_field_id(name_data)));
  if (!success)
    SETERRQ(mofemComm, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
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
    SETERRQ1(mofemComm, MOFEM_NOT_FOUND, "this < %s > is not there",
             fe_name.c_str());
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_row_change_bit_off(get_field_id(name_row)));
  if (!success)
    SETERRQ(mofemComm, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
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
    SETERRQ(mofemComm, MOFEM_NOT_FOUND, "this FiniteElement is there");
  bool success = finite_element_name_set.modify(
      it_fe, FiniteElement_col_change_bit_off(get_field_id(name_col)));
  if (!success)
    SETERRQ(mofemComm, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

BitFEId Core::getBitFEId(const std::string &name) const {
  auto &fe_by_id = finiteElements.get<FiniteElement_name_mi_tag>();
  auto miit = fe_by_id.find(name);
  if (miit == fe_by_id.end())
    THROW_MESSAGE(
        ("finite element < " + name + " > not found (top tip: check spelling)")
            .c_str());
  return (*miit)->getId();
}

std::string Core::getBitFEIdName(const BitFEId id) const {
  auto &fe_by_id = finiteElements.get<BitFEId_mi_tag>();
  auto miit = fe_by_id.find(id);
  if (miit == fe_by_id.end())
    THROW_MESSAGE("finite element not found");
  return (*miit)->getName();
}

EntityHandle Core::get_finite_element_meshset(const BitFEId id) const {
  auto &fe_by_id = finiteElements.get<BitFEId_mi_tag>();
  auto miit = fe_by_id.find(id);
  if (miit == fe_by_id.end())
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
  MoFEMFunctionBegin;
  for (auto &fe : finiteElements.get<FiniteElement_name_mi_tag>())
    MOFEM_LOG("SYNC", Sev::inform) << fe;

  MOFEM_LOG_SYNCHRONISE(mofemComm);
  MoFEMFunctionReturn(0);
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
  FECoreFunctionBegin;

  if (verb == -1)
    verb = verbose;
  *buildMoFEM &= 1 << 0;
  const BitFEId id = getBitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);

  auto &ref_MoFEMFiniteElement = refinedFiniteElements.get<Ent_mi_tag>();
  auto miit = ref_MoFEMFiniteElement.lower_bound(get_id_for_min_type(type));
  auto hi_miit = ref_MoFEMFiniteElement.upper_bound(get_id_for_max_type(type));

  int nb_add_fes = 0;
  for (; miit != hi_miit; miit++) {
    const auto &bit2 = miit->get()->getBitRefLevel();
    if ((bit2 & mask) != bit2)
      continue;
    if ((bit2 & bit).any()) {
      EntityHandle ent = miit->get()->getEnt();
      CHKERR get_moab().add_entities(idm, &ent, 1);
      nb_add_fes++;
    }
  }

  MOFEM_LOG("SYNC", Sev::inform)
      << "Finite element " << name << " added. Nb. of elements added "
      << nb_add_fes << " out of " << std::distance(miit, hi_miit);

  MOFEM_LOG_SYNCHRONISE(mofemComm)

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

MoFEMErrorCode
Core::buildFiniteElements(const boost::shared_ptr<FiniteElement> &fe,
                          const Range *ents_ptr, int verb) {
  FECoreFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  if (verb > QUIET)
    MOFEM_LOG("SYNC", Sev::verbose)
        << "Build Finite Elements " << fe->getName();

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
  VecOfWeakFEPtrs processed_fes;
  processed_fes.reserve(fe_ents.size());

  int last_data_field_ents_view_size = 0;
  int last_row_field_ents_view_size = 0;
  int last_col_field_ents_view_size = 0;

  // Entities adjacent to entities
  std::vector<EntityHandle> adj_ents;

  // Loop meshset finite element ents and add finite elements
  for (Range::const_pair_iterator peit = fe_ents.const_pair_begin();
       peit != fe_ents.const_pair_end(); peit++) {

    const auto first = peit->first;
    const auto second = peit->second;

    // Find range of ref entities that is sequence
    // note: iterator is a wrapper
    // check if is in refinedFiniteElements database
    auto ref_fe_miit =
        refinedFiniteElements.get<Ent_mi_tag>().lower_bound(first);
    if (ref_fe_miit == refinedFiniteElements.get<Ent_mi_tag>().end()) {
      std::ostringstream ss;
      ss << "refinedFiniteElements not in database ent = " << first << " type "
         << type_from_handle << " " << *fe;
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
      bool row_as_data = false, col_as_row = false;
      if (fe_fields[DATA] == fe_fields[ROW])
        row_as_data = true;
      if (fe_fields[ROW] == fe_fields[COL])
        col_as_row = true;

      fe_raw_ptr->getDataFieldEntsPtr()->reserve(
          last_data_field_ents_view_size);

      if (row_as_data) {
        fe_raw_ptr->getRowFieldEntsPtr() = fe_raw_ptr->getDataFieldEntsPtr();
      } else {
        // row and col are diffent
        if (fe_raw_ptr->getRowFieldEntsPtr() ==
            fe_raw_ptr->getDataFieldEntsPtr())
          fe_raw_ptr->getRowFieldEntsPtr() =
              boost::make_shared<FieldEntity_vector_view>();
        fe_raw_ptr->getRowFieldEntsPtr()->reserve(
            last_row_field_ents_view_size);
      }

      if (row_as_data && col_as_row) {
        fe_raw_ptr->getColFieldEntsPtr() = fe_raw_ptr->getDataFieldEntsPtr();
      } else if (col_as_row) {
        fe_raw_ptr->getColFieldEntsPtr() = fe_raw_ptr->getRowFieldEntsPtr();
      } else {
        if (

            fe_raw_ptr->getColFieldEntsPtr() ==
                fe_raw_ptr->getRowFieldEntsPtr() ||
            fe_raw_ptr->getColFieldEntsPtr() ==
                fe_raw_ptr->getDataFieldEntsPtr()

        )
          fe_raw_ptr->getColFieldEntsPtr() =
              boost::make_shared<FieldEntity_vector_view>();
        fe_raw_ptr->getColFieldEntsPtr()->reserve(
            last_col_field_ents_view_size);
      }

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
        if (miit == fields_by_id.end())
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Field not found");
        auto field_bit_number = (*miit)->getBitNumber();

        // Loop over adjacencies of element and find field entities on those
        // adjacencies, that create hash map map_uid_fe which is used later
        const std::string field_name = miit->get()->getName();
        const bool add_to_data = (field_id & fe_fields[DATA]).any();
        const bool add_to_row = (field_id & fe_fields[ROW]).any();
        const bool add_to_col = (field_id & fe_fields[COL]).any();

        // Resolve entities on element, those entities are used to build tag
        // with dof uids on finite element tag
        adj_ents.clear();
        CHKERR fe_raw_ptr->getElementAdjacency(*miit, adj_ents);

        for(auto ent : adj_ents) {

          auto dof_it = entsFields.get<Unique_mi_tag>().find(
              FieldEntity::getLocalUniqueIdCalculate(field_bit_number, ent));
          if(dof_it!=entsFields.get<Unique_mi_tag>().end()) {
              // Add entity to map with key entity uids pointers  and data
              // finite elements weak ptrs. I using pointers to uids instead
              // uids because this is faster.
              const UId *uid_ptr = &(dof_it->get()->getLocalUniqueId());
              auto &fe_vec = ent_uid_and_fe_vec[uid_ptr];
              if (add_to_data) {
                fe_raw_ptr->getDataFieldEntsPtr()->emplace_back(*dof_it);
              }
              if (add_to_row && !row_as_data) {
                fe_raw_ptr->getRowFieldEntsPtr()->emplace_back(*dof_it);
              }
              if (add_to_col && !col_as_row) {
                fe_raw_ptr->getColFieldEntsPtr()->emplace_back(*dof_it);
              }

              // add finite element to processed list
              fe_vec.emplace_back(*hint_p);
          }

        }


      }

      // Sort field ents by uid
      auto uid_comp = [](const auto &a, const auto &b) {
        return a.lock()->getLocalUniqueId() < b.lock()->getLocalUniqueId();
      };

      // Sort all views

      // Data
      sort(fe_raw_ptr->getDataFieldEntsPtr()->begin(),
           fe_raw_ptr->getDataFieldEntsPtr()->end(), uid_comp);
      last_data_field_ents_view_size =
          fe_raw_ptr->getDataFieldEntsPtr()->size();

      // Row
      if (!row_as_data) {
        sort(fe_raw_ptr->getRowFieldEntsPtr()->begin(),
             fe_raw_ptr->getRowFieldEntsPtr()->end(), uid_comp);
        last_row_field_ents_view_size =
            fe_raw_ptr->getRowFieldEntsPtr()->size();
      }

      // Column
      if (!col_as_row) {
        sort(fe_raw_ptr->getColFieldEntsPtr()->begin(),
             fe_raw_ptr->getColFieldEntsPtr()->end(), uid_comp);
        last_col_field_ents_view_size =
            fe_raw_ptr->getColFieldEntsPtr()->size();
      }
    }
  }

  auto &dofs_by_ent_uid = dofsField.get<Unique_mi_tag>();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_finite_elements(int verb) {
  FECoreFunctionBegin;

  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  // loop Finite Elements
  for (auto &fe : finiteElements)
    CHKERR buildFiniteElements(fe, NULL, verb);

  if (verb > QUIET) {

    auto &fe_ents = entsFiniteElements.get<FiniteElement_name_mi_tag>();
    for (auto &fe : finiteElements) {
      auto miit = fe_ents.lower_bound(fe->getName());
      auto hi_miit = fe_ents.upper_bound(fe->getName());
      const auto count = std::distance(miit, hi_miit);
      MOFEM_LOG("SYNC", Sev::inform)
          << "Finite element " << fe->getName()
          << " added. Nb. of elements added " << count;
      MOFEM_LOG("SYNC", Sev::noisy) << *fe;

      auto slg = MoFEM::LogManager::getLog("SYNC");
      for (auto &field : fIelds) {
        auto rec = slg.open_record(keywords::severity = Sev::verbose);
        if (rec) {
          logging::record_ostream strm(rec);
          strm << "Field " << field->getName() << " on finite element: ";
          if ((field->getId() & fe->getBitFieldIdRow()).any())
            strm << "row ";
          if ((field->getId() & fe->getBitFieldIdCol()).any())
            strm << "columns ";
          if ((field->getId() & fe->getBitFieldIdData()).any())
            strm << "data";
          strm.flush();
          slg.push_record(boost::move(rec));
        }
      }
    }

    MOFEM_LOG_SYNCHRONISE(mofemComm);
  }

  *buildMoFEM |= 1 << 1;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_finite_elements(const BitRefLevel &bit, int verb) {
  MoFEMFunctionBeginHot;
  SETERRQ(mofemComm, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
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
    SETERRQ1(mofemComm, MOFEM_NOT_FOUND, "Finite element <%s> not found",
             fe_name.c_str());

  CHKERR buildFiniteElements(*fe_miit, ents_ptr, verb);

  if (verb >= VERBOSE) {
    auto &fe_ents = entsFiniteElements.get<FiniteElement_name_mi_tag>();
    auto miit = fe_ents.lower_bound((*fe_miit)->getName());
    auto hi_miit = fe_ents.upper_bound((*fe_miit)->getName());
    const auto count = std::distance(miit, hi_miit);
    MOFEM_LOG("SYNC", Sev::inform) << "Finite element " << fe_name
                                   << " added. Nb. of elements added " << count;
    MOFEM_LOG_SYNCHRONISE(mofemComm);
  }

  *buildMoFEM |= 1 << 1;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::build_adjacencies(const Range &ents, int verb) {
  FECoreFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  if (!((*buildMoFEM) & BUILD_FIELD))
    SETERRQ(mofemComm, MOFEM_NOT_FOUND, "field not build");
  if (!((*buildMoFEM) & BUILD_FE))
    SETERRQ(mofemComm, MOFEM_NOT_FOUND, "fe not build");
  for (auto peit = ents.pair_begin(); peit != ents.pair_end(); ++peit) {
    auto fit = entsFiniteElements.get<Ent_mi_tag>().lower_bound(peit->first);
    auto hi_fit = entsFiniteElements.get<Ent_mi_tag>().upper_bound(peit->second);
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
      for (auto e : *(*fit)->getRowFieldEntsPtr()) {
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
        for (auto e : *(*fit)->getColFieldEntsPtr()) {
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
        for (auto &e : (*fit)->getDataFieldEnts()) {
          hint = entFEAdjacencies.emplace_hint(hint, e.lock(), *fit);
          bool success = entFEAdjacencies.modify(hint, modify_data);
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
        }
      }
    }
  }

  if (verb >= VERBOSE) {
    MOFEM_LOG("WORLD", Sev::inform)
        << "Number of adjacencies " << entFEAdjacencies.size();
    MOFEM_LOG_SYNCHRONISE(mofemComm)
  }

  *buildMoFEM |= 1 << 2;
  MoFEMFunctionReturn(0);
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
    SETERRQ1(mofemComm, 1, "finite element not found < %s >", name.c_str());
  }
  EntityHandle meshset = (*it)->getMeshset();

  int num_entities;
  CHKERR get_moab().get_number_entities_by_handle(meshset, num_entities);

  if (entsFiniteElements.get<FiniteElement_name_mi_tag>().count(
          (*it)->getName().c_str()) != (unsigned int)num_entities) {
    SETERRQ1(mofemComm, 1,
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
      SETERRQ1(mofemComm, 1,
               "not equal number of entities in meshset and finite elements "
               "multiindex < %s >",
               (*it)->getName().c_str());
    }
  }
  MoFEMFunctionReturn(0);
}
} // namespace MoFEM
