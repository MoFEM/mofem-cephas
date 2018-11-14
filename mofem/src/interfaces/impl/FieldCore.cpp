/** \file FieldCore.cpp
 * \brief Core interface methods for managing fields.
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

BitFieldId Core::getBitFieldId(const std::string &name) const {
  auto &set = fIelds.get<FieldName_mi_tag>();
  auto miit = set.find(name);
  if (miit == set.end()) {
    THROW_MESSAGE("field < " + name +
                  " > not in database (top tip: check spelling)");
  }
  return (*miit)->getId();
}

EntityHandle Core::get_field_meshset(const BitFieldId id) const {
  auto &set = fIelds.get<BitFieldId_mi_tag>();
  auto miit = set.find(id);
  if (miit == set.end())
    THROW_MESSAGE("field not in database (top tip: check spelling)");
  return (*miit)->meshSet;
}

EntityHandle Core::get_field_meshset(const std::string &name) const {
  return get_field_meshset(getBitFieldId(name));
}

bool Core::check_field(const std::string &name) const {
  auto &set = fIelds.get<FieldName_mi_tag>();
  auto miit = set.find(name);
  if (miit == set.end())
    return false;
  return true;
}

const Field *Core::get_field_structure(const std::string &name) {
  auto &set = fIelds.get<FieldName_mi_tag>();
  auto miit = set.find(name);
  if (miit == set.end()) {
    throw MoFEMException(
        MOFEM_NOT_FOUND,
        std::string("field < " + name +
                    " > not in database (top tip: check spelling)")
            .c_str());
  }
  return miit->get();
}

MoFEMErrorCode Core::get_field_entities_by_dimension(const std::string name,
                                                     int dim,
                                                     Range &ents) const {

  MoFEMFunctionBegin;
  EntityHandle meshset = get_field_meshset(name);
  CHKERR get_moab().get_entities_by_dimension(meshset, dim, ents, true);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::get_field_entities_by_type(const std::string name,
                                                EntityType type,
                                                Range &ents) const {
  MoFEMFunctionBegin;
  EntityHandle meshset = get_field_meshset(name);
  CHKERR get_moab().get_entities_by_type(meshset, type, ents, true);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::get_field_entities_by_handle(const std::string name,
                                                  Range &ents) const {
  MoFEMFunctionBegin;
  EntityHandle meshset = get_field_meshset(name);
  CHKERR get_moab().get_entities_by_handle(meshset, ents, true);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_field(const std::string &name, const FieldSpace space,
                               const FieldApproximationBase base,
                               const FieldCoefficientsNumber nb_of_coefficients,
                               const TagType tag_type, const enum MoFEMTypes bh,
                               int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  *buildMoFEM = 0;
  auto fit = fIelds.get<FieldName_mi_tag>().find(name);
  if (fit != fIelds.get<FieldName_mi_tag>().end()) {
    if (bh == MF_EXCL) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
               "field is <%s> in database", name.c_str());
    }
  } else {
    EntityHandle meshset;
    CHKERR get_moab().create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                                     meshset);
    // id
    BitFieldId id = getFieldShift();
    CHKERR get_moab().tag_set_data(th_FieldId, &meshset, 1, &id);
    // space
    CHKERR get_moab().tag_set_data(th_FieldSpace, &meshset, 1, &space);
    // base
    CHKERR get_moab().tag_set_data(th_FieldBase, &meshset, 1, &base);

    // name
    void const *tag_data[] = {name.c_str()};
    int tag_sizes[1];
    tag_sizes[0] = name.size();
    CHKERR get_moab().tag_set_by_ptr(th_FieldName, &meshset, 1, tag_data,
                                     tag_sizes);
    // name data prefix
    const std::string name_data_prefix("_App_Data");
    void const *tag_prefix_data[] = {name_data_prefix.c_str()};
    int tag_prefix_sizes[1];
    tag_prefix_sizes[0] = name_data_prefix.size();
    CHKERR get_moab().tag_set_by_ptr(th_FieldName_DataNamePrefix, &meshset, 1,
                                     tag_prefix_data, tag_prefix_sizes);
    Tag th_app_order, th_field_data, th_rank;
    // data
    std::string tag_data_name = name_data_prefix + name;
    const int def_len = 0;
    CHKERR get_moab().tag_get_handle(
        tag_data_name.c_str(), def_len, MB_TYPE_OPAQUE, th_field_data,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    // order
    ApproximationOrder def_ApproximationOrder = -1;
    const std::string Tag_ApproximationOrder_name = "_App_Order_" + name;
    CHKERR get_moab().tag_get_handle(
        Tag_ApproximationOrder_name.c_str(), sizeof(ApproximationOrder),
        MB_TYPE_OPAQUE, th_app_order, MB_TAG_CREAT | MB_TAG_BYTES | tag_type,
        &def_ApproximationOrder);
    // rank
    int def_rank = 1;
    const std::string tag_rank_name = "_Field_Rank_" + name;
    CHKERR get_moab().tag_get_handle(
        tag_rank_name.c_str(), sizeof(FieldCoefficientsNumber), MB_TYPE_OPAQUE,
        th_rank, MB_TAG_CREAT | MB_TAG_BYTES | tag_type, &def_rank);
    CHKERR get_moab().tag_set_data(th_rank, &meshset, 1, &nb_of_coefficients);

    // add meshset
    std::pair<Field_multiIndex::iterator, bool> p;
    CoordSystemsManager *cs_manger_ptr;
    CHKERR getInterface(cs_manger_ptr);
    boost::shared_ptr<CoordSys> undefined_cs_ptr;
    CHKERR cs_manger_ptr->getCoordSysPtr("UNDEFINED", undefined_cs_ptr);
    int sys_name_size[1];
    sys_name_size[0] = undefined_cs_ptr->getName().size();
    void const *sys_name[] = {&*undefined_cs_ptr->getNameRef().begin()};
    CHKERR get_moab().tag_set_by_ptr(cs_manger_ptr->get_th_CoordSysName(),
                                     &meshset, 1, sys_name, sys_name_size);
    EntityHandle coord_sys_id = undefined_cs_ptr->getMeshset();
    CHKERR get_moab().add_entities(coord_sys_id, &meshset, 1);
    p = fIelds.insert(
        boost::make_shared<Field>(moab, meshset, undefined_cs_ptr));
    if (bh == MF_EXCL) {
      if (!p.second)
        SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
                 "field not inserted %s (top tip, it could be already there)",
                 Field(moab, meshset, undefined_cs_ptr).getName().c_str());
    }
    if (verb > QUIET) {
      std::ostringstream ss;
      ss << "add: " << **p.first << std::endl;
      PetscPrintf(cOmm, ss.str().c_str());
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::addEntsToFieldByDim(const Range &ents, const int dim,
                                         const std::string &name, int verb) {

  *buildMoFEM = 0;
  EntityHandle idm = no_handle;
  if (verb == -1)
    verb = verbose;
  MoFEMFunctionBegin;
  idm = get_field_meshset(name);
  FieldSpace space;
  CHKERR get_moab().tag_get_data(th_FieldSpace, &idm, 1, &space);
  std::vector<int> nb_ents_on_dim(3, 0);
  switch (space) {
  case L2:
    CHKERR get_moab().add_entities(idm, ents);
    if (verb >= VERY_VERBOSE) {
      std::ostringstream ss;
      ss << "add entities to field " << name;
      ss << " nb. add ents " << ents.size();
      ss << std::endl;
      PetscSynchronizedPrintf(cOmm, ss.str().c_str());
    }
    break;
  case H1:
    CHKERR get_moab().add_entities(idm, ents);
    for (int dd = 0; dd != dim; ++dd) {
      Range adj_ents;
      CHKERR get_moab().get_adjacencies(ents, dd, false, adj_ents,
                                        moab::Interface::UNION);
      if (dd == 0) {
        Range topo_nodes;
        CHKERR get_moab().get_connectivity(ents, topo_nodes, true);
        Range mid_nodes;
        CHKERR get_moab().get_connectivity(ents, mid_nodes, false);
        mid_nodes = subtract(mid_nodes, topo_nodes);
        adj_ents = subtract(adj_ents, mid_nodes);
      }
      CHKERR get_moab().add_entities(idm, adj_ents);
      nb_ents_on_dim[dd] = adj_ents.size();
    }
    break;
  case HCURL:
    CHKERR get_moab().add_entities(idm, ents);
    for (int dd = 1; dd != dim; ++dd) {
      Range adj_ents;
      CHKERR get_moab().get_adjacencies(ents, dd, false, adj_ents,
                                        moab::Interface::UNION);
      CHKERR get_moab().add_entities(idm, adj_ents);
      nb_ents_on_dim[dd] = adj_ents.size();
    }
    break;
  case HDIV:
    CHKERR get_moab().add_entities(idm, ents);
    if (dim > 2) {
      Range adj_ents;
      CHKERR get_moab().get_adjacencies(ents, 2, false, adj_ents,
                                        moab::Interface::UNION);
      CHKERR get_moab().add_entities(idm, adj_ents);
      nb_ents_on_dim[2] = adj_ents.size();
    }
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "sorry, unknown space added to entity");
  }
  if (verb >= VERY_VERBOSE) {
    std::ostringstream ss;
    ss << "add entities to field " << name;
    ss << " nb. add ents " << ents.size();
    ss << " nb. add faces " << nb_ents_on_dim[2];
    ss << " nb. add edges " << nb_ents_on_dim[1];
    ss << " nb. add nodes " << nb_ents_on_dim[0];
    ss << std::endl;
    PetscSynchronizedPrintf(cOmm, ss.str().c_str());
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_ents_to_field_by_dim(const Range &ents, const int dim,
                                              const std::string &name,
                                              int verb) {
  Range ents_dim = ents.subset_by_dimension(dim);
  return addEntsToFieldByDim(ents_dim, dim, name, verb);
}

MoFEMErrorCode Core::add_ents_to_field_by_type(const Range &ents,
                                               const EntityType type,
                                               const std::string &name,
                                               int verb) {
  MoFEMFunctionBegin;
  Range ents_type = ents.subset_by_type(type);
  if (!ents_type.empty()) {
    const int dim = get_moab().dimension_from_handle(ents_type[0]);
    CHKERR addEntsToFieldByDim(ents_type, dim, name, verb);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_ents_to_field_by_dim(const EntityHandle meshset,
                                              const int dim,
                                              const std::string &name,
                                              const bool recursive, int verb) {
  MoFEMFunctionBeginHot;
  Range ents;
  CHKERR get_moab().get_entities_by_dimension(meshset, dim, ents, recursive);
  CHKERR addEntsToFieldByDim(ents, dim, name, verb);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::add_ents_to_field_by_type(const EntityHandle meshset,
                                               const EntityType type,
                                               const std::string &name,
                                               const bool recursive, int verb) {
  MoFEMFunctionBegin;
  Range ents;
  CHKERR get_moab().get_entities_by_type(meshset, type, ents, recursive);
  if (!ents.empty()) {
    const int dim = get_moab().dimension_from_handle(ents[0]);
    CHKERR addEntsToFieldByDim(ents, dim, name, verb);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::set_field_order(const Range &ents, const BitFieldId id,
                                     const ApproximationOrder order, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  *buildMoFEM = 0;

  // check field & meshset
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set_id.find(id);
  if (miit == set_id.end())
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "no filed found");

  EntityHandle idm = get_field_meshset(id);
  // intersection with field meshset
  Range ents_of_id_meshset;
  CHKERR get_moab().get_entities_by_handle(idm, ents_of_id_meshset, false);
  Range field_ents = intersect(ents, ents_of_id_meshset);
  if (verb > VERBOSE) {
    PetscSynchronizedPrintf(
        cOmm, "nb. of ents for order change in the field <%s> %d\n",
        miit->get()->getName().c_str(), field_ents.size());
  }

  // ent view by field id (in set all MoabEnts has the same FieldId)
  auto &set = entsFields.get<FieldName_mi_tag>();
  auto eiit = set.lower_bound(miit->get()->getNameRef());
  FieldEntity_multiIndex_ent_view ents_id_view;
  if (eiit != set.end()) {
    auto hi_eiit = set.upper_bound(miit->get()->getNameRef());
    std::copy(eiit, hi_eiit, std::back_inserter(ents_id_view));
  }
  if (verb > VERBOSE) {
    PetscSynchronizedPrintf(
        cOmm, "nb. of ents in the multi index field <%s> %d\n",
        miit->get()->getName().c_str(), ents_id_view.size());
  }

  // loop over ents
  int nb_ents_set_order_up = 0;
  int nb_ents_set_order_down = 0;
  int nb_ents_set_order_new = 0;

  FieldEntity_change_order modify_order(order);

  for (auto pit = field_ents.const_pair_begin();
       pit != field_ents.const_pair_end(); pit++) {
    EntityHandle first = pit->first;
    EntityHandle second = pit->second;

    // Sanity check
    switch ((*miit)->getSpace()) {
    case H1:
      if (get_moab().type_from_handle(first) == MBVERTEX) {
        if (order >= 0 && order != 1) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "approximation order for H1 space and vertex different than "
                  "1 makes not sense");
        }
      }
      break;
    case HCURL:
      if (get_moab().type_from_handle(first) == MBVERTEX) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "HDIV space on vertices makes no sense");
      }
      break;
    case HDIV:
      if (get_moab().type_from_handle(first) == MBVERTEX) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "HDIV space on vertices makes no sense");
      }
      if (get_moab().type_from_handle(first) == MBEDGE) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "HDIV space on edges makes no sense");
      }
      break;
    default:
      break;
    }

    // Entity is in database, change order only if needed
    Range ents_in_database;
    auto vit = ents_id_view.get<1>().lower_bound(first);
    auto hi_vit = ents_id_view.get<1>().upper_bound(second);
    for (; vit != hi_vit; ++vit) {
      ents_in_database.insert(vit->get()->getEnt());
      if (order >= 0) {
        // entity is in database and order is changed or reset
        const ApproximationOrder old_approximation_order =
            (*vit)->getMaxOrder();
        if (old_approximation_order == order)
          continue;
        FieldEntity_multiIndex::iterator miit =
            entsFields.get<Unique_mi_tag>().find((*vit)->getGlobalUniqueId());

        if ((*miit)->getMaxOrder() < order)
          nb_ents_set_order_up++;
        if ((*miit)->getMaxOrder() > order)
          nb_ents_set_order_down++;

        // set dofs inactive if order is reduced, and set new order to entity
        // if order is increased (note that dofs are not build if order is
        // increased)

        DofEntityByNameAndEnt &dofs_by_name =
            dofsField.get<Composite_Name_And_Ent_mi_tag>();
        DofEntityByNameAndEnt::iterator dit = dofs_by_name.lower_bound(
            boost::make_tuple((*miit)->getNameRef(), (*miit)->getEnt()));
        if (dit != dofs_by_name.end()) {
          DofEntityByNameAndEnt::iterator hi_dit = dofs_by_name.upper_bound(
              boost::make_tuple((*miit)->getNameRef(), (*miit)->getEnt()));
          for (; dit != hi_dit; dit++) {
            if ((*dit)->getDofOrder() <= order)
              continue;
            bool success = dofsField.modify(dofsField.project<0>(dit),
                                            DofEntity_active_change(false));
            if (!success)
              SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                      "modification unsuccessful");
          }
        }
        bool success =
            entsFields.modify(entsFields.project<0>(miit), modify_order);
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
      }
    }

    Range new_ents = subtract(Range(first, second), ents_in_database);
    for (Range::const_pair_iterator pit = new_ents.const_pair_begin();
         pit != new_ents.const_pair_end(); ++pit) {
      EntityHandle first = pit->first;
      EntityHandle second = pit->second;

      // reserve memory for field  dofs
      boost::shared_ptr<std::vector<FieldEntity>> ents_array =
          boost::make_shared<std::vector<FieldEntity>>(
              std::vector<FieldEntity>());

      // Add sequence to field data structure. Note that entities are allocated
      // once into vector. This vector is passed into sequence as a weak_ptr.
      // Vector is destroyed at the point last entity inside that vector is
      // destroyed.
      miit->get()->getEntSequenceContainer()->push_back(ents_array);
      ents_array->reserve(second - first + 1);

      // Entity is not in database and order is changed or reset
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent,
          hi_miit_ref_ent;
      miit_ref_ent = refinedEntities.get<Ent_mi_tag>().lower_bound(first);
      hi_miit_ref_ent = refinedEntities.get<Ent_mi_tag>().upper_bound(second);
      Range ents_in_ref_ent;
      for (; miit_ref_ent != hi_miit_ref_ent; ++miit_ref_ent) {
        const EntityHandle ent = miit_ref_ent->get()->getRefEnt();
        ents_in_ref_ent.insert(ent);
        ents_array->emplace_back(*miit, *miit_ref_ent);
        if (order >= 0) {
          modify_order(&(ents_array->back()));
        }
        nb_ents_set_order_new++;
      }

      Range ents_not_in_database =
          subtract(Range(first, second), ents_in_ref_ent);
      for (Range::iterator eit = ents_not_in_database.begin();
           eit != ents_not_in_database.end(); ++eit) {
        RefEntity ref_ent(basicEntityDataPtr, *eit);
        // FIXME: need some consistent policy in that case
        if (ref_ent.getBitRefLevel().any()) {
          std::cerr << ref_ent << std::endl;
          SETERRQ(
              PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Try to add entities which are not seeded or added to database");
        }
      }

      // Add entities to database
      auto hint = entsFields.end();
      for (auto &v : *ents_array) {
        hint = entsFields.emplace_hint(hint, ents_array, &v);
      }
    }
  }

  if (verb >= VERY_VERBOSE) {
    PetscSynchronizedPrintf(cOmm,
                            "nb. of entities in field <%s> for which order was "
                            "increased %d (order %d)\n",
                            miit->get()->getName().c_str(),
                            nb_ents_set_order_up, order);
    PetscSynchronizedPrintf(cOmm,
                            "nb. of entities in field <%s> for which order was "
                            "reduced %d (order %d)\n",
                            miit->get()->getName().c_str(),
                            nb_ents_set_order_down, order);
    PetscSynchronizedPrintf(
        cOmm,
        "nb. of entities in field <%s> for which order set %d (order %d)\n",
        miit->get()->getName().c_str(), nb_ents_set_order_new, order);
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  }

  MoFEMFunctionReturn(0);
}
MoFEMErrorCode Core::set_field_order(const EntityHandle meshset,
                                     const EntityType type, const BitFieldId id,
                                     const ApproximationOrder order, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  CHKERR get_moab().get_entities_by_type(meshset, type, ents);
  if (verb > VERBOSE) {
    PetscSynchronizedPrintf(cOmm, "nb. of ents for order change %d\n",
                            ents.size());
  }
  CHKERR set_field_order(ents, id, order, verb);
  if (verb > VERBOSE) {
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  }
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode Core::set_field_order(const EntityHandle meshset,
                                     const EntityType type,
                                     const std::string &name,
                                     const ApproximationOrder order, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  *buildMoFEM = 0;
  CHKERR set_field_order(meshset, type, getBitFieldId(name), order, verb);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode Core::set_field_order(const Range &ents, const std::string &name,
                                     const ApproximationOrder order, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  *buildMoFEM = 0;
  CHKERR set_field_order(ents, getBitFieldId(name), order, verb);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode Core::set_field_order_by_entity_type_and_bit_ref(
    const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
    const BitFieldId id, const ApproximationOrder order, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByTypeAndRefLevel(bit, mask, type,
                                                           ents, verb);
  CHKERR set_field_order(ents, id, order, verb);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::set_field_order_by_entity_type_and_bit_ref(
    const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
    const std::string &name, const ApproximationOrder order, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByTypeAndRefLevel(bit, mask, type,
                                                           ents, verb);
  CHKERR set_field_order(ents, getBitFieldId(name), order, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::buildFieldForNoField(const BitFieldId id,
                           std::map<EntityType, int> &dof_counter, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  // field it
  auto &set_id = fIelds.get<BitFieldId_mi_tag>();
  // find fiels
  auto miit = set_id.find(id);
  if (miit == set_id.end()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "field not found");
  }

  // ents in the field meshset
  Range ents_of_id_meshset;
  CHKERR get_moab().get_entities_by_handle((*miit)->meshSet, ents_of_id_meshset,
                                           false);
  if (verb > VERY_NOISY) {
    PetscSynchronizedPrintf(cOmm, "ents in field %s meshset %d\n",
                            (*miit)->getName().c_str(),
                            ents_of_id_meshset.size());
  }
  for (Range::iterator eit = ents_of_id_meshset.begin();
       eit != ents_of_id_meshset.end(); eit++) {
    // serch if field meshset is in database
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent;
    miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if (miit_ref_ent == refinedEntities.get<Ent_mi_tag>().end()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Entity is not in MoFEM database, entities in field meshset need "
              "to be seeded (i.e. bit ref level add to them)");
    }
    std::pair<FieldEntity_multiIndex::iterator, bool> e_miit;
    // create database entity
    e_miit = entsFields.insert(
        boost::make_shared<FieldEntity>(*miit, *miit_ref_ent));
    // this is nor real field in space (set order to zero)
    bool success = entsFields.modify(e_miit.first, FieldEntity_change_order(0));
    if (!success)
      SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
              "modification unsuccessful");
    FieldCoefficientsNumber rank = 0;
    // create dofs on this entity (nb. of dofs is equal to rank)
    for (; rank < (*e_miit.first)->getNbOfCoeffs(); rank++) {
      std::pair<DofEntity_multiIndex::iterator, bool> d_miit;
      // check if dof is in darabase
      d_miit.first = dofsField.project<0>(dofsField.get<Unique_mi_tag>().find(
          DofEntity::getGlobalUniqueIdCalculate(rank, *(e_miit.first))));
      // if dof is not in database
      if (d_miit.first == dofsField.end()) {
        // insert dof
        d_miit = dofsField.insert(
            boost::make_shared<DofEntity>(*(e_miit.first), 0, rank, rank));
        if (d_miit.second) {
          dof_counter[MBENTITYSET]++; // Count entities in the meshset
        }
        bool success =
            dofsField.modify(d_miit.first, DofEntity_active_change(true));
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "modification unsuccessful");
      }
    }
  }
  if (verb > VERY_VERBOSE) {
    typedef DofEntity_multiIndex::index<FieldName_mi_tag>::type DofsByName;
    DofsByName &set = dofsField.get<FieldName_mi_tag>();
    DofsByName::iterator miit2 = set.lower_bound(miit->get()->getNameRef());
    DofsByName::iterator hi_miit2 = set.upper_bound(miit->get()->getNameRef());
    assert(miit2 != hi_miit2);
    for (; miit2 != hi_miit2; miit2++) {
      std::ostringstream ss;
      ss << *miit2 << std::endl;
      ;
      PetscSynchronizedPrintf(cOmm, ss.str().c_str());
    }
    PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::buildFieldForL2H1HcurlHdiv(
    const BitFieldId id, std::map<EntityType, int> &dof_counter,
    std::map<EntityType, int> &inactive_dof_counter, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  // Find field
  auto &set_id = fIelds.get<BitFieldId_mi_tag>();
  auto field_it = set_id.find(id);
  if (field_it == set_id.end()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "Field not found");
  }
  const int rank = field_it->get()->getNbOfCoeffs();
  const boost::string_ref &field_name = field_it->get()->getNameRef();

  // Ents in the field meshset
  Range ents_of_id_meshset;
  CHKERR get_moab().get_entities_by_handle((*field_it)->meshSet,
                                           ents_of_id_meshset, false);
  if (verb > VERY_NOISY) {
    PetscSynchronizedPrintf(PETSC_COMM_SELF, "Ents in field %s meshset %d\n",
                            (*field_it)->getName().c_str(),
                            ents_of_id_meshset.size());
  }

  for (auto p_eit = ents_of_id_meshset.pair_begin();
       p_eit != ents_of_id_meshset.pair_end(); ++p_eit) {

    const EntityHandle first = p_eit->first;
    const EntityHandle second = p_eit->second;

    auto feit = entsFields.get<Composite_Name_And_Ent_mi_tag>().lower_bound(
        boost::make_tuple(field_name, first));
    if (feit == entsFields.get<Composite_Name_And_Ent_mi_tag>().end())
      continue;
    auto hi_feit = entsFields.get<Composite_Name_And_Ent_mi_tag>().upper_bound(
        boost::make_tuple(field_name, second));

    // If there are DOFs in that range is more pragmatic to remove them rather
    // than to find sub-ranges or make them inactive
    auto dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(
        boost::make_tuple(field_name, first));
    auto hi_dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(
        boost::make_tuple(field_name, second));
    dofsField.get<Composite_Name_And_Ent_mi_tag>().erase(dit, hi_dit);

    // Add vertices DOFs by bulk
    boost::shared_ptr<std::vector<DofEntity>> dofs_array =
        boost::make_shared<std::vector<DofEntity>>(std::vector<DofEntity>());
    // Add Sequence of DOFs to sequence container as weak_ptr
    int nb_dofs_on_ents = 0;
    for (auto tmp_feit = feit; tmp_feit != hi_feit; ++tmp_feit) {
      nb_dofs_on_ents += rank * tmp_feit->get()->getOrderNbDofs(
                                    tmp_feit->get()->getMaxOrder());
    }
    // Add Sequence of DOFs to sequence container as weak_ptr
    dofs_array->reserve(nb_dofs_on_ents);
    for (; feit != hi_feit; ++feit) {
      // Create dofs instances and shared pointers
      int DD = 0;
      // Loop orders (loop until max entity order is set)
      for (int oo = 0; oo <= feit->get()->getMaxOrder(); ++oo) {
        // Loop nb. dofs at order oo
        for (int dd = 0; dd < feit->get()->getOrderNbDofsDiff(oo); ++dd) {
          // Loop rank
          for (int rr = 0; rr < rank; ++rr, ++DD) {
            dofs_array->emplace_back(*feit, oo, rr, DD, true);
            ++dof_counter[feit->get()->getEntType()];
          }
        }
      }
      if (DD > feit->get()->getNbDofsOnEnt()) {
        std::ostringstream ss;
        ss << "rank " << rAnk << " ";
        ss << **feit << std::endl;
        SETERRQ3(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Expected number of DOFs on entity not equal to number added "
                 "to database (DD = %d != %d = "
                 "feit->get()->getNbDofsOnEnt())\n"
                 "%s",
                 DD, feit->get()->getNbDofsOnEnt(), ss.str().c_str());
      }
    }
    // Insert into Multi-Index container
    int dofs_field_size0 = dofsField.size();
    auto hint = dofsField.end();
    for (auto &v : *dofs_array) {
      hint = dofsField.emplace_hint(hint, dofs_array, &v);
    }
    field_it->get()->getDofSequenceContainer()->push_back(dofs_array);
    if (PetscUnlikely(static_cast<int>(dofs_array.use_count()) !=
                      static_cast<int>(dofs_array->size() + 1))) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong use count %d != %d", dofs_array.use_count(),
               dofs_array->size() + 1);
    }
    if (dofs_field_size0 + dofs_array->size() != dofsField.size()) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong number of inserted DOFs %d != %d", dofs_array->size(),
               dofsField.size() - dofs_field_size0);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::buildField(const boost::shared_ptr<Field> &field,
                                int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  if (verb > QUIET) {
    PetscSynchronizedPrintf(cOmm, "Build Field %s (rank %d)\n",
                            field->getName().c_str(), rAnk);
  }
  std::map<EntityType, int> dof_counter;
  std::map<EntityType, int> inactive_dof_counter;

  switch (field->getSpace()) {
  case NOFIELD:
    CHKERR buildFieldForNoField(field->getId(), dof_counter, verb);
    break;
  case L2:
  case H1:
  case HCURL:
  case HDIV:
    CHKERR buildFieldForL2H1HcurlHdiv(field->getId(), dof_counter,
                                      inactive_dof_counter, verb);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }

  if (verb > QUIET) {
    int nb_added_dofs = 0;
    int nb_inactive_added_dofs = 0;
    for (auto const &it : dof_counter) {
      switch (it.first) {
      case MBVERTEX:
        PetscSynchronizedPrintf(cOmm,
                                "nb added dofs (vertices) %d (inactive %d)\n",
                                it.second, inactive_dof_counter[it.first]);
        break;
      case MBEDGE:
        PetscSynchronizedPrintf(cOmm,
                                "nb added dofs (edges) %d (inactive %d)\n",
                                it.second, inactive_dof_counter[it.first]);
        break;
      case MBTRI:
        PetscSynchronizedPrintf(cOmm,
                                "nb added dofs (triangles) %d (inactive %d)\n",
                                it.second, inactive_dof_counter[it.first]);
        break;
      case MBQUAD:
        PetscSynchronizedPrintf(cOmm,
                                "nb added dofs (quads) %d (inactive %d)\n",
                                it.second, inactive_dof_counter[it.first]);
        break;
      case MBTET:
        PetscSynchronizedPrintf(cOmm, "nb added dofs (tets) %d (inactive %d)\n",
                                it.second, inactive_dof_counter[it.first]);
        break;
      case MBPRISM:
        PetscSynchronizedPrintf(cOmm,
                                "nb added dofs (prisms) %d (inactive %d)\n",
                                it.second, inactive_dof_counter[it.first]);
        break;
      case MBENTITYSET:
        PetscSynchronizedPrintf(cOmm,
                                "nb added dofs (meshsets) %d (inactive %d)\n",
                                it.second, inactive_dof_counter[it.first]);
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
      }
      nb_added_dofs += it.second;
      nb_inactive_added_dofs += inactive_dof_counter[it.first];
    }
    if (verb > QUIET) {
      PetscSynchronizedPrintf(cOmm,
                              "nb added dofs %d (number of inactive dofs %d)\n",
                              nb_added_dofs, nb_inactive_added_dofs);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_field(const std::string field_name, int verb) {
  MoFEMFunctionBegin;
  auto miit = fIelds.get<FieldName_mi_tag>().find(field_name);
  if (miit == fIelds.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "Field < %s > not found",
             field_name.c_str());
  }
  CHKERR buildField((*miit), verb);
  PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_fields(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  auto &set_id = fIelds.get<BitFieldId_mi_tag>();
  for (auto miit = set_id.begin(); miit != set_id.end(); miit++) {
    CHKERR buildField((*miit), verb);
  }
  *buildMoFEM = 1 << 0;
  if (verb > QUIET) {
    PetscSynchronizedPrintf(cOmm, "Nb. dofs %u\n", dofsField.size());
  }
  PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  MoFEMFunctionReturn(0);
  // return 0;
}

MoFEMErrorCode
Core::list_dofs_by_field_name(const std::string &field_name) const {
  MoFEMFunctionBeginHot;
  auto dit = dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
  auto hi_dit = dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
  for (; dit != hi_dit; dit++) {
    std::ostringstream ss;
    ss << "rank " << rAnk << " ";
    ss << *dit << std::endl;
    PetscSynchronizedPrintf(cOmm, ss.str().c_str());
  }
  PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode Core::list_fields() const {
  MoFEMFunctionBeginHot;
  auto &set_id = fIelds.get<BitFieldId_mi_tag>();
  for (auto &miit : set_id) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscSynchronizedPrintf(cOmm, ss.str().c_str());
  }
  PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::list_adjacencies() const {
  MoFEMFunctionBeginHot;
  for (auto &miit : entFEAdjacencies) {
    std::ostringstream ss;
    ss << miit << std::endl;
    PetscSynchronizedPrintf(cOmm, ss.str().c_str());
  }
  PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::get_problem_finite_elements_entities(const std::string &problem_name,
                                           const std::string &fe_name,
                                           const EntityHandle meshset) {
  MoFEMFunctionBegin;
  auto &prb = pRoblems.get<Problem_mi_tag>();
  auto p_miit = prb.find(problem_name);
  if (p_miit == prb.end())
    SETERRQ1(PETSC_COMM_SELF, 1, "no such problem like < %s >",
             problem_name.c_str());
  auto &numered_finite_elements =
      const_cast<NumeredEntFiniteElement_multiIndex &>(
          p_miit->numeredFiniteElements);
  auto miit =
      numered_finite_elements.get<FiniteElement_name_mi_tag>().lower_bound(
          fe_name);
  auto hi_miit =
      numered_finite_elements.get<FiniteElement_name_mi_tag>().upper_bound(
          fe_name);
  for (; miit != hi_miit; miit++) {
    EntityHandle ent = (*miit)->getEnt();
    CHKERR get_moab().add_entities(meshset, &ent, 1);
    int part = (*miit)->getPart();
    CHKERR get_moab().tag_set_data(th_Part, &ent, 1, &part);
  }
  MoFEMFunctionReturn(0);
}

FieldEntityByFieldName::iterator
Core::get_ent_field_by_name_begin(const std::string &field_name) const {
  return entsFields.get<FieldName_mi_tag>().lower_bound(field_name);
}
FieldEntityByFieldName::iterator
Core::get_ent_field_by_name_end(const std::string &field_name) const {
  return entsFields.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofEntityByFieldName::iterator
Core::get_dofs_by_name_begin(const std::string &field_name) const {
  return dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
}
DofEntityByFieldName::iterator
Core::get_dofs_by_name_end(const std::string &field_name) const {
  return dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofEntityByNameAndEnt::iterator
Core::get_dofs_by_name_and_ent_begin(const std::string &field_name,
                                     const EntityHandle ent) const {
  return dofsField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(
      boost::make_tuple(field_name, ent));
}
DofEntityByNameAndEnt::iterator
Core::get_dofs_by_name_and_ent_end(const std::string &field_name,
                                   const EntityHandle ent) const {
  return dofsField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(
      boost::make_tuple(field_name, ent));
}
DofEntityByNameAndType::iterator
Core::get_dofs_by_name_and_type_begin(const std::string &field_name,
                                      const EntityType type) const {
  return dofsField.get<Composite_Name_And_Type_mi_tag>().lower_bound(
      boost::make_tuple(field_name, type));
}
DofEntityByNameAndType::iterator
Core::get_dofs_by_name_and_type_end(const std::string &field_name,
                                    const EntityType type) const {
  return dofsField.get<Composite_Name_And_Type_mi_tag>().upper_bound(
      boost::make_tuple(field_name, type));
}
MoFEMErrorCode
Core::check_number_of_ents_in_ents_field(const std::string &name) const {
  MoFEMFunctionBeginHot;
  auto it = fIelds.get<FieldName_mi_tag>().find(name);
  if (it == fIelds.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "field not found < %s >", name.c_str());
  }
  EntityHandle meshset = (*it)->getMeshset();
  int num_entities;
  CHKERR get_moab().get_number_entities_by_handle(meshset, num_entities);
  if (entsFields.get<FieldName_mi_tag>().count((*it)->getName()) >
      (unsigned int)num_entities) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "not equal number of entities in meshset and field multiindex "
             "< %s >",
             name.c_str());
  }
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode Core::check_number_of_ents_in_ents_field() const {
  MoFEMFunctionBegin;
  for (auto &it : fIelds.get<FieldName_mi_tag>()) {
    if (it->getSpace() == NOFIELD)
      continue; // FIXME: should be treated properly, not test is just
                // skipped for this NOFIELD space
    EntityHandle meshset = it->getMeshset();
    int num_entities;
    CHKERR get_moab().get_number_entities_by_handle(meshset, num_entities);
    if (entsFields.get<FieldName_mi_tag>().count(it->getName()) >
        (unsigned int)num_entities) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "not equal number of entities in meshset and field "
               "multiindex < %s >",
               it->getName().c_str());
    }
  }
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
