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

#include <MoFEM.hpp>

#define FieldCoreFunctionBegin                                                 \
  MoFEMFunctionBegin;                                                          \
  MOFEM_LOG_CHANNEL("WORLD");                                                  \
  MOFEM_LOG_CHANNEL("SYNC");                                                   \
  MOFEM_LOG_FUNCTION();                                                        \
  MOFEM_LOG_TAG("WORLD", "FieldCore");                                         \
  MOFEM_LOG_TAG("SYNC", "FieldCore");

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

FieldBitNumber Core::get_field_bit_number(const std::string name) const {
  auto &set = fIelds.get<FieldName_mi_tag>();
  auto miit = set.find(name);
  if (miit == set.end())
    THROW_MESSAGE("field not in database (top tip: check spelling)");
  return (*miit)->getBitNumber();
}

EntityHandle Core::get_field_meshset(const BitFieldId id) const {
  auto &set = fIelds.get<BitFieldId_mi_tag>();
  auto miit = set.find(id);
  if (miit == set.end())
    THROW_MESSAGE("field not in database (top tip: check spelling)");
  return (*miit)->meshSet;
}

EntityHandle Core::get_field_meshset(const std::string name) const {
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

MoFEMErrorCode Core::addField(const std::string &name, const FieldSpace space,
                              const FieldApproximationBase base,
                              const FieldCoefficientsNumber nb_of_coefficients,
                              const TagType tag_type, const enum MoFEMTypes bh,
                              int verb) {
  MOFEM_LOG_CHANNEL("WORLD");
  MOFEM_LOG_TAG("WORLD", "FieldCore");
  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_TAG("SYNC", "FieldCore");
  MOFEM_LOG_FUNCTION();
  MoFEMFunctionBegin;

  if (verb == -1)
    verb = verbose;
  *buildMoFEM = 0;
  auto fit = fIelds.get<FieldName_mi_tag>().find(name);
  if (fit != fIelds.get<FieldName_mi_tag>().end()) {
    if (bh == MF_EXCL)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
               "field is <%s> in database", name.c_str());

  } else {

    EntityHandle meshset;
    CHKERR get_moab().create_meshset(MESHSET_SET, meshset);

    // Add field mesh set to partion meshset. In case of no elements
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

    // id
    BitFieldId id = getFieldShift();

    auto create_tags = [&]() {
      MoFEMFunctionBegin;
      CHKERR
      get_moab().tag_set_data(th_FieldId, &meshset, 1, &id);
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
      Tag th_app_order, th_field_data, th_field_data_vert, th_rank;
      // data
      std::string tag_data_name = name_data_prefix + name;
      const int def_len = 0;
      CHKERR get_moab().tag_get_handle(
          tag_data_name.c_str(), def_len, MB_TYPE_DOUBLE, th_field_data,
          MB_TAG_CREAT | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
      std::string tag_data_name_verts = name_data_prefix + name + "V";
      VectorDouble def_vert_data(nb_of_coefficients);
      def_vert_data.clear();
      CHKERR get_moab().tag_get_handle(
          tag_data_name_verts.c_str(), nb_of_coefficients, MB_TYPE_DOUBLE,
          th_field_data_vert, MB_TAG_CREAT | tag_type, &*def_vert_data.begin());
      // order
      ApproximationOrder def_ApproximationOrder = -1;
      const std::string Tag_ApproximationOrder_name = "_App_Order_" + name;
      CHKERR get_moab().tag_get_handle(
          Tag_ApproximationOrder_name.c_str(), 1, MB_TYPE_INTEGER, th_app_order,
          MB_TAG_CREAT | tag_type, &def_ApproximationOrder);
      // rank
      int def_rank = 1;
      const std::string tag_rank_name = "_Field_Rank_" + name;
      CHKERR get_moab().tag_get_handle(tag_rank_name.c_str(), 1,
                                       MB_TYPE_INTEGER, th_rank,
                                       MB_TAG_CREAT | MB_TAG_SPARSE, &def_rank);
      CHKERR get_moab().tag_set_data(th_rank, &meshset, 1, &nb_of_coefficients);

      MoFEMFunctionReturn(0);
    };

    CHKERR create_tags();

    auto p = fIelds.insert(boost::make_shared<Field>(moab, meshset));
    if (verb > QUIET) {
      MOFEM_LOG("WORLD", Sev::inform) << "Add field " << **p.first;
      MOFEM_LOG("WORLD", Sev::noisy)
          << "Field " << (*p.first)->getName() << " core value < "
          << this->getValue() << " > field value ) "
          << (*p.first)->getBitNumber() << " )";
    }

    if (!p.second)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
               "field not inserted %s (top tip, it could be already "
               "there)",
               Field(moab, meshset).getName().c_str());
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_field(const std::string &name, const FieldSpace space,
                               const FieldApproximationBase base,
                               const FieldCoefficientsNumber nb_of_coefficients,
                               const TagType tag_type, const enum MoFEMTypes bh,
                               int verb) {
  return this->addField(name, space, base, nb_of_coefficients, tag_type, bh,
                        verb);
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

MoFEMErrorCode Core::create_vertices_and_add_to_field(const std::string name,
                                                      const double coords[],
                                                      int size, int verb) {
  MoFEMFunctionBegin;

  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  Range verts;

  auto create_vertices = [&]() {
    MoFEMFunctionBegin;

    vector<double *> arrays_coord;
    EntityHandle startv = 0;
    ReadUtilIface *iface;
    CHKERR get_moab().query_interface(iface);
    CHKERR iface->get_node_coords(3, size, 0, startv, arrays_coord);
    verts.insert(startv, startv + size - 1);
    for (int n = 0; n != size; ++n)
      for (auto d : {0, 1, 2})
        arrays_coord[d][n] = coords[3 * n + d];

    MoFEMFunctionReturn(0);
  };

  auto add_verts_to_field = [&]() {
    MoFEMFunctionBegin;
    EntityHandle field_meshset = get_field_meshset(name);
    CHKERR get_moab().add_entities(field_meshset, verts);
    MoFEMFunctionReturn(0);
  };

  CHKERR create_vertices();
  CHKERR add_verts_to_field();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::setFieldOrderImpl(boost::shared_ptr<Field> field_ptr,
                                       const Range &ents,
                                       const ApproximationOrder order,
                                       int verb) {
  FieldCoreFunctionBegin;

  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;
  *buildMoFEM = 0;

  MOFEM_LOG("WORLD", Sev::noisy) << "Test field " << *field_ptr;

  const auto field_meshset = field_ptr->getMeshset();
  const auto bit_number = field_ptr->getBitNumber();

  // intersection with field meshset
  Range ents_of_id_meshset;
  CHKERR get_moab().get_entities_by_handle(field_meshset, ents_of_id_meshset,
                                           false);
  Range field_ents = intersect(ents, ents_of_id_meshset);
  if (verb > QUIET)
    MOFEM_LOG_C("SYNC", Sev::noisy,
                "change nb. of ents for order in the field <%s> %d",
                field_ptr->getName().c_str(), field_ents.size());

  // ent view by field id (in set all MoabEnts has the same FieldId)
  auto eiit = entsFields.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId(field_ptr->getBitNumber()));
  FieldEntity_multiIndex_ent_view ents_id_view;
  if (eiit != entsFields.get<Unique_mi_tag>().end()) {
    auto hi_eiit = entsFields.get<Unique_mi_tag>().upper_bound(
        FieldEntity::getHiBitNumberUId(field_ptr->getBitNumber()));
    std::copy(eiit, hi_eiit, std::back_inserter(ents_id_view));
  }

  if (verb > QUIET)
    MOFEM_LOG_C("SYNC", Sev::noisy,
                "current nb. of ents in the multi index field <%s> %d",
                field_ptr->getName().c_str(), ents_id_view.size());

  // loop over ents
  int nb_ents_set_order_up = 0;
  int nb_ents_set_order_down = 0;
  int nb_ents_set_order_new = 0;

  FieldEntity_change_order modify_order_no_size_change(order, false);
  FieldEntity_change_order modify_order_size_change(order, true);

  for (auto pit = field_ents.const_pair_begin();
       pit != field_ents.const_pair_end(); pit++) {
    EntityHandle first = pit->first;
    EntityHandle second = pit->second;

    // Sanity check
    switch (field_ptr->getSpace()) {
    case H1:
      break;
    case HCURL:
      if (get_moab().type_from_handle(first) == MBVERTEX)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Hcurl space on vertices makes no sense");

      break;
    case HDIV:
      if (get_moab().type_from_handle(first) == MBVERTEX)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Hdiv space on vertices makes no sense");

      if (get_moab().type_from_handle(first) == MBEDGE)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Hdiv space on edges makes no sense");

      break;
    default:
      break;
    }

    // Entity is in database, change order only if needed
    Range ents_in_database;
    auto vit = ents_id_view.get<1>().lower_bound(first);
    auto hi_vit = ents_id_view.get<1>().upper_bound(second);
    if (order >= 0) {
      for (; vit != hi_vit; ++vit) {
        ents_in_database.insert(vit->get()->getEnt());
        // entity is in database and order is changed or reset
        const ApproximationOrder old_approximation_order =
            (*vit)->getMaxOrder();

        if (old_approximation_order != order) {

          FieldEntity_multiIndex::iterator miit =
              entsFields.get<Unique_mi_tag>().find((*vit)->getLocalUniqueId());

          if ((*miit)->getMaxOrder() < order)
            nb_ents_set_order_up++;
          if ((*miit)->getMaxOrder() > order)
            nb_ents_set_order_down++;

          // set dofs inactive if order is reduced, and set new order to entity
          // if order is increased (note that dofs are not build if order is
          // increased)

          bool can_change_size = true;
          auto dit = dofsField.get<Unique_mi_tag>().lower_bound(
              FieldEntity::getLoLocalEntityBitNumber(bit_number,
                                                     (*miit)->getEnt()));
          if (dit != dofsField.get<Unique_mi_tag>().end()) {
            auto hi_dit = dofsField.get<Unique_mi_tag>().upper_bound(
                FieldEntity::getHiLocalEntityBitNumber(bit_number,
                                                       (*miit)->getEnt()));

            if (dit != hi_dit)
              can_change_size = false;
            for (; dit != hi_dit; dit++) {
              if ((*dit)->getDofOrder() > order) {
                bool success = dofsField.modify(dofsField.project<0>(dit),
                                                DofEntity_active_change(false));
                if (!success)
                  SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                          "modification unsuccessful");
              }
            }
          }

          bool success =
              entsFields.modify(entsFields.project<0>(miit),
                                can_change_size ? modify_order_size_change
                                                : modify_order_no_size_change);
          if (!success)
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
        }
      }
    }

    Range new_ents = subtract(Range(first, second), ents_in_database);
    for (Range::const_pair_iterator pit = new_ents.const_pair_begin();
         pit != new_ents.const_pair_end(); ++pit) {
      EntityHandle first = pit->first;
      EntityHandle second = pit->second;
      const EntityType ent_type = get_moab().type_from_handle(first);
      auto get_nb_dofs_on_order = [&](const int order) {
        return order >= 0 ? (field_ptr->getFieldOrderTable()[ent_type])(order)
                          : 0;
      };
      const int nb_dofs_on_order = get_nb_dofs_on_order(order);
      if (nb_dofs_on_order || order == -1) {

        const int field_rank = field_ptr->getNbOfCoeffs();
        const int nb_dofs = nb_dofs_on_order * field_rank;

        // Entity is not in database and order is changed or reset
        auto miit_ref_ent =
            refinedEntities.get<Ent_mi_tag>().lower_bound(first);

        auto create_tags_for_max_order = [&](const Range &ents) {
          MoFEMFunctionBegin;
          if (order >= 0) {
            std::vector<ApproximationOrder> o_vec(ents.size(), order);
            CHKERR get_moab().tag_set_data(field_ptr->th_AppOrder, ents,
                                           &*o_vec.begin());
          }
          MoFEMFunctionReturn(0);
        };

        auto create_tags_for_data = [&](const Range &ents) {
          MoFEMFunctionBegin;
          if (order >= 0) {

            if (nb_dofs > 0) {
              if (ent_type == MBVERTEX) {
                std::vector<FieldData> d_vec(nb_dofs * ents.size(), 0);
                CHKERR get_moab().tag_set_data(field_ptr->th_FieldDataVerts,
                                               ents, &*d_vec.begin());
              } else {
                std::vector<int> tag_size(ents.size(), nb_dofs);
                std::vector<FieldData> d_vec(nb_dofs, 0);
                std::vector<void const *> d_vec_ptr(ents.size(),
                                                    &*d_vec.begin());
                CHKERR get_moab().tag_set_by_ptr(field_ptr->th_FieldData, ents,
                                                 &*d_vec_ptr.begin(),
                                                 &*tag_size.begin());
              }
            }
          }

          MoFEMFunctionReturn(0);
        };

        auto get_ents_in_ref_ent = [&](auto miit_ref_ent) {
          auto hi = refinedEntities.get<Ent_mi_tag>().upper_bound(second);
          Range in;
          for (; miit_ref_ent != hi; ++miit_ref_ent)
            in.insert(miit_ref_ent->get()->getEnt());
          return in;
        };

        auto get_ents_max_order = [&](const Range &ents) {
          boost::shared_ptr<std::vector<const void *>> vec(
              new std::vector<const void *>());
          vec->resize(ents.size());
          CHKERR get_moab().tag_get_by_ptr(field_ptr->th_AppOrder, ents,
                                           &*vec->begin());
          return vec;
        };

        auto get_ents_field_data_vector_adaptor =
            [&](const Range &ents,
                boost::shared_ptr<std::vector<const void *>> &ents_max_orders) {
              // create shared pointer and reserve memory
              boost::shared_ptr<std::vector<double *>> vec(
                  new std::vector<double *>());
              vec->reserve(ents.size());

              if (order >= 0 && nb_dofs == 0) {
                // set empty vector adaptor
                for (int i = 0; i != ents.size(); ++i)
                  vec->emplace_back(nullptr);
              } else {
                moab::ErrorCode rval;
                std::vector<int> tag_size(ents.size());
                std::vector<const void *> d_vec_ptr(ents.size());

                // get tags data
                if (ent_type == MBVERTEX)
                  rval = get_moab().tag_get_by_ptr(field_ptr->th_FieldDataVerts,
                                                   ents, &*d_vec_ptr.begin(),
                                                   &*tag_size.begin());
                else
                  rval = get_moab().tag_get_by_ptr(field_ptr->th_FieldData,
                                                   ents, &*d_vec_ptr.begin(),
                                                   &*tag_size.begin());

                auto cast = [](auto p) {
                  return const_cast<FieldData *const>(
                      static_cast<const FieldData *>(p));
                };

                // some of entities has tag not set or zero dofs on entity
                if (rval == MB_SUCCESS) {
                  // all is ok, all entities has tag set
                  auto tit = d_vec_ptr.begin();
                  auto oit = ents_max_orders->begin();
                  for (auto sit = tag_size.begin(); sit != tag_size.end();
                       ++sit, ++tit, ++oit)
                    vec->emplace_back(cast(*tit));

                } else {
                  // set empty vector adaptor
                  for (int i = 0; i != ents.size(); ++i)
                    vec->emplace_back(nullptr);

                  // check order on all entities, and if for that order non zero
                  // dofs are expected get pointer to tag data and reset vector
                  // adaptor
                  auto oit = ents_max_orders->begin();
                  auto dit = vec->begin();
                  for (auto eit = ents.begin(); eit != ents.end();
                       ++eit, ++oit, ++dit) {

                    const int ent_order =
                        *static_cast<const ApproximationOrder *>(*oit);
                    const int ent_nb_dofs = get_nb_dofs_on_order(ent_order);

                    if (ent_nb_dofs) {
                      int tag_size;
                      const void *ret_val;
                      if (ent_type == MBVERTEX) {
                        rval = get_moab().tag_get_by_ptr(
                            field_ptr->th_FieldDataVerts, &*eit, 1, &ret_val,
                            &tag_size);
                      } else {
                        rval = get_moab().tag_get_by_ptr(
                            field_ptr->th_FieldData, &*eit, 1, &ret_val,
                            &tag_size);
                        if (rval != MB_SUCCESS) {

                          const int set_tag_size[] = {ent_nb_dofs * field_rank};
                          std::array<FieldData, MAX_DOFS_ON_ENTITY> set_d_vec;
                          std::fill(set_d_vec.begin(),
                                    &set_d_vec[set_tag_size[0]], 0);
                          const void *set_d_vec_ptr[] = {set_d_vec.data()};
                          CHKERR get_moab().tag_set_by_ptr(
                              field_ptr->th_FieldData, &*eit, 1, set_d_vec_ptr,
                              set_tag_size);
                          rval = get_moab().tag_get_by_ptr(
                              field_ptr->th_FieldData, &*eit, 1, &ret_val,
                              &tag_size);

                          if (rval != MB_SUCCESS) {
                            MOFEM_LOG_ATTRIBUTES("SELF",
                                                 LogManager::BitLineID |
                                                     LogManager::BitScope);
                            MOFEM_LOG("SELF", Sev::error)
                                << "Error is triggered in MOAB, field tag data "
                                   "for same reason can not be for accessed.";
                            MOFEM_LOG("SELF", Sev::error)
                                << "Set order: " << order;
                            MOFEM_LOG("SELF", Sev::error)
                                << "Nb. dofs on entity for given order: "
                                << set_tag_size[0];
                            MOFEM_LOG("SELF", Sev::error)
                                << "Entity type: "
                                << moab::CN::EntityTypeName(ent_type);
                            MOFEM_LOG("SELF", Sev::error)
                                << "Field: " << *field_ptr;
                            CHKERR rval;
                          }
                        }
                      }
                      const_cast<FieldData *&>(*dit) = cast(ret_val);
                    }
                  }
                }
              }
              return vec;
            };

        auto ents_in_ref_ent = get_ents_in_ref_ent(miit_ref_ent);

        CHKERR create_tags_for_max_order(ents_in_ref_ent);
        CHKERR create_tags_for_data(ents_in_ref_ent);
        auto ents_max_order = get_ents_max_order(ents_in_ref_ent);
        auto ent_field_data =
            get_ents_field_data_vector_adaptor(ents_in_ref_ent, ents_max_order);

        // reserve memory for field  dofs
        auto ents_array = boost::make_shared<std::vector<FieldEntity>>();
        // Add sequence to field data structure. Note that entities are
        // allocated once into vector. This vector is passed into sequence as a
        // weak_ptr. Vector is destroyed at the point last entity inside that
        // vector is destroyed.
        ents_array->reserve(second - first + 1);
        auto vit_max_order = ents_max_order->begin();
        auto vit_field_data = ent_field_data->begin();
        for (auto ent : ents_in_ref_ent) {
          ents_array->emplace_back(
              field_ptr, *miit_ref_ent,
              boost::shared_ptr<double *const>(ent_field_data,
                                               &*vit_field_data),
              boost::shared_ptr<const int>(
                  ents_max_order, static_cast<const int *>(*vit_max_order)));
          ++vit_max_order;
          ++vit_field_data;
          ++miit_ref_ent;
        }
        if (!ents_array->empty())
          if ((*ents_array)[0].getFieldRawPtr() != field_ptr.get())
            SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Get field ent poiter and field pointer do not match for "
                     "field %s",
                     field_ptr->getName().c_str());
        nb_ents_set_order_new += ents_array->size();

        // Check if any of entities in the range has bit level but is not added
        // to database. That generate data inconsistency and error.
        if (ents_in_ref_ent.size() < (second - first + 1)) {
          Range ents_not_in_database =
              subtract(Range(first, second), ents_in_ref_ent);
          std::vector<const void *> vec_bits(ents_not_in_database.size());
          CHKERR get_moab().tag_get_by_ptr(
              get_basic_entity_data_ptr()->th_RefBitLevel, ents_not_in_database,
              &*vec_bits.begin());
          auto cast = [](auto p) {
            return static_cast<const BitRefLevel *>(p);
          };
          for (auto v : vec_bits)
            if (cast(v)->any())
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "Try to add entities which are not seeded or added to "
                      "database");
        }

        // Add entities to database
        auto hint = entsFields.end();
        for (auto &v : *ents_array)
          hint = entsFields.emplace_hint(hint, ents_array, &v);
      }
    }
  }

  if (verb > QUIET) {
    MOFEM_LOG_C("SYNC", Sev::noisy,
                "nb. of entities in field <%s> for which order was "
                "increased %d (order %d)",
                field_ptr->getName().c_str(), nb_ents_set_order_up, order);
    MOFEM_LOG_C("SYNC", Sev::noisy,
                "nb. of entities in field <%s> for which order was "
                "reduced %d (order %d)",
                field_ptr->getName().c_str(), nb_ents_set_order_down, order);
    MOFEM_LOG_C(
        "SYNC", Sev::noisy,
        "nb. of entities in field <%s> for which order set %d (order %d)",
        field_ptr->getName().c_str(), nb_ents_set_order_new, order);
    MOFEM_LOG_SYNCHRONISE(cOmm);
  }

  if (verb > QUIET) {
    auto eiit = entsFields.get<Unique_mi_tag>().lower_bound(
        FieldEntity::getLoBitNumberUId(field_ptr->getBitNumber()));
    auto hi_eiit = entsFields.get<Unique_mi_tag>().upper_bound(
        FieldEntity::getHiBitNumberUId(field_ptr->getBitNumber()));
    MOFEM_LOG_C("SYNC", Sev::noisy,
                "nb. of ents in the multi index field <%s> %d",
                field_ptr->getName().c_str(), std::distance(eiit, hi_eiit));
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::set_field_order(const Range &ents, const BitFieldId id,
                                     const ApproximationOrder order, int verb) {
  MOFEM_LOG_CHANNEL("WORLD");
  MOFEM_LOG_TAG("WORLD", "FieldCore");
  MOFEM_LOG_FUNCTION();
  MoFEMFunctionBegin;

  // check field & meshset
  auto field_it = fIelds.get<BitFieldId_mi_tag>().find(id);
  if (field_it == fIelds.get<BitFieldId_mi_tag>().end())
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "no filed found");

  MOFEM_LOG("WORLD", Sev::noisy)
      << "Field " << (*field_it)->getName() << " core value < "
      << this->getValue() << " > field value ( " << (*field_it)->getBitNumber()
      << " )";

  CHKERR this->setFieldOrderImpl(*field_it, ents, order, verb);

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
  CHKERR this->set_field_order(ents, id, order, verb);
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
  CHKERR this->set_field_order(meshset, type, getBitFieldId(name), order, verb);
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode Core::set_field_order(const Range &ents, const std::string &name,
                                     const ApproximationOrder order, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  *buildMoFEM = 0;
  CHKERR this->set_field_order(ents, getBitFieldId(name), order, verb);
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
  CHKERR this->set_field_order(ents, id, order, verb);
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
  CHKERR this->set_field_order(ents, getBitFieldId(name), order, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::buildFieldForNoFieldImpl(boost::shared_ptr<Field> field_ptr,
                               std::map<EntityType, int> &dof_counter,
                               int verb) {
  FieldCoreFunctionBegin;

  const auto bit_number = field_ptr->getBitNumber();

  // ents in the field meshset
  Range ents;
  CHKERR get_moab().get_entities_by_handle(field_ptr->getMeshset(), ents,
                                           false);
  if (verb > VERBOSE)
    MOFEM_LOG_C("SYNC", Sev::noisy, "Ents in field %s meshset %d\n",
                field_ptr->getName().c_str(), ents.size());

  // ent view by field id (in set all MoabEnts has the same FieldId)
  auto eiit = entsFields.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId(field_ptr->getBitNumber()));
  FieldEntity_multiIndex_ent_view ents_id_view;
  if (eiit != entsFields.get<Unique_mi_tag>().end()) {
    auto hi_eiit = entsFields.get<Unique_mi_tag>().upper_bound(
        FieldEntity::getHiBitNumberUId(field_ptr->getBitNumber()));
    std::copy(eiit, hi_eiit, std::back_inserter(ents_id_view));
  }

  boost::shared_ptr<const int> zero_order(new const int(0));

  for (auto ent : ents) {
    // search if field meshset is in database
    auto ref_ent_it = refinedEntities.get<Ent_mi_tag>().find(ent);
    if (ref_ent_it == refinedEntities.get<Ent_mi_tag>().end())
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Entity is not in MoFEM database, entities in field meshset need "
              "to be seeded (i.e. bit ref level add to them)");

    auto add_dofs = [&](auto field_eit) {
      MoFEMFunctionBegin;
      // create dofs on this entity (nb. of dofs is equal to rank)
      for (FieldCoefficientsNumber rank = 0; rank < field_ptr->getNbOfCoeffs();
           rank++) {
        // insert dof
        auto p = dofsField.insert(
            boost::make_shared<DofEntity>(field_eit, 0, rank, rank));
        if (p.second) {
          dof_counter[MBENTITYSET]++; // Count entities in the meshset
        } else
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "Dof expected to be created");
      }
      MoFEMFunctionReturn(0);
    };

    // create database entity
    auto field_ent_it = ents_id_view.get<1>().find(ent);
    if (field_ent_it == ents_id_view.get<1>().end()) {

      auto p = entsFields.insert(

          boost::make_shared<FieldEntity>(
              field_ptr, *ref_ent_it,
              FieldEntity::makeSharedFieldDataAdaptorPtr(field_ptr,
                                                         *ref_ent_it),
              boost::shared_ptr<const int>(zero_order, zero_order.get()))

      );

      if ((*p.first)->getFieldRawPtr() != field_ptr.get())
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Get field ent poiter and field pointer do not match for "
                 "field %s",
                 field_ptr->getName().c_str());

      if (!p.second)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Entity should be created");

      CHKERR add_dofs(*(p.first));

    } else {

      // If there are DOFs in that range is more pragmatic to remove them
      // rather than to find sub-ranges or make them inactive
      auto dit = dofsField.get<Unique_mi_tag>().lower_bound(
          FieldEntity::getLoLocalEntityBitNumber(bit_number, ent));
      auto hi_dit = dofsField.get<Unique_mi_tag>().upper_bound(
          FieldEntity::getHiLocalEntityBitNumber(bit_number, ent));
      dofsField.get<Unique_mi_tag>().erase(dit, hi_dit);
      CHKERR add_dofs(*field_ent_it);
    }
  }

  if (verb > VERBOSE) {
    auto lo_dof = dofsField.get<Unique_mi_tag>().lower_bound(
        FieldEntity::getLoBitNumberUId(field_ptr->getBitNumber()));
    auto hi_dof = dofsField.get<Unique_mi_tag>().upper_bound(
        FieldEntity::getHiBitNumberUId(field_ptr->getBitNumber()));
    for (; lo_dof != hi_dof; lo_dof++)
      MOFEM_LOG("SYNC", Sev::noisy) << **lo_dof;
    MOFEM_LOG_SYNCHRONISE(cOmm);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::buildFieldForNoField(const BitFieldId id,
                           std::map<EntityType, int> &dof_counter, int verb) {
  FieldCoreFunctionBegin;

  if (verb == -1)
    verb = verbose;

  // find fields
  auto field_it = fIelds.get<BitFieldId_mi_tag>().find(id);
  if (field_it == fIelds.get<BitFieldId_mi_tag>().end())
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "Field not found");

  if (verb > QUIET)
    MOFEM_LOG("WORLD", Sev::noisy)
        << "Field " << (*field_it)->getName() << " core value < "
        << this->getValue() << " > field value () "
        << (*field_it)->getBitNumber() << " )";

  CHKERR this->buildFieldForNoFieldImpl(*field_it, dof_counter, verb);

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
  const int bit_number = field_it->get()->getBitNumber();
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
    const auto lo_uid =
        FieldEntity::getLoLocalEntityBitNumber(bit_number, first);
    const auto hi_uid =
        FieldEntity::getHiLocalEntityBitNumber(bit_number, second);

    auto feit = entsFields.get<Unique_mi_tag>().lower_bound(lo_uid);
    if (feit == entsFields.get<Unique_mi_tag>().end())
      continue;
    auto hi_feit = entsFields.get<Unique_mi_tag>().upper_bound(hi_uid);

    // If there are DOFs in that range is more pragmatic to remove them
    // rather than to find sub-ranges or make them inactive
    auto dit = dofsField.get<Unique_mi_tag>().lower_bound(lo_uid);
    auto hi_dit = dofsField.get<Unique_mi_tag>().upper_bound(hi_uid);
    dofsField.get<Unique_mi_tag>().erase(dit, hi_dit);

    // Add vertices DOFs by bulk
    boost::shared_ptr<std::vector<DofEntity>> dofs_array =
        boost::make_shared<std::vector<DofEntity>>(std::vector<DofEntity>());
    // Add Sequence of DOFs to sequence container as weak_ptr
    int nb_dofs_on_ents = 0;
    for (auto tmp_feit = feit; tmp_feit != hi_feit; ++tmp_feit) {
      nb_dofs_on_ents += rank * tmp_feit->get()->getOrderNbDofs(
                                    tmp_feit->get()->getMaxOrder());
    }

    // Reserve memory
    dofs_array->reserve(nb_dofs_on_ents);

    // Create DOFs
    for (; feit != hi_feit; ++feit) {
      // Create dofs instances and shared pointers
      int DD = 0;
      // Loop orders (loop until max entity order is set)
      for (int oo = 0; oo <= feit->get()->getMaxOrder(); ++oo) {
        // Loop nb. dofs at order oo
        for (int dd = 0; dd < feit->get()->getOrderNbDofsDiff(oo); ++dd) {
          // Loop rank
          for (int rr = 0; rr < rank; ++rr, ++DD) {
            dofs_array->emplace_back(*feit, oo, rr, DD);
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
    for (auto &v : *dofs_array)
      hint = dofsField.emplace_hint(hint, dofs_array, &v);

    // Add Sequence of DOFs to sequence container as weak_ptr
    field_it->get()->getDofSequenceContainer().push_back(dofs_array);

    // Check data consistency
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
  FieldCoreFunctionBegin;
  if (verb == -1)
    verb = verbose;
  if (verb > QUIET)
    MOFEM_LOG("SYNC", Sev::verbose) << "Build field " << field->getName();

  std::map<EntityType, int> dof_counter;
  std::map<EntityType, int> inactive_dof_counter;

  // Need to rebuild order table since number of dofs on each order when
  // field was created.
  if (field->getApproxBase() == USER_BASE)
    CHKERR field->rebuildDofsOrderMap();

  switch (field->getSpace()) {
  case NOFIELD:
    CHKERR this->buildFieldForNoField(field->getId(), dof_counter, verb);
    break;
  case L2:
  case H1:
  case HCURL:
  case HDIV:
    CHKERR this->buildFieldForL2H1HcurlHdiv(field->getId(), dof_counter,
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
        MOFEM_LOG("SYNC", Sev::verbose)
            << "Nb. of dofs (vertices) " << it.second << " (inactive "
            << inactive_dof_counter[it.first] << ")";
        break;
      case MBEDGE:
        MOFEM_LOG("SYNC", Sev::verbose)
            << "Nb. of dofs (edge) " << it.second << " (inactive "
            << inactive_dof_counter[it.first] << ")";
        break;
      case MBTRI:
        MOFEM_LOG("SYNC", Sev::verbose)
            << "Nb. of dofs (triangles) " << it.second << " (inactive "
            << inactive_dof_counter[it.first] << ")";
        break;
      case MBQUAD:
        MOFEM_LOG("SYNC", Sev::verbose)
            << "Nb. of dofs (quads) " << it.second << " (inactive "
            << inactive_dof_counter[it.first] << ")";
        break;
      case MBTET:
        MOFEM_LOG("SYNC", Sev::verbose)
            << "Nb. of dofs (tetrahedra) " << it.second << " (inactive "
            << inactive_dof_counter[it.first] << ")";
        break;
      case MBPRISM:
        MOFEM_LOG("SYNC", Sev::verbose)
            << "Nb. of dofs (prisms) " << it.second << " (inactive "
            << inactive_dof_counter[it.first] << ")";
        break;
      case MBENTITYSET:
        MOFEM_LOG("SYNC", Sev::verbose)
            << "Nb. of dofs (meshsets) " << it.second << " (inactive "
            << inactive_dof_counter[it.first] << ")";
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
      }
      nb_added_dofs += it.second;
      nb_inactive_added_dofs += inactive_dof_counter[it.first];
    }
    if (verb > QUIET) {
      MOFEM_LOG("SYNC", Sev::verbose)
          << "Nb. added dofs " << nb_added_dofs << " (number of inactive dofs "
          << nb_inactive_added_dofs << " )";
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_field(const std::string field_name, int verb) {
  FieldCoreFunctionBegin;
  auto field_it = fIelds.get<FieldName_mi_tag>().find(field_name);
  if (field_it == fIelds.get<FieldName_mi_tag>().end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "Field < %s > not found",
             field_name.c_str());

  CHKERR this->buildField(*field_it, verb);
  if (verb > QUIET)
    MOFEM_LOG_SYNCHRONISE(cOmm);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_fields(int verb) {
  FieldCoreFunctionBegin;
  if (verb == -1)
    verb = verbose;

  for (auto field : fIelds.get<BitFieldId_mi_tag>())
    CHKERR this->buildField(field, verb);

  *buildMoFEM = 1 << 0;
  if (verb > QUIET) {
    MOFEM_LOG("SYNC", Sev::inform) << "Number of dofs " << dofsField.size();
    MOFEM_LOG_SYNCHRONISE(cOmm);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::list_dofs_by_field_name(const std::string &field_name) const {
  FieldCoreFunctionBegin;
  auto dit = dofsField.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId(get_field_bit_number((field_name))));
  auto hi_dit = dofsField.get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId(get_field_bit_number(field_name)));
  MOFEM_LOG("SYNC", Sev::inform) << "List DOFs:";
  for (; dit != hi_dit; dit++)
    MOFEM_LOG("SYNC", Sev::inform) << *dit;

  MOFEM_LOG_SYNCHRONISE(cOmm);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::list_fields() const {
  FieldCoreFunctionBegin;
  MOFEM_LOG("SYNC", Sev::inform) << "List Fields:";
  for (auto &miit : fIelds.get<BitFieldId_mi_tag>())
    MOFEM_LOG("SYNC", Sev::inform) << *miit;

  MOFEM_LOG_SYNCHRONISE(cOmm);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::get_problem_finite_elements_entities(const std::string &problem_name,
                                           const std::string &fe_name,
                                           const EntityHandle meshset) {
  MoFEMFunctionBegin;
  auto &prb = pRoblems.get<Problem_mi_tag>();
  auto p_miit = prb.find(problem_name);
  if (p_miit == prb.end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "No such problem like < %s >", problem_name.c_str());
  auto miit = p_miit->numeredFiniteElementsPtr->get<FiniteElement_name_mi_tag>()
                  .lower_bound(fe_name);
  auto hi_miit = p_miit->numeredFiniteElementsPtr->get<FiniteElement_name_mi_tag>()
                     .upper_bound(fe_name);
  for (; miit != hi_miit; miit++) {
    EntityHandle ent = (*miit)->getEnt();
    CHKERR get_moab().add_entities(meshset, &ent, 1);
    const int part = (*miit)->getPart();
    CHKERR get_moab().tag_set_data(th_Part, &ent, 1, &part);
  }
  MoFEMFunctionReturn(0);
}

FieldEntityByUId::iterator
Core::get_ent_field_by_name_begin(const std::string &field_name) const {
  return entsFields.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId(get_field_bit_number(field_name)));
}
FieldEntityByUId::iterator
Core::get_ent_field_by_name_end(const std::string &field_name) const {
  return entsFields.get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId(get_field_bit_number(field_name)));
}
DofEntityByUId::iterator
Core::get_dofs_by_name_begin(const std::string &field_name) const {
  return dofsField.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId(get_field_bit_number(field_name)));
}
DofEntityByUId::iterator
Core::get_dofs_by_name_end(const std::string &field_name) const {
  return dofsField.get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId(get_field_bit_number(field_name)));
}
DofEntityByUId::iterator
Core::get_dofs_by_name_and_ent_begin(const std::string &field_name,
                                     const EntityHandle ent) const {
  return dofsField.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoLocalEntityBitNumber(get_field_bit_number(field_name),
                                             ent));
}
DofEntityByUId::iterator
Core::get_dofs_by_name_and_ent_end(const std::string &field_name,
                                   const EntityHandle ent) const {
  return dofsField.get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiLocalEntityBitNumber(get_field_bit_number(field_name),
                                             ent));
}
DofEntityByUId::iterator
Core::get_dofs_by_name_and_type_begin(const std::string &field_name,
                                      const EntityType type) const {
  return dofsField.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoLocalEntityBitNumber(get_field_bit_number(field_name),
                                             get_id_for_min_type(type)));
}
DofEntityByUId::iterator
Core::get_dofs_by_name_and_type_end(const std::string &field_name,
                                    const EntityType type) const {
  return dofsField.get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiLocalEntityBitNumber(get_field_bit_number(field_name),
                                             get_id_for_max_type(type)));
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

  auto count_field_ents = [&]() {
    auto bit_number = (*it)->getBitNumber();
    auto low_eit = entsFields.get<Unique_mi_tag>().lower_bound(
        FieldEntity::getLoBitNumberUId(bit_number));
    auto hi_eit = entsFields.get<Unique_mi_tag>().upper_bound(
        FieldEntity::getHiBitNumberUId(bit_number));
    return std::distance(low_eit, hi_eit);
  };

  if (count_field_ents() > (unsigned int)num_entities) {
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

    auto count_field_ents = [&]() {
      auto bit_number = it->getBitNumber();
      auto low_eit = entsFields.get<Unique_mi_tag>().lower_bound(
          FieldEntity::getLoBitNumberUId(bit_number));
      auto hi_eit = entsFields.get<Unique_mi_tag>().upper_bound(
          FieldEntity::getHiBitNumberUId(bit_number));
      return std::distance(low_eit, hi_eit);
    };

    if (count_field_ents() > (unsigned int)num_entities) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "not equal number of entities in meshset and field "
               "multiindex < %s >",
               it->getName().c_str());
    }
  }
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
