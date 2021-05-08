/** \file DeleteCore.cpp
 * \brief Core interface methods for managing deletions and insertion dofs
 *
 * \note If entity/dof/finite element is cleared it mean that is erased from
 * multi-index database
 *
 * \note If entity/dof/finite element is removed is is clearded and
 * removed from filed or finite element meshset
 *
 * \note If entity if deleted is cleated, removed and deleted from MoAB
 * database.
 *
 * \todo Implement tag_field_delete, tag_field_delete_by_bit_ref to remove field
 * tags from entities. This is for entities which are as well removed, thus
 * cleared as well.
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

#define DeleteCoreFunctionBegin                                                \
  MoFEMFunctionBegin;                                                          \
  MOFEM_LOG_CHANNEL("WORLD");                                                  \
  MOFEM_LOG_FUNCTION();                                                        \
  MOFEM_LOG_TAG("WORLD", "DeleteCore");

namespace MoFEM {

MoFEMErrorCode Core::clear_inactive_dofs(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  for (DofEntity_multiIndex::iterator dit = dofsField.begin();
       dit != dofsField.end();) {
    if (!dit->get()->getActive()) {
      dit = dofsField.erase(dit);
    } else {
      ++dit;
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_dofs_fields_by_bit_ref(const BitRefLevel bit,
                                       const BitRefLevel mask, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_dofs_fields(ents, verb);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_dofs_fields(const Range ents, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;

  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first  = p_eit->first;
    EntityHandle second = p_eit->second;
    // get dofs range
    DofEntityByEnt::iterator dit, hi_dit;
    dit = dofsField.get<Ent_mi_tag>().lower_bound(first);
    if (dit == dofsField.get<Ent_mi_tag>().end())
      continue;
    hi_dit = dofsField.get<Ent_mi_tag>().upper_bound(second);
    // finally clear dofs
    dofsField.get<Ent_mi_tag>().erase(dit,hi_dit);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_dofs_fields(const std::string name,
                                       const Range ents, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;

  const auto bit_number = get_field_bit_number(name);

  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    const auto first  = p_eit->first;
    const auto second = p_eit->second;
    const auto lo_uid =
        DofEntity::getLoFieldEntityUId(bit_number, first);
    const auto hi_uid =
        DofEntity::getHiFieldEntityUId(bit_number, second);
    auto dit = dofsField.get<Unique_mi_tag>().lower_bound(lo_uid);
    auto hi_dit = dofsField.get<Unique_mi_tag>().upper_bound(hi_uid);
    dofsField.get<Unique_mi_tag>().erase(dit, hi_dit);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_ents_fields_by_bit_ref(const BitRefLevel bit,
                                       const BitRefLevel mask, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_dofs_fields(ents, verb);
  CHKERR clear_adjacencies_entities(ents, verb);
  CHKERR clear_ents_fields(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_ents_fields(const Range ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_dofs_fields(ents, verb);
  CHKERR clear_adjacencies_entities(ents, verb);
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first  = p_eit->first;
    EntityHandle second = p_eit->second;
    FieldEntity_multiIndex::index<Ent_mi_tag>::type::iterator dit, hi_dit;
    dit = entsFields.get<Ent_mi_tag>().lower_bound(first);
    hi_dit = entsFields.get<Ent_mi_tag>().upper_bound(second);
    entsFields.get<Ent_mi_tag>().erase(dit, hi_dit);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_ents_fields(const std::string name,
                                       const Range ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  const auto bit_number = get_field_bit_number(name);
  CHKERR clear_dofs_fields(name, ents, verb);
  CHKERR clear_adjacencies_entities(name, ents, verb);
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    const auto first = p_eit->first;
    const auto second = p_eit->second;
    auto dit = entsFields.get<Unique_mi_tag>().lower_bound(
        FieldEntity::getLocalUniqueIdCalculate(bit_number, first));
    auto hi_dit = entsFields.get<Unique_mi_tag>().upper_bound(
        FieldEntity::getLocalUniqueIdCalculate(bit_number, second));
    entsFields.get<Unique_mi_tag>().erase(dit, hi_dit);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_field(const std::string name,
                                            const EntityHandle meshset,
                                            const EntityType type, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR get_moab().get_entities_by_type(meshset, type, ents);
  CHKERR remove_ents_from_field(name, ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_field(const std::string name,
                                            const Range ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  EntityHandle meshset;
  meshset = get_field_meshset(name);
  CHKERR clear_ents_fields(name, ents, verb);
  CHKERR get_moab().remove_entities(meshset, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_field(const Range ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_ents_fields(ents, verb);
  for (Field_multiIndex::iterator fit = fIelds.begin(); fit != fIelds.end();
       fit++) {
    EntityHandle meshset = fit->get()->getMeshset();
    CHKERR get_moab().remove_entities(meshset, ents);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_field_by_bit_ref(const BitRefLevel bit,
                                                       const BitRefLevel mask,
                                                       int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR remove_ents_from_field(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_adjacencies_entities(const BitRefLevel bit,
                                                const BitRefLevel mask,
                                                int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_adjacencies_entities(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_adjacencies_entities(const Range ents, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); ++p_eit) {
    const EntityHandle first  = p_eit->first;
    const EntityHandle second = p_eit->second;
    FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
        Ent_mi_tag>::type::iterator ait,
        hi_ait;
    ait = entFEAdjacencies.get<Ent_mi_tag>().lower_bound(first);
    hi_ait = entFEAdjacencies.get<Ent_mi_tag>().upper_bound(second);
    entFEAdjacencies.get<Ent_mi_tag>().erase(ait, hi_ait);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_adjacencies_entities(const std::string name,
                                                const Range ents, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;

  const Field *field_ptr   = get_field_structure(name);
  int field_bit_number     = field_ptr->getBitNumber();
  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&get_moab(), basicEntityDataPtr->pcommID);

  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {

    // First and last handle
    const EntityHandle first  = p_eit->first;
    const EntityHandle second = p_eit->second;

    // Get UId
    UId first_uid =
        FieldEntity::getLocalUniqueIdCalculate(field_bit_number, first);
    UId second_uid =
        FieldEntity::getLocalUniqueIdCalculate(field_bit_number, second);

    // Find adjacencies
    auto ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(first_uid);
    auto hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(second_uid);
    entFEAdjacencies.get<Unique_mi_tag>().erase(ait, hi_ait);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_finite_elements_by_bit_ref(const BitRefLevel bit,
                                                      const BitRefLevel mask,
                                                      int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_finite_elements(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_finite_elements(const Range ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_adjacencies_finite_elements(ents, verb);
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first = p_eit->first;
    EntityHandle second = p_eit->second;
    EntFiniteElement_multiIndex::index<Ent_mi_tag>::type::iterator fit, hi_fit;
    fit = entsFiniteElements.get<Ent_mi_tag>().lower_bound(first);
    hi_fit = entsFiniteElements.get<Ent_mi_tag>().upper_bound(second);
    entsFiniteElements.get<Ent_mi_tag>().erase(fit, hi_fit);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_finite_elements(const std::string name,
                                           const Range ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_adjacencies_finite_elements(name, ents, verb);
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first  = p_eit->first;
    EntityHandle second = p_eit->second;
    EntFiniteElement_multiIndex::index<
        Composite_Name_And_Ent_mi_tag>::type::iterator fit,
        hi_fit;
    fit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().lower_bound(
        boost::make_tuple(name, first));
    hi_fit =
        entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().upper_bound(
            boost::make_tuple(name, second));
    fit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().erase(fit,
                                                                        hi_fit);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_adjacencies_finite_elements(const BitRefLevel bit,
                                                       const BitRefLevel mask,
                                                       int verb) {
  MoFEMFunctionBegin;
  Range ents;                                                       
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_adjacencies_finite_elements(ents,verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_adjacencies_finite_elements(const Range ents,
                                                       int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first = p_eit->first;
    EntityHandle second = p_eit->second;
    FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
        FEEnt_mi_tag>::type::iterator ait,
        hi_ait;
    ait = entFEAdjacencies.get<FEEnt_mi_tag>().lower_bound(first);
    hi_ait = entFEAdjacencies.get<FEEnt_mi_tag>().upper_bound(second);
    entFEAdjacencies.get<FEEnt_mi_tag>().erase(ait, hi_ait);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_adjacencies_finite_elements(const std::string name,
                                                       const Range ents,
                                                       int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;

  FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator
      it_fe = finiteElements.get<FiniteElement_name_mi_tag>().find(name);
  if (it_fe != finiteElements.get<FiniteElement_name_mi_tag>().end()) {

    const auto fe_uid = (*it_fe)->getFEUId();

    for (Range::const_pair_iterator p_eit = ents.pair_begin();
         p_eit != ents.pair_end(); p_eit++) {

      // First and last handle
      const EntityHandle first  = p_eit->first;
      const EntityHandle second = p_eit->second;

      // Get UId
      UId first_uid =
          EntFiniteElement::getLocalUniqueIdCalculate(first, fe_uid);
      UId second_uid =
          EntFiniteElement::getLocalUniqueIdCalculate(second, fe_uid);

      // Find and remove adjacencies
      FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
          FE_Unique_mi_tag>::type::iterator ait,
          hi_ait;
      ait = entFEAdjacencies.get<FE_Unique_mi_tag>().lower_bound(first_uid);
      hi_ait = entFEAdjacencies.get<FE_Unique_mi_tag>().upper_bound(second_uid);
      entFEAdjacencies.get<FE_Unique_mi_tag>().erase(ait, hi_ait);
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::remove_ents_from_finite_element(const std::string name,
                                                     const EntityHandle meshset,
                                                     const EntityType type,
                                                     int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR get_moab().get_entities_by_type(meshset, type, ents, false);
  CHKERR remove_ents_from_finite_element(name, ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_finite_element(const std::string name,
                                                     const Range ents,
                                                     int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_finite_elements(name, ents, verb);
  const EntityHandle idm = get_finite_element_meshset(name);
  CHKERR get_moab().remove_entities(idm, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_finite_element(const Range ents,
                                                     int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_finite_elements(ents, verb);
  for (FiniteElement_multiIndex::iterator fe_it = finiteElements.begin();
       fe_it != finiteElements.end(); fe_it++) {
    EntityHandle meshset = fe_it->get()->getMeshset();
    CHKERR get_moab().remove_entities(meshset, ents);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_finite_element_by_bit_ref(
    const BitRefLevel bit, const BitRefLevel mask, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;                                                       
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR remove_ents_from_finite_element(ents,verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents(const Range ents, int verb) {
   MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR remove_ents_from_finite_element(ents, verb);
  CHKERR remove_ents_from_field(ents, verb);

  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); ++p_eit) {

    RefElement_multiIndex::index<Ent_mi_tag>::type::iterator frit, hi_frit;
    frit = refinedFiniteElements.get<Ent_mi_tag>().lower_bound(p_eit->first);
    hi_frit =
        refinedFiniteElements.get<Ent_mi_tag>().upper_bound(p_eit->second);
    refinedFiniteElements.get<Ent_mi_tag>().erase(frit, hi_frit);

    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator rit, hi_rit;
    rit = refinedEntities.get<Ent_mi_tag>().lower_bound(p_eit->first);
    hi_rit = refinedEntities.get<Ent_mi_tag>().upper_bound(p_eit->second);
    refinedEntities.get<Ent_mi_tag>().erase(rit, hi_rit);

  }

  MoFEMFunctionReturn(0); 
}

MoFEMErrorCode Core::remove_ents_by_bit_ref(const BitRefLevel bit,
                                            const BitRefLevel mask, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;                                                       
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR remove_ents(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_parents_by_bit_ref(const BitRefLevel bit,
                                               const BitRefLevel mask,
                                               int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR remove_parents_by_ents(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_parents_by_ents(const Range &ents, int verb) {
  MoFEMFunctionBegin;
  
  std::vector<EntityHandle> leftovers_ents;
  leftovers_ents.reserve(ents.size());

  for (auto pit = ents.pair_begin(); pit != ents.pair_end(); ++pit) {

    EntityHandle f = pit->first;
    const EntityHandle s = pit->second;
    auto lo = refinedEntities.lower_bound(f);
    for (; f <= s; ++f) {

      auto check = [this](auto lo, auto f) {
        if (lo == refinedEntities.end())
          return false;
        if ((*lo)->getEnt() == f)
          return true;
        return false;
      };

      if (check(lo, f)) {
        bool success = refinedEntities.modify(lo, RefEntity_change_parent(0));
        if (!success)
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "Operation of removing parent unsuccessful");
        ++lo;
      } else
        leftovers_ents.emplace_back(f);
    }
  }

  if (!leftovers_ents.empty()) {
    std::vector<EntityHandle> zero_parents(leftovers_ents.size());
    CHKERR get_moab().tag_set_data(th_RefParentHandle, &leftovers_ents[0],
                                   leftovers_ents.size(),
                                   &*zero_parents.begin());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_parents_by_parents(const Range &ents, int verb) {
  MoFEMFunctionBegin;
  for (Range::iterator eit = ents.begin(); eit != ents.end(); ++eit) {
    RefEntity_multiIndex::index<Ent_Ent_mi_tag>::type::iterator it;
    while ((it = refinedEntities.get<Ent_Ent_mi_tag>().find(*eit)) !=
           refinedEntities.get<Ent_Ent_mi_tag>().end()) {
      bool success = refinedEntities.get<Ent_Ent_mi_tag>().modify(
          it, RefEntity_change_parent(0));
      if (!success) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                "Operation of removing parent unsuccessful");
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::delete_ents_by_bit_ref(const BitRefLevel bit,
                                            const BitRefLevel mask,
                                            const bool remove_parent,
                                            int verb) {
  DeleteCoreFunctionBegin;
  if (verb == -1)
    verb = verbose;

  Range ents;                                                       
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  Range ents_meshsets = ents.subset_by_type(MBENTITYSET);
  ents = subtract(ents,ents_meshsets);

  CHKERR remove_ents(ents, verb);

  // remove parent
  if (remove_parent) {
    CHKERR remove_parents_by_parents(ents);
  }

  Range meshsets;
  CHKERR get_moab().get_entities_by_type(0, MBENTITYSET, meshsets, true);
  for (Range::iterator mit = meshsets.begin(); mit != meshsets.end(); mit++) {
    CHKERR get_moab().remove_entities(*mit, ents);
  }

  rval = get_moab().delete_entities(ents);
  if(rval != MB_SUCCESS) {
    if(verb >= VERY_VERBOSE) {
      for (Range::iterator eit = ents.begin(); eit != ents.end(); ++eit) {
        try {
          RefEntity ref_ent(basicEntityDataPtr, *eit);
          MOFEM_LOG("WORLD", Sev::error)
              << "Error: " << RefEntity(basicEntityDataPtr, *eit) << " "
              << ref_ent.getBitRefLevel();
        } catch (std::exception const &ex) {
        }
      };
    }
    EntityHandle out_meshset;
    CHKERR get_moab().create_meshset(MESHSET_SET, out_meshset);
    CHKERR get_moab().add_entities(out_meshset,ents);
    CHKERR get_moab().write_file("error.vtk", "VTK", "", &out_meshset, 1);
    THROW_MESSAGE("Can not delete entities from MoAB database (see error.vtk)");
  }

  if (verb >= VERBOSE)
    MOFEM_LOG_C("WORLD", Sev::verbose, "Nb. of deleted entities %d",
                ents.size());

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::delete_finite_element(const std::string name, int verb) {
  MoFEMFunctionBegin;
  auto &fe = finiteElements.get<FiniteElement_name_mi_tag>();
  auto mit = fe.find(name);
  if (mit == fe.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
             "Finite element <%s> not found", name.c_str());
  }
  EntityHandle meshset = mit->get()->getMeshset();
  Range ents;
  CHKERR get_moab().get_entities_by_handle(meshset, ents, false);
  CHKERR remove_ents_from_finite_element(name, ents, verb);
  fe.erase(mit);
  CHKERR get_moab().delete_entities(&meshset, 1);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::delete_field(const std::string name, int verb) {
  MoFEMFunctionBegin;
  auto &f = fIelds.get<FieldName_mi_tag>();
  auto mit = f.find(name);
  if (mit == f.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
             "Finite element <%s> not found", name.c_str());
  }
  EntityHandle meshset = mit->get()->getMeshset();
  Range ents;
  CHKERR get_moab().get_entities_by_handle(meshset, ents, false);
  CHKERR remove_ents_from_field(name, ents, verb);
  CHKERR get_moab().tag_delete((*mit)->th_FieldDataVerts);
  CHKERR get_moab().tag_delete((*mit)->th_FieldData);
  CHKERR get_moab().tag_delete((*mit)->th_AppOrder);
  f.erase(mit);
  CHKERR get_moab().delete_entities(&meshset, 1);
  MoFEMFunctionReturn(0);
}
}
