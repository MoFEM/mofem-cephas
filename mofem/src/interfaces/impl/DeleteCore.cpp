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

MoFEMErrorCode Core::clear_inactive_dofs(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  for (DofEntity_multiIndex::iterator dit = dofsField.begin();
       dit != dofsField.end();++dit) {
    if (!dit->get()->getActive()) {
      ents.insert(dit->get()->getEnt());
    }
  }
  CHKERR clear_dofs_fields(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_dofs_fields(const BitRefLevel &bit,
                                       const BitRefLevel &mask, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_dofs_fields(ents, verb);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_dofs_fields(const Range &ents, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first  = p_eit->first;
    EntityHandle second = p_eit->second;
    DofEntityByEnt::iterator dit, hi_dit;
    dit    = dofsField.get<Ent_mi_tag>().lower_bound(first);
    hi_dit = dofsField.get<Ent_mi_tag>().upper_bound(second);
    for (; dit != hi_dit;) {
      FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
          Unique_mi_tag>::type::iterator ait,
          hi_ait;
      ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(
          (*dit)->getFieldEntityPtr()->getGlobalUniqueId());
      hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(
          (*dit)->getFieldEntityPtr()->getGlobalUniqueId());
      for (; ait != hi_ait; ait++) {
        boost::shared_ptr<EntFiniteElement> ent_fe_ptr;
        ent_fe_ptr = ait->entFePtr;
        ent_fe_ptr->row_dof_view->erase((*dit)->getGlobalUniqueId());
        if (ent_fe_ptr->row_dof_view != ent_fe_ptr->col_dof_view) {
          ent_fe_ptr->col_dof_view->erase((*dit)->getGlobalUniqueId());
        }
        ent_fe_ptr->data_dofs->get<Unique_mi_tag>().erase(
            (*dit)->getGlobalUniqueId());
      }
      dit = dofsField.get<Ent_mi_tag>().erase(dit);
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_dofs_fields(const std::string &name,
                                       const Range &ents, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first  = p_eit->first;
    EntityHandle second = p_eit->second;
    DofEntityByNameAndEnt::iterator dit, hi_dit;
    dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(
        boost::make_tuple(name, first));
    hi_dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(
        boost::make_tuple(name, second));
    for (; dit != hi_dit;) {
      FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
          Unique_mi_tag>::type::iterator ait,
          hi_ait;
      ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(
          (*dit)->getFieldEntityPtr()->getGlobalUniqueId());
      hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(
          (*dit)->getFieldEntityPtr()->getGlobalUniqueId());
      for (; ait != hi_ait; ait++) {
        boost::shared_ptr<EntFiniteElement> ent_fe_ptr;
        ent_fe_ptr = ait->entFePtr;
        ent_fe_ptr->row_dof_view->erase((*dit)->getGlobalUniqueId());
        if (ent_fe_ptr->row_dof_view != ent_fe_ptr->col_dof_view) {
          ent_fe_ptr->col_dof_view->erase((*dit)->getGlobalUniqueId());
        }
        ent_fe_ptr->data_dofs->get<Unique_mi_tag>().erase(
            (*dit)->getGlobalUniqueId());
      }
      dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().erase(dit);
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_ents_fields(const BitRefLevel &bit,
                                       const BitRefLevel &mask, int verb) {
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

MoFEMErrorCode Core::clear_ents_fields(const Range &ents, int verb) {
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
    dit    = entsFields.get<Ent_mi_tag>().lower_bound(first);
    hi_dit = entsFields.get<Ent_mi_tag>().upper_bound(second);
    entsFields.get<Ent_mi_tag>().erase(dit, hi_dit);
  }
  for (Field_multiIndex::iterator fit = fIelds.begin(); fit != fIelds.end();
      fit++) {
    rval = moab.tag_delete_data(fit->get()->th_AppOrder, ents);
    if (rval != MB_SUCCESS && rval != MB_TAG_NOT_FOUND) {
      CHKERRG(rval);
    } else {
      rval = MB_SUCCESS;
    }
    rval = moab.tag_delete_data(fit->get()->th_FieldData, ents);
    if (rval != MB_SUCCESS && rval != MB_TAG_NOT_FOUND) {
      CHKERRG(rval);
    } else {
      rval = MB_SUCCESS;
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_ents_fields(const std::string &name,
                                       const Range &ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_dofs_fields(name, ents, verb);
  CHKERR clear_adjacencies_entities(name, ents, verb);
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first  = p_eit->first;
    EntityHandle second = p_eit->second;
    FieldEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator
        dit,
        hi_dit, last_dit;
    dit = entsFields.get<Composite_Name_And_Ent_mi_tag>().lower_bound(
        boost::make_tuple(name, first));
    hi_dit = entsFields.get<Composite_Name_And_Ent_mi_tag>().upper_bound(
        boost::make_tuple(name, second));
    entsFields.get<Composite_Name_And_Ent_mi_tag>().erase(dit, hi_dit);
  }
  const Field *field_ptr = get_field_structure(name);
  rval                   = moab.tag_delete_data(field_ptr->th_AppOrder, ents);
  if (rval != MB_SUCCESS && rval != MB_TAG_NOT_FOUND) {
    CHKERRG(rval);
  } else {
    rval = MB_SUCCESS;
  }
  rval = moab.tag_delete_data(field_ptr->th_FieldData, ents);
  if (rval != MB_SUCCESS && rval != MB_TAG_NOT_FOUND) {
    CHKERRG(rval);
  } else {
    rval = MB_SUCCESS;
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_adjacencies_entities(const BitRefLevel &bit,
                                                const BitRefLevel &mask,
                                                int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_adjacencies_entities(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_adjacencies_entities(const Range &ents, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    const EntityHandle first  = p_eit->first;
    const EntityHandle second = p_eit->second;
    FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
        Ent_mi_tag>::type::iterator ait,
        hi_ait;
    ait    = entFEAdjacencies.get<Ent_mi_tag>().lower_bound(first);
    hi_ait = entFEAdjacencies.get<Ent_mi_tag>().upper_bound(second);
    entFEAdjacencies.get<Ent_mi_tag>().erase(ait, hi_ait);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_adjacencies_entities(const std::string &name,
                                                const Range &ents, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;

  const Field *field_ptr   = get_field_structure(name);
  int field_bit_number     = field_ptr->getBitNumber();
  bool is_distributed_mesh = basicEntityDataPtr->trueIfDistributedMesh();
  ParallelComm *pcomm      = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);

  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {

    // First and last handle
    const EntityHandle first  = p_eit->first;
    const EntityHandle second = p_eit->second;

    // Get owner proc and owner handle
    int f_owner_proc;
    EntityHandle f_moab_owner_handle;
    CHKERR pcomm->get_owner_handle(first, f_owner_proc, f_moab_owner_handle);
    int s_owner_proc;
    EntityHandle s_moab_owner_handle;
    CHKERR pcomm->get_owner_handle(second, s_owner_proc, s_moab_owner_handle);

    // Get UId
    UId first_uid = FieldEntity::getGlobalUniqueIdCalculate(
        f_owner_proc, field_bit_number, f_moab_owner_handle,
        is_distributed_mesh);
    UId second_uid = FieldEntity::getGlobalUniqueIdCalculate(
        s_owner_proc, field_bit_number, s_moab_owner_handle,
        is_distributed_mesh);

    // Find adjacencies
    FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
        Unique_mi_tag>::type::iterator ait,
        hi_ait;
    ait    = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(first_uid);
    hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(second_uid);
    entFEAdjacencies.get<Unique_mi_tag>().erase(ait, hi_ait);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_finite_elements(const BitRefLevel &bit,
                                           const BitRefLevel &mask, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_finite_elements(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_finite_elements(const Range &ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_adjacencies_finite_elements(ents, verb);
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first  = p_eit->first;
    EntityHandle second = p_eit->second;
    EntFiniteElement_multiIndex::index<Ent_mi_tag>::type::iterator fit, hi_fit;
    fit    = entsFiniteElements.get<Ent_mi_tag>().lower_bound(first);
    hi_fit = entsFiniteElements.get<Ent_mi_tag>().upper_bound(second);
    entsFiniteElements.get<Ent_mi_tag>().erase(fit, hi_fit);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_finite_elements(const std::string &name,
                                           const Range &ents, int verb) {
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

MoFEMErrorCode Core::clear_adjacencies_finite_elements(const BitRefLevel &bit,
                                                       const BitRefLevel &mask,
                                                       int verb) {
  MoFEMFunctionBegin;
  Range ents;                                                       
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_adjacencies_finite_elements(ents,verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_adjacencies_finite_elements(const Range &ents,
                                                       int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  for (Range::const_pair_iterator p_eit = ents.pair_begin();
       p_eit != ents.pair_end(); p_eit++) {
    EntityHandle first  = p_eit->first;
    EntityHandle second = p_eit->second;
    FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
        FEEnt_mi_tag>::type::iterator ait,
        hi_ait;
    ait    = entFEAdjacencies.get<FEEnt_mi_tag>().lower_bound(first);
    hi_ait = entFEAdjacencies.get<FEEnt_mi_tag>().upper_bound(second);
    entFEAdjacencies.get<FEEnt_mi_tag>().erase(ait, hi_ait);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_adjacencies_finite_elements(const std::string &name,
                                                       const Range &ents,
                                                       int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;

  FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator
      it_fe = finiteElements.get<FiniteElement_name_mi_tag>().find(name);
  if (it_fe != finiteElements.get<FiniteElement_name_mi_tag>().end()) {

    const int fe_bit_number = it_fe->get()->getBitNumber();

    for (Range::const_pair_iterator p_eit = ents.pair_begin();
         p_eit != ents.pair_end(); p_eit++) {

      // First and last handle
      const EntityHandle first  = p_eit->first;
      const EntityHandle second = p_eit->second;

      // Get UId
      UId first_uid =
          EntFiniteElement::getGlobalUniqueIdCalculate(first, fe_bit_number);
      UId second_uid =
          EntFiniteElement::getGlobalUniqueIdCalculate(second, fe_bit_number);

      // Find and remove adjacencies
      FieldEntityEntFiniteElementAdjacencyMap_multiIndex::index<
          FE_Unique_mi_tag>::type::iterator ait,
          hi_ait;
      ait    = entFEAdjacencies.get<FE_Unique_mi_tag>().lower_bound(first_uid);
      hi_ait = entFEAdjacencies.get<FE_Unique_mi_tag>().upper_bound(second_uid);
      entFEAdjacencies.get<FE_Unique_mi_tag>().erase(ait, hi_ait);
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::remove_ents_from_finite_element(const std::string &name,
                                                     const EntityHandle meshset,
                                                     const EntityType type,
                                                     int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR moab.get_entities_by_type(meshset, type, ents, false);
  CHKERR remove_ents_from_finite_element(name, ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_finite_element(const std::string &name,
                                                     const Range &ents,
                                                     int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_finite_elements(name, ents, verb);
  const EntityHandle idm = get_finite_element_meshset(name);
  CHKERR moab.remove_entities(idm, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_finite_element(const Range &ents,
                                                     int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_finite_elements(ents, verb);
  for (FiniteElement_multiIndex::iterator fe_it = finiteElements.begin();
       fe_it != finiteElements.end(); fe_it++) {
    EntityHandle meshset = fe_it->get()->getMeshset();
    CHKERR moab.remove_entities(meshset, ents);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_finite_element_by_bit_ref(
    const BitRefLevel &bit, const BitRefLevel &mask, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;                                                       
  CHKERR BitRefManager(*this).getEntitiesByRefLevel(bit, mask, ents, verb);
  CHKERR clear_finite_elements(ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_field_by_bit_ref(const BitRefLevel &bit,
                                                       const BitRefLevel &mask,
                                                       int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clear_ents_fields(bit, mask, verb);
  Field_multiIndex::iterator f_it = fIelds.begin();
  for (; f_it != fIelds.end(); f_it++) {
    EntityHandle meshset = (*f_it)->getMeshset();
    Range ents_to_remove;
    CHKERR moab.get_entities_by_handle(meshset, ents_to_remove, false);
    Range::iterator eit = ents_to_remove.begin();
    for (; eit != ents_to_remove.end();) {
      if (moab.type_from_handle(*eit) == MBENTITYSET) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      BitRefLevel bit2;
      CHKERR moab.tag_get_data(th_RefBitLevel, &*eit, 1, &bit2);
      if ((bit2 & mask) != bit2) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      if ((bit2 & bit).none()) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      FieldEntity_multiIndex::index<
          Composite_Name_And_Ent_mi_tag>::type::iterator iit;
      iit = entsFields.get<Composite_Name_And_Ent_mi_tag>().find(
          boost::make_tuple((*f_it)->getName(), *eit));
      if (iit != entsFields.get<Composite_Name_And_Ent_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Entity in the filed found, should not be there");
      }
      ++eit;
    }
    CHKERR moab.remove_entities(meshset, ents_to_remove);
    if (verb > 0) {
      PetscPrintf(cOmm, "number of removed entities = %u from field %s\n",
                  ents_to_remove.size(), (*f_it)->getName().c_str());
      if (verb > 1) {
        int num_entities;
        CHKERR moab.get_number_entities_by_handle(meshset, num_entities);
        PetscPrintf(
            cOmm, "\tnumber of entities in database = %u and meshset = %u\n",
            entsFields.get<FieldName_mi_tag>().count((*f_it)->getNameRef()),
            num_entities);
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_field(const std::string &name,
                                            const EntityHandle meshset,
                                            const EntityType type, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Range ents;
  CHKERR moab.get_entities_by_type(meshset, type, ents);
  CHKERR remove_ents_from_field(name, ents, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::remove_ents_from_field(const std::string &name,
                                            const Range &ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  EntityHandle meshset;
  meshset = get_field_meshset(name);
  CHKERR moab.remove_entities(meshset, ents);
  CHKERR clear_ents_fields(name, ents, verb);
  MoFEMFunctionReturn(0);
}

  

  MoFEMErrorCode Core::remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
    MoFEMFunctionBeginHot;
    if(verb==-1) verb = verbose;
    ierr = delete_finite_elements_by_bit_ref(bit,mask,verb); CHKERRG(ierr);
    ierr = remove_ents_from_field_by_bit_ref(bit,mask,verb); CHKERRG(ierr);
    RefEntity_multiIndex::iterator ent_it = refinedEntities.begin();
    for(;ent_it!=refinedEntities.end();) {
      BitRefLevel bit2 = (*ent_it)->getBitRefLevel();
      if((*ent_it)->getEntType()==MBENTITYSET) {
        ent_it++;
        continue;
      }
      if((bit2&mask)!=bit2) {
        ent_it++;
        continue;
      }
      if((bit2&bit).none()) {
        ent_it++;
        continue;
      }
      ent_it = refinedEntities.erase(ent_it);
    }
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode Core::delete_ents_by_bit_ref(const BitRefLevel &bit,
                                              const BitRefLevel &mask,
                                              const bool remove_parent,
                                              int verb) {
    MoFEMFunctionBegin;
    Range ents_to_delete;
    CHKERR moab.get_entities_by_handle(0, ents_to_delete, false);
    {
      Range::iterator eit = ents_to_delete.begin();
      for (; eit != ents_to_delete.end();) {
        if (moab.type_from_handle(*eit) == MBENTITYSET) {
          eit = ents_to_delete.erase(eit);
          continue;
        }
        BitRefLevel bit2;
        CHKERR moab.tag_get_data(th_RefBitLevel, &*eit, 1, &bit2);
        if ((bit2 & mask) != bit2) {
          eit = ents_to_delete.erase(eit);
          continue;
        }
        if ((bit2 & bit).none()) {
          eit = ents_to_delete.erase(eit);
          continue;
        }
        eit++;
      }
    }
    if (remove_parent) { // remove parent
      Range::iterator eit = ents_to_delete.begin();
      for (; eit != ents_to_delete.end(); eit++) {
        RefEntity_multiIndex::index<Ent_Ent_mi_tag>::type::iterator pit, hi_pit;
        pit = refinedEntities.get<Ent_Ent_mi_tag>().lower_bound(*eit);
        hi_pit = refinedEntities.get<Ent_Ent_mi_tag>().upper_bound(*eit);
        for (; pit != hi_pit; pit++) {
          EntityHandle ent = (*pit)->getRefEnt();
          if (ents_to_delete.find(ent) != ents_to_delete.end()) {
            continue;
          }
          bool success =
              refinedEntities.modify(refinedEntities.project<0>(pit),
                                     RefEntity_change_remove_parent());
          if (!success) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                    "modification unsuccessful");
          }
        }
      }
    }
    { // remove deleted entities form cubit meshsets
      for (CubitMeshSet_multiIndex::iterator cubit_it =
               getInterface<MeshsetsManager>()->getBegin();
           cubit_it != getInterface<MeshsetsManager>()->getEnd(); cubit_it++) {
        EntityHandle cubit_meshset = cubit_it->meshset;
        rval = moab.remove_entities(cubit_meshset, ents_to_delete);
        Range meshsets;
        CHKERR moab.get_entities_by_type(cubit_meshset, MBENTITYSET, meshsets);
        for (Range::iterator mit = meshsets.begin(); mit != meshsets.end();
             mit++) {
          CHKERR moab.remove_entities(*mit, ents_to_delete);
        }
      }
    }
    CHKERR remove_ents_by_bit_ref(bit, mask, verb);
    if (verb >= VERBOSE) {
      PetscSynchronizedPrintf(cOmm, "number of deleted entities = %u\n",
                              ents_to_delete.size());
      PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
    }
    if (verb >= VERY_VERBOSE) {
      EntityHandle out_meshset;
      CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER,
                                 out_meshset);
      CHKERR moab.add_entities(out_meshset, ents_to_delete);
      CHKERR moab.write_file("debug_ents_to_delete.vtk", "VTK", "",
                             &out_meshset, 1);
      CHKERR moab.delete_entities(&out_meshset, 1);
    }
    Range meshsets;
    CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshsets, true);
    for (Range::iterator mit = meshsets.begin(); mit != meshsets.end(); mit++) {
      CHKERR moab.remove_entities(*mit, ents_to_delete);
    }
    CHKERR moab.delete_entities(ents_to_delete); CHKERRQ_MOAB(rval);
    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode Core::delete_finite_elements_by_bit_ref(
    const BitRefLevel &bit,const BitRefLevel &mask,int verb
  ) {
    MoFEMFunctionBeginHot;
    if(verb==-1) verb = verbose;
    ierr = remove_ents_from_finite_element_by_bit_ref(bit,mask,verb); CHKERRG(ierr);
    RefElement_multiIndex::iterator fe_it = refinedFiniteElements.begin();
    for(;fe_it!=refinedFiniteElements.end();) {
      BitRefLevel bit2 = fe_it->getBitRefLevel();
      if(fe_it->getEntType()==MBENTITYSET) {
        fe_it++;
        continue;
      }
      if((bit2&mask)!=bit2) {
        fe_it++;
        continue;
      }
      if((bit2&bit).none()) {
        fe_it++;
        continue;
      }
      fe_it = refinedFiniteElements.erase(fe_it);
    }
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode Core::delete_finite_element(const std::string name,int verb) {
    MoFEMFunctionBeginHot;
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    FiniteElements_by_name& fe = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator miit = fe.find(name);
    if(miit==fe.end()) {
      SETERRQ1(
        PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,
        "finite element <%s> not found",name.c_str()
      );
    }
    EntityHandle meshset = (*miit)->getMeshset();
    Range ents;
    rval = moab.get_entities_by_handle(meshset,ents,false); CHKERRQ_MOAB(rval);
    ierr = remove_ents_from_finite_element(name,ents,verb); CHKERRG(ierr);
    fe.erase(miit);
    rval = moab.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
    MoFEMFunctionReturnHot(0);
  }

}
