/** \file DeleteCore.cpp
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

#include <version.h>
#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

namespace MoFEM {

  PetscErrorCode Core::clear_inactive_dofs(int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    DofEntity_multiIndex::iterator dit;
    dit = dofsField.begin();
    for(;dit!=dofsField.end();dit++) {
      if(!(*dit)->getActive()) {
        MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type::iterator ait,hi_ait;
        ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
        hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
        for(;ait!=hi_ait;ait++) {
          boost::shared_ptr<EntFiniteElement> ent_fe_ptr;
          ent_fe_ptr = ait->entFePtr;
          ent_fe_ptr->row_dof_view->erase((*dit)->getGlobalUniqueId());
          if(ent_fe_ptr->row_dof_view!=ent_fe_ptr->col_dof_view) {
            ent_fe_ptr->col_dof_view->erase((*dit)->getGlobalUniqueId());
          }
          // if(
          //   ent_fe_ptr->row_dof_view!=ent_fe_ptr->data_dof_view||
          //   ent_fe_ptr->col_dof_view!=ent_fe_ptr->data_dof_view
          // ) {
          //   ent_fe_ptr->data_dof_view->erase((*dit)->getGlobalUniqueId());
          // }
          ent_fe_ptr->data_dofs.get<Unique_mi_tag>().erase((*dit)->getGlobalUniqueId());
        }
        dit = dofsField.erase(dit);
        if(dit==dofsField.end()) break;
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::clear_dofs_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    {
      DofEntity_multiIndex::iterator dit;
      dit = dofsField.begin();
      for(;dit!=dofsField.end();) {
        BitRefLevel bit2 = (*dit)->getBitRefLevel();
        if((*dit)->getEntType()==MBENTITYSET) {
          dit++;
          continue;
        }
        if((bit2&mask)!=bit2) {
          dit++;
          continue;
        }
        if((bit2&bit).none()) {
          dit++;
          continue;
        }
        MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type::iterator ait,hi_ait;
        ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
        hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
        for(;ait!=hi_ait;ait++) {
          boost::shared_ptr<EntFiniteElement> ent_fe_ptr;
          ent_fe_ptr = ait->entFePtr;
          ent_fe_ptr->row_dof_view->erase((*dit)->getGlobalUniqueId());
          if(ent_fe_ptr->row_dof_view!=ent_fe_ptr->col_dof_view) {
            ent_fe_ptr->col_dof_view->erase((*dit)->getGlobalUniqueId());
          }
          // if(
          //   ent_fe_ptr->row_dof_view!=ent_fe_ptr->data_dof_view||
          //   ent_fe_ptr->col_dof_view!=ent_fe_ptr->data_dof_view
          // ) {
          //   ent_fe_ptr->data_dof_view->erase((*dit)->getGlobalUniqueId());
          // }
          ent_fe_ptr->data_dofs.get<Unique_mi_tag>().erase((*dit)->getGlobalUniqueId());
        }
        dit = dofsField.erase(dit);
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::clear_dofs_fields(const std::string &name,const Range ents,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();eit++) {
      DofEntityByNameAndEnt::iterator dit,hi_dit;
      dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
      hi_dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
      for(;dit!=hi_dit;) {
        MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type::iterator ait,hi_ait;
        ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
        hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
        for(;ait!=hi_ait;ait++) {
          boost::shared_ptr<EntFiniteElement> ent_fe_ptr;
          ent_fe_ptr = ait->entFePtr;
          ent_fe_ptr->row_dof_view->erase((*dit)->getGlobalUniqueId());
          if(ent_fe_ptr->row_dof_view!=ent_fe_ptr->col_dof_view) {
            ent_fe_ptr->col_dof_view->erase((*dit)->getGlobalUniqueId());
          }
          // if(
          //   ent_fe_ptr->row_dof_view!=ent_fe_ptr->data_dof_view||
          //   ent_fe_ptr->col_dof_view!=ent_fe_ptr->data_dof_view
          // ) {
          //   ent_fe_ptr->data_dof_view->erase((*dit)->getGlobalUniqueId());
          // }
          // ent_fe_ptr->data_dofs.get<Unique_mi_tag>().erase((*dit)->getGlobalUniqueId());
        }
        dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().erase(dit);
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::clear_ents_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = clear_dofs_fields(bit,mask,verb); CHKERRQ(ierr);
    ierr = clear_adjacencies_entities(bit,mask,verb); CHKERRQ(ierr);
    MoFEMEntity_multiIndex::iterator eit;
    eit = entsFields.begin();
    for(;eit!=entsFields.end();) {
      if((*eit)->getEntType()==MBENTITYSET) {
        eit++;
        continue;
      }
      BitRefLevel bit2 = (*eit)->getBitRefLevel();
      if((bit2&mask)!=bit2) {
        eit++;
        continue;
      }
      if((bit2&bit).none()) {
        eit++;
        continue;
      }
      EntityHandle ent = (*eit)->getEnt();
      rval = moab.tag_delete_data((*eit)->sFieldPtr->th_AppOrder,&ent,1); CHKERRQ_MOAB(rval);
      if((*eit)->tag_FieldData_size>0) {
        rval = moab.tag_delete_data((*eit)->sFieldPtr->th_FieldData,&ent,1); CHKERRQ_MOAB(rval);
      }
      eit = entsFields.erase(eit);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::clear_ents_fields(const std::string &name,const Range ents,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = clear_dofs_fields(name,ents,verb); CHKERRQ(ierr);
    ierr = clear_adjacencies_entities(name,ents,verb); CHKERRQ(ierr);
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();eit++) {
      MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
      dit = entsFields.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
      hi_dit = entsFields.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
      for(;dit!=hi_dit;) {
        dit = entsFields.get<Composite_Name_And_Ent_mi_tag>().erase(dit);
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::clear_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = clear_adjacencies_finite_elements(bit,mask,verb); CHKERRQ(ierr);
    EntFiniteElement_multiIndex::iterator fe_it = entsFiniteElements.begin();
    for(;fe_it!=entsFiniteElements.end();) {
      BitRefLevel bit2 = (*fe_it)->getBitRefLevel();
      if((*fe_it)->getEntType()==MBENTITYSET) {
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
      fe_it = entsFiniteElements.erase(fe_it);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::clear_finite_elements(const std::string &name,const Range &ents,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = clear_adjacencies_finite_elements(name,ents,verb); CHKERRQ(ierr);
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();eit++) {
      EntFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator fit,hi_fit;
      fit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
      hi_fit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
      for(;fit!=hi_fit;) {
        fit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().erase(fit);
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::clear_adjacencies_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
    PetscFunctionBegin;
    MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator ait;
    ait = entFEAdjacencies.begin();
    for(;ait!=entFEAdjacencies.end();) {
      BitRefLevel bit2 = ait->entFePtr->getBitRefLevel();
      if(ait->entFePtr->getEntType()==MBENTITYSET) {
        ait++;
        continue;
      }
      if((bit2&mask)!=bit2) {
        ait++;
        continue;
      }
      if((bit2&bit).none()) {
        ait++;
        continue;
      }
      ait = entFEAdjacencies.erase(ait);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::clear_adjacencies_finite_elements(const std::string &name,const Range &ents,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();eit++) {
      MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<FEEnt_mi_tag>::type::iterator ait,hi_ait;
      ait = entFEAdjacencies.get<FEEnt_mi_tag>().lower_bound(*eit);
      hi_ait = entFEAdjacencies.get<FEEnt_mi_tag>().upper_bound(*eit);
      for(;ait!=hi_ait;) {
        if(ait->entFePtr->getName() == name) {
          ait = entFEAdjacencies.get<FEEnt_mi_tag>().erase(ait);
        } else {
          ait++;
        }
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::clear_adjacencies_entities(const BitRefLevel &bit,const BitRefLevel &mask,int verb ) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator ait;
    ait = entFEAdjacencies.begin();
    for(;ait!=entFEAdjacencies.end();) {
      BitRefLevel bit2 = ait->mofemEntPtr->getBitRefLevel();
      if(ait->mofemEntPtr->getEntType()==MBENTITYSET) {
        ait++;
        continue;
      }
      if((bit2&mask)!=bit2) {
        ait++;
        continue;
      }
      if((bit2&bit).none()) {
        ait++;
        continue;
      }
      ait = entFEAdjacencies.erase(ait);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::clear_adjacencies_entities(const std::string &name,const Range &ents,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();eit++) {
      MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<Ent_mi_tag>::type::iterator ait,hi_ait;
      ait = entFEAdjacencies.get<Ent_mi_tag>().lower_bound(*eit);
      hi_ait = entFEAdjacencies.get<Ent_mi_tag>().upper_bound(*eit);
      for(;ait!=hi_ait;) {
        if(ait->mofemEntPtr->getName() == name) {
          ait = entFEAdjacencies.get<Ent_mi_tag>().erase(ait);
        } else {
          ait++;
        }
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::remove_ents_from_field_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = clear_ents_fields(bit,mask,verb); CHKERRQ(ierr);
    Field_multiIndex::iterator f_it = fIelds.begin();
    for(;f_it!=fIelds.end();f_it++) {
      EntityHandle meshset = (*f_it)->getMeshset();
      Range ents_to_remove;
      rval = moab.get_entities_by_handle(
        meshset,ents_to_remove,false); CHKERRQ_MOAB(rval);
      Range::iterator eit = ents_to_remove.begin();
      for(;eit!=ents_to_remove.end();) {
        if(moab.type_from_handle(*eit)==MBENTITYSET) {
          eit = ents_to_remove.erase(eit);
          continue;
        }
        BitRefLevel bit2;
        rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
        if((bit2&mask)!=bit2) {
          eit = ents_to_remove.erase(eit);
          continue;
        }
        if((bit2&bit).none()) {
          eit = ents_to_remove.erase(eit);
          continue;
        }
        MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator iit;
        iit = entsFields.get<Composite_Name_And_Ent_mi_tag>().find(boost::make_tuple((*f_it)->getName(),*eit));
        if(iit != entsFields.get<Composite_Name_And_Ent_mi_tag>().end()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        eit++;
      }
      rval = moab.remove_entities(meshset,ents_to_remove); CHKERRQ_MOAB(rval);
      if(verb>0) {
        PetscPrintf(comm,
          "number of removed entities = %u from field %s\n",
          ents_to_remove.size(),
          (*f_it)->getName().c_str()
        );
        if(verb>1) {
          int num_entities;
          rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
          PetscPrintf(comm,"\tnumber of entities in database = %u and meshset = %u\n",
          entsFields.get<BitFieldId_mi_tag>().count((*f_it)->getId()),num_entities);
        }
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::remove_ents_from_field(const std::string& name,const EntityHandle meshset,const EntityType type,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    Range ents;
    rval = moab.get_entities_by_type(meshset,type,ents); CHKERRQ_MOAB(rval);
    ierr = remove_ents_from_field(name,ents,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::remove_ents_from_field(const std::string& name,const Range &ents,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    EntityHandle meshset;
    try {
      meshset = get_field_meshset(name);
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
    rval = moab.remove_entities(meshset,ents); CHKERRQ_MOAB(rval);
    ierr = clear_ents_fields(name,ents,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::remove_ents_from_finite_element_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = clear_finite_elements(bit,mask,verb); CHKERRQ(ierr);
    FiniteElement_multiIndex::iterator fe_it = finiteElements.begin();
    for(;fe_it!=finiteElements.end();fe_it++) {
      EntityHandle meshset = (*fe_it)->getMeshset();
      Range ents_to_remove;
      rval = moab.get_entities_by_handle(
        meshset,ents_to_remove,false
      ); CHKERRQ_MOAB(rval);
      Range::iterator eit = ents_to_remove.begin();
      for(;eit!=ents_to_remove.end();) {
        if(moab.type_from_handle(*eit)==MBENTITYSET) {
          eit = ents_to_remove.erase(eit);
          continue;
        }
        BitRefLevel bit2;
        rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
        if((bit2&mask)!=bit2) {
          eit = ents_to_remove.erase(eit);
          continue;
        }
        if((bit2&bit).none()) {
          eit = ents_to_remove.erase(eit);
          continue;
        }
        EntFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator iit;
        iit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().find(
          boost::make_tuple((*fe_it)->getName(),*eit)
        );
        if(iit != entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().end()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        eit++;
      }
      rval = moab.remove_entities(meshset,ents_to_remove); CHKERRQ_MOAB(rval);
      if(verb>0) {
        PetscPrintf(comm,
          "number of removed entities = %u from finite element %s\n",
          ents_to_remove.size(),(*fe_it)->getName().c_str()
        );
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::remove_ents_from_finite_element(const std::string &name,const EntityHandle meshset,const EntityType type,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    Range ents;
    rval = moab.get_entities_by_type(meshset,type,ents,false); CHKERRQ_MOAB(rval);
    ierr = remove_ents_from_finite_element(name,ents,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::remove_ents_from_finite_element(const std::string &name,const Range &ents,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = clear_finite_elements(name,ents,verb); CHKERRQ(ierr);
    const BitFEId id = getBitFEId(name);
    const EntityHandle idm = get_finite_element_meshset(id);
    rval = moab.remove_entities(idm,ents); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = delete_finite_elements_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
    ierr = remove_ents_from_field_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
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
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::delete_ents_by_bit_ref(
    const BitRefLevel &bit,const BitRefLevel &mask,const bool remove_parent,int verb
  ) {
    PetscFunctionBegin;
    Range ents_to_delete;
    rval = moab.get_entities_by_handle(0,ents_to_delete,false); CHKERRQ_MOAB(rval);
    {
      Range::iterator eit = ents_to_delete.begin();
      for(;eit!=ents_to_delete.end();) {
        if(moab.type_from_handle(*eit)==MBENTITYSET) {
          eit = ents_to_delete.erase(eit);
          continue;
        }
        BitRefLevel bit2;
        rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
        if((bit2&mask)!=bit2) {
          eit = ents_to_delete.erase(eit);
          continue;
        }
        if((bit2&bit).none()) {
          eit = ents_to_delete.erase(eit);
          continue;
        }
        eit++;
      }
    }
    if(remove_parent) { //remove parent
      Range::iterator eit = ents_to_delete.begin();
      for(;eit != ents_to_delete.end();eit++) {
        RefEntity_multiIndex::index<Ent_Ent_mi_tag>::type::iterator pit,hi_pit;
        pit = refinedEntities.get<Ent_Ent_mi_tag>().lower_bound(*eit);
        hi_pit = refinedEntities.get<Ent_Ent_mi_tag>().upper_bound(*eit);
        for(;pit!=hi_pit;pit++) {
          EntityHandle ent = (*pit)->getRefEnt();
          if(ents_to_delete.find(ent) != ents_to_delete.end()) {
            continue;
          }
          /*if(rAnk==0) {
          EntityHandle out_meshset;
          rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,out_meshset); CHKERRQ_MOAB(rval);
          rval = moab.add_entities(out_meshset,&ent,1); CHKERRQ_MOAB(rval);
          rval = moab.add_entities(out_meshset,&*eit,1); CHKERRQ_MOAB(rval);
          rval = moab.write_file("error.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
          rval = moab.delete_entities(&out_meshset,1); CHKERRQ_MOAB(rval);
          }
          std::ostringstream ss;
          ss << "child:\n" << *pit << std::endl;
          ss << "parent:\n" << RefEntity(moab,*eit) << std::endl;
          SETERRQ1(PETSC_COMM_SELF,1,
          "entity can not be removed, it is parent for some other entity\n%s",ss.str().c_str());*/
          bool success = refinedEntities.modify(
            refinedEntities.project<0>(pit),RefEntity_change_remove_parent()
          );
          if(!success) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          }
        }
      }
    }
    { //remove deleted entities form cubit meshsets
      CubitMeshSet_multiIndex::iterator cubit_it;
      cubit_it = get_meshsets_manager_ptr()->getBegin();
      for(;cubit_it!=get_meshsets_manager_ptr()->getEnd();cubit_it++) {
        EntityHandle cubit_meshset = cubit_it->meshset;
        rval = moab.remove_entities(cubit_meshset,ents_to_delete); CHKERRQ_MOAB(rval);
        Range meshsets;
        rval = moab.get_entities_by_type(cubit_meshset,MBENTITYSET,meshsets);  CHKERRQ_MOAB(rval);
        for(Range::iterator mit = meshsets.begin();mit!=meshsets.end();mit++) {
          rval = moab.remove_entities(*mit,ents_to_delete); CHKERRQ_MOAB(rval);
        }
      }
    }
    ierr = remove_ents_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
    if(verb>0) {
      PetscSynchronizedPrintf(comm,"number of deleted entities = %u\n",ents_to_delete.size());
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }
    if(verb>2) {
      EntityHandle out_meshset;
      rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,out_meshset); CHKERRQ_MOAB(rval);
      rval = moab.add_entities(out_meshset,ents_to_delete.subset_by_type(MBTET)); CHKERRQ_MOAB(rval);
      rval = moab.write_file("debug_ents_to_delete.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
      rval = moab.delete_entities(&out_meshset,1); CHKERRQ_MOAB(rval);
    }
    Range meshsets;
    rval = moab.get_entities_by_type(0,MBENTITYSET,meshsets,true);
    for(Range::iterator mit = meshsets.begin();mit!=meshsets.end();mit++) {
      rval = moab.remove_entities(*mit,ents_to_delete);
    }
    // rval = moab.delete_entities(ents_to_delete); CHKERRQ_MOAB(rval);
    for(int dd = 3;dd>=0;dd--) {
      rval = moab.delete_entities(ents_to_delete.subset_by_dimension(dd)); CHKERRQ_MOAB(rval);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::delete_finite_elements_by_bit_ref(
    const BitRefLevel &bit,const BitRefLevel &mask,int verb
  ) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = remove_ents_from_finite_element_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
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
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::delete_finite_element(const std::string name,int verb) {
    PetscFunctionBegin;
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
    ierr = remove_ents_from_finite_element(name,ents,verb); CHKERRQ(ierr);
    fe.erase(miit);
    rval = moab.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

}
