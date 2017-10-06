/** \file UpdateMeshsetsAndRanges.cpp
 * \brief Implementation of UpdateMeshsetsAndRanges interface
 * \ingroup mofem_update_meshsets_and_ranges
 *
 * MoFEM is free software: you can redistribute it and/or modify it under
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

#include <Includes.hpp>
#include <version.h>
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
#include <FEMultiIndices.hpp>
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

#include <UpdateMeshsetsAndRanges.hpp>

namespace MoFEM {

  PetscErrorCode UpdateMeshsetsAndRanges::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMUpdateMeshsetsAndRanges) {
      *iface = dynamic_cast<UpdateMeshsetsAndRanges*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }


  PetscErrorCode UpdateMeshsetsAndRanges::updateMeshsetByEntitiesChildren(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,
    const bool recursive,int verb
  ) {
    MoFEM::Interface& m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    const RefEntity_multiIndex *ref_ents_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_ref_ents(&ref_ents_ptr); CHKERRQ(ierr);
    Range ents;
    rval = moab.get_entities_by_handle(parent,ents,recursive);
    if(rval != MB_SUCCESS) {
      std::cerr << parent << std::endl;
      std::cerr << moab.type_from_handle(parent) <<  " " << MBENTITYSET << std::endl;
    } CHKERRQ_MOAB(rval);

    typedef RefEntity_multiIndex::index<Composite_ParentEnt_And_EntType_mi_tag>::type RefEntsByComposite;
    RefEntsByComposite &ref_ents = const_cast<RefEntity_multiIndex*>(ref_ents_ptr)->get<Composite_ParentEnt_And_EntType_mi_tag>();
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();eit++) {
      if(verb>2) {
        std::ostringstream ss;
        ss << "ent " << *eit << std::endl;;
        PetscPrintf(m_field.get_comm(),ss.str().c_str());
      }
      RefEntsByComposite::iterator miit = ref_ents.lower_bound(boost::make_tuple(*eit,child_type));
      RefEntsByComposite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(*eit,child_type));
      for(;miit!=hi_miit;miit++) {
        if(verb>2) {
          std::ostringstream ss;
          ss << "any bit " << *miit << std::endl;;
          PetscPrintf(m_field.get_comm(),ss.str().c_str());
        }
        if(((*miit)->getBitRefLevel()&child_bit).any()) {
          EntityHandle ref_ent = (*miit)->getRefEnt();
          if(ref_ent == *eit) continue;
          if(ref_ent == 0) {
            SETERRQ(m_field.get_comm(),MOFEM_IMPOSIBLE_CASE,"this should not happen");
          }
          if(moab.type_from_handle(*eit)==MBENTITYSET) {
            SETERRQ(m_field.get_comm(),MOFEM_IMPOSIBLE_CASE,"this should not happen");
          }
          rval = moab.add_entities(child,&ref_ent,1); CHKERRQ_MOAB(rval);
          if(verb>1) {
            std::ostringstream ss;
            ss << "good bit " << *miit << std::endl;
            PetscPrintf(m_field.get_comm(),ss.str().c_str());
          }
        }
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode UpdateMeshsetsAndRanges::updateFieldMeshsetByEntitiesChildren(const BitRefLevel &child_bit,int verb) {
    MoFEM::Interface& m_field = cOre;
    const Field_multiIndex *fields_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_fields(&fields_ptr); CHKERRQ(ierr);
    Field_multiIndex::iterator fit = fields_ptr->begin();
    for(;fit!=fields_ptr->end();fit++) {
      EntityHandle meshset = (*fit)->getMeshset();
      ierr = updateMeshsetByEntitiesChildren(meshset,child_bit,meshset,MBTET,false,verb);  CHKERRQ(ierr);
      ierr = updateMeshsetByEntitiesChildren(meshset,child_bit,meshset,MBTRI,false,verb);  CHKERRQ(ierr);
      ierr = updateMeshsetByEntitiesChildren(meshset,child_bit,meshset,MBEDGE,false,verb);  CHKERRQ(ierr);
      ierr = updateMeshsetByEntitiesChildren(meshset,child_bit,meshset,MBVERTEX,false,verb);  CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode UpdateMeshsetsAndRanges::updateFieldMeshsetByEntitiesChildren(
    const std::string name,const BitRefLevel &child_bit,int verb
  ) {
    MoFEM::Interface& m_field = cOre;
    PetscFunctionBegin;
    try {
      EntityHandle meshset = m_field.get_field_structure(name)->getMeshset();
      ierr = updateMeshsetByEntitiesChildren(meshset,child_bit,meshset,MBTET,false,verb);  CHKERRQ(ierr);
      ierr = updateMeshsetByEntitiesChildren(meshset,child_bit,meshset,MBTRI,false,verb);  CHKERRQ(ierr);
      ierr = updateMeshsetByEntitiesChildren(meshset,child_bit,meshset,MBEDGE,false,verb);  CHKERRQ(ierr);
      ierr = updateMeshsetByEntitiesChildren(meshset,child_bit,meshset,MBVERTEX,false,verb);  CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(m_field.get_comm(),e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode UpdateMeshsetsAndRanges::updateFiniteElementMeshsetByEntitiesChildren(
    const std::string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb
  ) {
    MoFEM::Interface& m_field = cOre;
    PetscFunctionBegin;
    try {
      EntityHandle meshset = m_field.get_finite_element_meshset(name);
      ierr = updateMeshsetByEntitiesChildren(meshset,child_bit,meshset,fe_ent_type,false,verb);  CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(m_field.get_comm(),e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

}
