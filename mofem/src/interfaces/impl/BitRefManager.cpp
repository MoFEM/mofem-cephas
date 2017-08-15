/** \file BitRefManager.cpp
 * \brief Managing BitRefLevels
 * \mofem_bit_ref
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
#include <BitRefManager.hpp>

namespace MoFEM {

  PetscErrorCode BitRefManager::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMBitRefManager) {
      *iface = dynamic_cast<BitRefManager*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }

  BitRefManager::BitRefManager(const MoFEM::Core& core):
  cOre(const_cast<MoFEM::Core&>(core)),
  dEbug(false) {
  }
  BitRefManager::~BitRefManager() {
  }

  PetscErrorCode BitRefManager::setBitRefLevel(
    const Range &ents,const BitRefLevel &bit,const bool only_tets,int verb
  ) const {
    MoFEM::Interface &m_field = cOre;
    const RefEntity_multiIndex *ref_ents_ptr;
    const RefElement_multiIndex *ref_fe_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_ref_ents(&ref_ents_ptr); CHKERRQ(ierr);
    ierr = m_field.get_ref_finite_elements(&ref_fe_ptr); CHKERRQ(ierr);
    Range seeded_ents;
    try {
      if(verb > 1) {
        PetscSynchronizedPrintf(
          m_field.get_comm(),
          "nb. entities for seed %d\n",ents.size()
        );
      }
      Range::iterator tit = ents.begin();
      for(;tit!=ents.end();tit++) {
        boost::shared_ptr<RefEntity> ref_ent(
          new RefEntity(m_field.get_basic_entity_data_ptr(),*tit)
        );
        std::bitset<8> ent_pstat(ref_ent->getPStatus());
        ent_pstat.flip(0);
        std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
        const_cast<RefEntity_multiIndex*>(ref_ents_ptr)->insert(ref_ent);
        if(dEbug > 0) {
          ierr = test_moab(m_field.get_moab(),*tit); CHKERRQ(ierr);
        }
        bool success = const_cast<RefEntity_multiIndex*>(ref_ents_ptr)->
        modify(p_ent.first,RefEntity_change_add_bit(bit));
        if(!success) {
          SETERRQ(
            m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful"
          );
        }
        if(verb>2) {
          std::ostringstream ss;
          ss << **p_ent.first;
          PetscSynchronizedPrintf(m_field.get_comm(),"%s\n",ss.str().c_str());
        }
        std::pair<RefElement_multiIndex::iterator,bool> p_fe;
        switch((*p_ent.first)->getEntType()) {
          case MBVERTEX:
          p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
            boost::shared_ptr<RefElement>(new RefElement_VERTEX(*p_ent.first)))
          );
          seeded_ents.insert(*tit);
          break;
          case MBEDGE:
          p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
            boost::shared_ptr<RefElement>(new RefElement_EDGE(*p_ent.first)))
          );
          seeded_ents.insert(*tit);
          break;
          case MBTRI:
          p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
            boost::shared_ptr<RefElement>(new RefElement_TRI(*p_ent.first)))
          );
          seeded_ents.insert(*tit);
          break;
          case MBTET:
          p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
            boost::shared_ptr<RefElement>(new RefElement_TET(*p_ent.first)))
          );
          seeded_ents.insert(*tit);
          break;
          case MBPRISM:
          p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
            boost::shared_ptr<RefElement>(new RefElement_PRISM(*p_ent.first)))
          );
          if(!only_tets) {
            seeded_ents.insert(*tit);
          }
          break;
          case MBENTITYSET:
          p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
            boost::shared_ptr<RefElement>(new RefElement_MESHSET(*p_ent.first)))
          );
          break;
          default:
          SETERRQ(m_field.get_comm(),MOFEM_NOT_IMPLEMENTED,"not implemented");
        }
        if(verb>3) {
          std::ostringstream ss;
          ss << *(p_fe.first->getRefElement());
          PetscSynchronizedPrintf(m_field.get_comm(),"%s\n",ss.str().c_str());
        }
      }
      if(!seeded_ents.empty()) {
        int dim = m_field.get_moab().dimension_from_handle(seeded_ents[0]);
        for(int dd = 0;dd<dim;dd++) {
          Range ents;
          rval = m_field.get_moab().get_adjacencies(seeded_ents,dd,true,ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
          if(dd == 2 && only_tets) {
            // currently only works with triangles
            ents = ents.subset_by_type(MBTRI);
          }
          Range::iterator eit = ents.begin();
          for(;eit!=ents.end();eit++) {
            std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
            const_cast<RefEntity_multiIndex*>(ref_ents_ptr)->insert(
              boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),*eit))
            );
            bool success = const_cast<RefEntity_multiIndex*>(ref_ents_ptr)->modify(p_ent.first,RefEntity_change_add_bit(bit));
            if(!success) {
              SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
            }
            if(verb>2) {
              std::ostringstream ss;
              ss << *(*p_ent.first);
              PetscSynchronizedPrintf(m_field.get_comm(),"%s\n",ss.str().c_str());
            }
          }
        }
      }
      if(verb>2) {
        PetscSynchronizedPrintf(m_field.get_comm(),"\n");
        PetscSynchronizedFlush(m_field.get_comm(),PETSC_STDOUT);
      }
    } catch (MoFEMException const &e) {
      SETERRQ(m_field.get_comm(),e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode BitRefManager::setBitLevelToMeshset(
    const EntityHandle meshset,const BitRefLevel &bit,int verb
  ) const {
    MoFEM::Interface &m_field = cOre;
    const RefEntity_multiIndex *ref_ents_ptr;
    const RefElement_multiIndex *ref_fe_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_ref_ents(&ref_ents_ptr); CHKERRQ(ierr);
    ierr = m_field.get_ref_finite_elements(&ref_fe_ptr); CHKERRQ(ierr);
    std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
    const_cast<RefEntity_multiIndex*>(ref_ents_ptr)->insert(
      boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),meshset))
    );
    const_cast<RefEntity_multiIndex*>(ref_ents_ptr)->modify(p_ent.first,RefEntity_change_add_bit(bit));
    ptrWrapperRefElement pack_fe(
      boost::shared_ptr<RefElement>(new RefElement_MESHSET(*p_ent.first))
    );
    std::pair<RefElement_multiIndex::iterator,bool> p_fe =
    const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(pack_fe);
    if(verb > 0) {
      std::ostringstream ss;
      ss << "add meshset as ref_ent " << *(p_fe.first->getRefElement()) << std::endl;
      PetscPrintf(m_field.get_comm(),ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode BitRefManager::setBitRefLevelByDim(
    const EntityHandle meshset,const int dim,const BitRefLevel &bit,int verb
  ) const {
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    Range ents;
    rval = m_field.get_moab().get_entities_by_dimension(
      meshset,dim,ents,false
    ); CHKERRQ_MOAB(rval);
    ierr = setBitRefLevel(ents,bit,false,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode BitRefManager::setBitRefLevelByType(
    const EntityHandle meshset,const EntityType type,const BitRefLevel &bit,int verb
  ) const {
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    Range ents;
    rval = m_field.get_moab().get_entities_by_type(
      meshset,type,ents,false
    ); CHKERRQ_MOAB(rval);
    ierr = setBitRefLevel(ents,bit,false,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

}
