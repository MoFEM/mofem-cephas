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

      /// tool class with methods used more than twp times
      struct Tool {

        const BitRefLevel &bIt;                         ///< bit to set
        const RefEntity_multiIndex *refEntsPtr;          ///< access to databse
        boost::shared_ptr<BasicEntityData>& baseEntData; ///< base entity data

        /// constrictor
        Tool(
          MoFEM::Interface &m_field,
          const BitRefLevel &bit,
          const RefEntity_multiIndex *ref_ents_ptr
        ):
        bIt(bit),
        refEntsPtr(ref_ents_ptr),
        baseEntData(m_field.get_basic_entity_data_ptr()) {
        }

        /// find entities and change entity bit if in databse
        PetscErrorCode findEntsToAdd(
          EntityHandle f,EntityHandle s,
          std::vector<EntityHandle>& seed_ents_vec,
          std::vector<boost::shared_ptr<RefEntity> > *shared_ref_ents_vec_for_fe = NULL
        ) const {
          PetscFunctionBegin;
          RefEntity_multiIndex::iterator rit,hi_rit;
          // get lower bound of multi-index
          rit = refEntsPtr->lower_bound(f);
          if(rit==refEntsPtr->end()) {
            // all enties in range are added to databse
            seed_ents_vec.reserve(s-f+1);
            for(;f<=s;f++) {
              seed_ents_vec.push_back(f);
            }
          } else {
            // some entities from range are in databse
            hi_rit = refEntsPtr->upper_bound(s);
            for(;f<=s;f++) {
              if(f==rit->get()->getRefEnt()) {
                // entitity is in databse, change bit levele only
                bool success = const_cast<RefEntity_multiIndex*>(refEntsPtr)->
                modify(rit,RefEntity_change_add_bit(bIt));
                if(!success) {
                  SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
                }
                if(shared_ref_ents_vec_for_fe!=NULL) {
                  shared_ref_ents_vec_for_fe->push_back(*rit);
                }
                rit++; // move to next one
                if(rit==hi_rit) {
                  // break loop, rest of the enetities in range are not in databse
                  break;
                }
              } else {
                // this entitity added to databse
                seed_ents_vec.push_back(f);
              }
            }
            // add rest entitites to vector of entitites going to be added to databse
            for(;f<=s;f++) {
              seed_ents_vec.push_back(f);
            }
          }
          PetscFunctionReturn(0);
        }

        /// add entities to databse
        PetscErrorCode addToDatabase(
          std::vector<EntityHandle>& seed_ents_vec,
          std::vector<boost::shared_ptr<RefEntity> > *shared_ref_ents_vec_for_fe = NULL
        ) const {
          PetscFunctionBegin;
          // add entitites to databse
          boost::shared_ptr<std::vector<RefEntity> > ref_ents_vec =
          boost::make_shared<std::vector<RefEntity> >();
          ref_ents_vec->reserve(seed_ents_vec.size());
          // create ref entitity instances
          for(
            std::vector<EntityHandle>::const_iterator vit=seed_ents_vec.begin();
            vit!=seed_ents_vec.end();vit++
          ) {
            ref_ents_vec->push_back(RefEntity(baseEntData,*vit));
            RefEntity_change_add_bit(bIt).operator()(ref_ents_vec->back());
          }
          std::vector<boost::shared_ptr<RefEntity> > shared_ref_ents_vec;
          shared_ref_ents_vec.reserve(ref_ents_vec->size());
          // create aliased shared pointers to ref entitity instances
          for(
            std::vector<RefEntity>::iterator vit = ref_ents_vec->begin();
            vit!=ref_ents_vec->end();vit++
          ) {
            shared_ref_ents_vec.push_back(boost::shared_ptr<RefEntity>(ref_ents_vec,&*vit));
          }
          if(shared_ref_ents_vec_for_fe) {
            shared_ref_ents_vec_for_fe->insert(
              shared_ref_ents_vec_for_fe->end(),
              shared_ref_ents_vec.begin(),shared_ref_ents_vec.end()
            );
          }
          // add shared pointers to databse
          const_cast<RefEntity_multiIndex*>(refEntsPtr)->insert(
            shared_ref_ents_vec.begin(),shared_ref_ents_vec.end()
          );
          PetscFunctionReturn(0);
        }

      };

      for(Range::const_pair_iterator pit = ents.pair_begin();pit!=ents.pair_end();pit++) {
        // get first and last element of range
        EntityHandle f = pit->first;
        EntityHandle s = pit->second;
        std::vector<EntityHandle> seed_ents_vec; // entities seeded not in database
        std::vector<boost::shared_ptr<RefEntity> > shared_ref_ents_vec_for_fe;
        // find ents to add
        ierr = Tool(m_field,bit,ref_ents_ptr).findEntsToAdd(
          f,s,seed_ents_vec,&shared_ref_ents_vec_for_fe
        ); CHKERRQ(ierr);
        // add elements
        if(!seed_ents_vec.empty()) {
          ierr = Tool(m_field,bit,ref_ents_ptr).addToDatabase(
            seed_ents_vec,&shared_ref_ents_vec_for_fe
          ); CHKERRQ(ierr);
        }

        // create finite elements
        for(
          std::vector<boost::shared_ptr<RefEntity> >::iterator vit = shared_ref_ents_vec_for_fe.begin();
          vit != shared_ref_ents_vec_for_fe.end();vit++
        ) {
          std::pair<RefElement_multiIndex::iterator,bool> p_fe;
          switch((*vit)->getEntType()) {
            case MBVERTEX:
            p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_VERTEX(*vit)))
            );
            seeded_ents.insert((*vit)->getRefEnt());
            break;
            case MBEDGE:
            p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_EDGE(*vit)))
            );
            seeded_ents.insert((*vit)->getRefEnt());
            break;
            case MBTRI:
            p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_TRI(*vit)))
            );
            seeded_ents.insert((*vit)->getRefEnt());
            break;
            case MBTET:
            p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_TET(*vit)))
            );
            seeded_ents.insert((*vit)->getRefEnt());
            break;
            case MBPRISM:
            p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_PRISM(*vit)))
            );
            if(!only_tets) {
              seeded_ents.insert((*vit)->getRefEnt());
            }
            break;
            case MBENTITYSET:
            p_fe = const_cast<RefElement_multiIndex*>(ref_fe_ptr)->insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_MESHSET(*vit)))
            );
            break;
            default:
            SETERRQ(m_field.get_comm(),MOFEM_NOT_IMPLEMENTED,"not implemented");
          }
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
          std::vector<EntityHandle> seed_ents_vec; // entities seeded not in database
          for(Range::pair_iterator pit = ents.pair_begin();pit!=ents.pair_end();pit++) {
            seed_ents_vec.clear();
            // get first and last element of range
            EntityHandle f = pit->first;
            EntityHandle s = pit->second;
            ierr = Tool(m_field,bit,ref_ents_ptr).findEntsToAdd(f,s,seed_ents_vec); CHKERRQ(ierr);
            if(!seed_ents_vec.empty()) {
              ierr = Tool(m_field,bit,ref_ents_ptr).addToDatabase(seed_ents_vec); CHKERRQ(ierr);
            }
          }
        }
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

  PetscErrorCode BitRefManager::addBitRefLevel(
    const Range &ents,const BitRefLevel &bit,int verb
  ) const {
    MoFEM::Interface &m_field = cOre;
    const RefEntity_multiIndex *ref_ent_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_ref_ents(&ref_ent_ptr);
    for(
      Range::const_pair_iterator pit = ents.const_pair_begin();
      pit!=ents.const_pair_end();pit++
    ) {
      EntityHandle first = pit->first;
      EntityHandle second = pit->second;
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator dit;
      dit = ref_ent_ptr->get<Ent_mi_tag>().lower_bound(first);
      if(dit==ref_ent_ptr->get<Ent_mi_tag>().end()) continue;
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator hi_dit;
      hi_dit = ref_ent_ptr->get<Ent_mi_tag>().upper_bound(second);
      for(;dit!=hi_dit;dit++) {
        bool success = const_cast<RefEntity_multiIndex*>(ref_ent_ptr)->modify(
          ref_ent_ptr->project<0>(dit),RefEntity_change_add_bit(bit)
        );
        if(!success) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"operation unsuccessful");
        };
        if(verb>0) {
          cerr << **dit << endl;
        }
      }
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode BitRefManager::setNthBitRefLevel(
    const Range &ents,const int n,const bool b,int verb
  ) const {
    MoFEM::Interface &m_field = cOre;
    const RefEntity_multiIndex *ref_ent_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_ref_ents(&ref_ent_ptr);
    for(
      Range::const_pair_iterator pit = ents.const_pair_begin();
      pit!=ents.const_pair_end();pit++
    ) {
      EntityHandle first = pit->first;
      EntityHandle second = pit->second;
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator dit;
      dit = ref_ent_ptr->get<Ent_mi_tag>().lower_bound(first);
      if(dit==ref_ent_ptr->get<Ent_mi_tag>().end()) continue;
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator hi_dit;
      hi_dit = ref_ent_ptr->get<Ent_mi_tag>().upper_bound(second);
      for(;dit!=hi_dit;dit++) {
        bool success = const_cast<RefEntity_multiIndex*>(ref_ent_ptr)->modify(
          ref_ent_ptr->project<0>(dit),RefEntity_change_set_nth_bit(n,b)
        );
        if(!success) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"operation unsuccessful");
        };
        if(verb>0) {
          cerr << **dit << endl;
        }
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode BitRefManager::setNthBitRefLevel(
    const int n,const bool b,int verb
  ) const {
    MoFEM::Interface &m_field = cOre;
    const RefEntity_multiIndex *ref_ent_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_ref_ents(&ref_ent_ptr);
    RefEntity_multiIndex::iterator dit,hi_dit;
    dit = ref_ent_ptr->begin();
    hi_dit = ref_ent_ptr->end();
    for(;dit!=hi_dit;dit++) {
      bool success = const_cast<RefEntity_multiIndex*>(ref_ent_ptr)->modify(
        dit,RefEntity_change_set_nth_bit(n,b)
      );
      if(!success) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"operation unsuccessful");
      };
      if(verb>0) {
        cerr << **dit << endl;
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode BitRefManager::shiftLeftBitRef(const int shift,const BitRefLevel mask,int verb) const {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    PetscFunctionReturn(0);
  }

  PetscErrorCode BitRefManager::shiftRightBitRef(const int shift,const BitRefLevel mask,int verb) const {
    MoFEM::Interface &m_field = cOre;
    const RefEntity_multiIndex *ref_ent_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_ref_ents(&ref_ent_ptr);
    for(int ii = 0;ii<shift;ii++) {
      // delete bits on the right which are shifted to zero
      BitRefLevel delete_bits = BitRefLevel().set(ii)&mask;
      if(delete_bits.none()) continue;
      ierr = m_field.delete_ents_by_bit_ref(delete_bits,delete_bits,6); CHKERRQ(ierr);
    }
    RefEntity_multiIndex::iterator ent_it = ref_ent_ptr->begin();
    for(;ent_it!=ref_ent_ptr->end();ent_it++) {
      if(verb>5) {
        std::cerr << (*ent_it)->getBitRefLevel() << " : ";
      }
      bool success = const_cast<RefEntity_multiIndex*>(ref_ent_ptr)->modify(
        ent_it,RefEntity_change_right_shift(shift,mask)
      );
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"inconsistency in data");
      if(verb>5) {
        std::cerr << (*ent_it)->getBitRefLevel() << std::endl;
      }
    }
    PetscFunctionReturn(0);

  }



}
