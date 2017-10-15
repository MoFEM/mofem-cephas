/** \file Tools.cpp
 * \brief Auxilairy tools
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

namespace MoFEM {

  PetscErrorCode Tools::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    MoFEMFunctionBeginHot;
    *iface = NULL;
    if(uuid == IDD_MOFEMNodeMerger) {
      *iface = dynamic_cast<Tools*>(this);
      MoFEMFunctionReturnHot(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      MoFEMFunctionReturnHot(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Tools::writeBitLevelByType(
    const BitRefLevel& bit,
    const BitRefLevel& mask,
    const EntityType type,
    const char * 	file_name,
    const char * 	file_type,
    const char * 	options
  ) const {
    MoFEM::Interface& m_field = cOre;
    moab::Interface& moab(m_field.get_moab());
    MoFEMFunctionBeginHot;
    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET,meshset); CHKERRQ_MOAB(rval);
    Range ents;
    ierr = getEntitiesByTypeAndRefLevel(bit,mask,type,meshset); CHKERRQ(ierr);
    rval = moab.write_file(file_name,file_type,options,&meshset,1); CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Tools::getEntitiesByTypeAndRefLevel(
    const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb
  ) const {
    MoFEM::Interface& m_field = cOre;
    moab::Interface& moab(m_field.get_moab());
    MoFEMFunctionBeginHot;
    Range ents;
    ierr = getEntitiesByTypeAndRefLevel(bit,mask,type,ents,verb); CHKERRQ(ierr);
    rval = moab.add_entities(meshset,ents); CHKERRQ_MOAB(rval);
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Tools::getEntitiesByTypeAndRefLevel(
    const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb
  ) const {
    MoFEM::Interface& m_field = cOre;
    moab::Interface& moab(m_field.get_moab());
    MoFEMFunctionBeginHot;
    ierr = moab.get_entities_by_type(0,type,ents,false); CHKERRQ(ierr);
    const BitRefLevel* tag_bit;
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();tag_bit++) {
      rval = moab.tag_get_by_ptr(
        cOre.get_th_RefBitLevel(),&*eit,1,(const void **)(&tag_bit)
      ); CHKERRQ_MOAB(rval);
      if(mask.any()&&tag_bit->none()) {
        eit = ents.erase(eit);
        continue;
      }
      // Not masked
      if(((*tag_bit)&mask) != (*tag_bit)) {
        eit = ents.erase(eit);
        continue;
      }
      // Not in bit
      if(((*tag_bit)&bit).none()) {
        eit = ents.erase(eit);
        continue;
      }
      eit++;
    }
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Tools::getEntitiesByRefLevel(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) const {
    MoFEM::Interface& m_field = cOre;
    moab::Interface& moab(m_field.get_moab());
    MoFEMFunctionBeginHot;
    Range ents;
    ierr = getEntitiesByRefLevel(bit,mask,ents); CHKERRQ(ierr);
    rval = moab.add_entities(meshset,ents); CHKERRQ_MOAB(rval);
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Tools::getEntitiesByRefLevel(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents) const {
    MoFEM::Interface& m_field = cOre;
    moab::Interface& moab(m_field.get_moab());
    MoFEMFunctionBeginHot;
    Range meshset_ents;
    rval = moab.get_entities_by_type(0,MBENTITYSET,meshset_ents,false); CHKERRQ_MOAB(rval);
    rval = moab.get_entities_by_handle(0,ents,false); CHKERRQ_MOAB(rval);
    ents.merge(meshset_ents);
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();) {
      switch (moab.type_from_handle(*eit)) {
        case MBVERTEX:
        case MBEDGE:
        case MBTRI:
        case MBQUAD:
        case MBTET:
        case MBPRISM:
        break;
        case MBENTITYSET:
        break;
        default:
        eit = ents.erase(eit);
        continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(cOre.get_th_RefBitLevel(),&*eit,1,&bit2); CHKERRQ_MOAB(rval);
      if(mask.any()&&bit2.none()) {
        eit = ents.erase(eit);
        continue;
      }
      if((bit2&mask) != bit2) {
        eit = ents.erase(eit);
        continue;
      }
      if((bit2&bit).none()) {
        eit = ents.erase(eit);
        continue;
      }
      eit++;
    }
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Tools::getAdjacenciesEquality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    MoFEM::Interface& m_field = cOre;
    moab::Interface& moab(m_field.get_moab());
    MoFEMFunctionBeginHot;
    RefEntity from_ref_entiti(m_field.get_basic_entity_data_ptr(),from_entiti);
    //std::cerr << "from:\n";
    //std::cerr << from_ref_entiti << std::endl;
    rval = moab.get_adjacencies(&from_entiti,1,to_dimension,false,adj_entities); CHKERRQ_MOAB(rval);
    std::vector<BitRefLevel> bit_levels(adj_entities.size());
    rval = moab.tag_get_data(cOre.get_th_RefBitLevel(),adj_entities,&*bit_levels.begin()); CHKERRQ_MOAB(rval);
    std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
    Range::iterator eit = adj_entities.begin();
    //std::cerr << "to:\n";
    for(;eit!=adj_entities.end();b_it++) {
      //RefEntity adj_entiti(moab,*eit);
      //std::cerr << "\t" << adj_entiti << std::endl;
      if(from_ref_entiti.getBitRefLevel() != *b_it/*adj_entiti.getBitRefLevel()*/) {
        eit = adj_entities.erase(eit);
      } else {
        eit++;
      }
    }
    if(b_it!=bit_levels.end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
    }
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Tools::getAdjacenciesAny(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    MoFEM::Interface& m_field = cOre;
    moab::Interface& moab(m_field.get_moab());
    MoFEMFunctionBeginHot;
    RefEntity from_ref_entiti(m_field.get_basic_entity_data_ptr(),from_entiti);
    //std::cerr << "from:\n";
    //std::cerr << from_ref_entiti << std::endl;
    rval = moab.get_adjacencies(&from_entiti,1,to_dimension,false,adj_entities); CHKERRQ_MOAB(rval);
    std::vector<BitRefLevel> bit_levels(adj_entities.size());
    rval = moab.tag_get_data(cOre.get_th_RefBitLevel(),adj_entities,&*bit_levels.begin());
    std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
    Range::iterator eit = adj_entities.begin();
    //std::cerr << "to:\n";
    for(;eit!=adj_entities.end();b_it++) {
      // RefEntity adj_entiti(moab,*eit);
      //std::cerr << "\t" << adj_entiti << std::endl;
      if(!(from_ref_entiti.getBitRefLevel()&(*b_it)).any()/*adj_entiti.getBitRefLevel()).any()*/) {
        eit = adj_entities.erase(eit);
      } else {
        eit++;
      }
    }
    if(b_it!=bit_levels.end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
    }
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Tools::getAdjacencies(
    const Problem *problem_ptr,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type,
    const int verb
  ) const {
    MoFEMFunctionBeginHot;
    BitRefLevel bit = problem_ptr->getBitRefLevel();
    ierr = getAdjacencies(bit,from_entities,num_netities,to_dimension,adj_entities,operation_type); CHKERRQ(ierr);
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Tools::getAdjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type,
    const int verb
  ) const {
    MoFEM::Interface& m_field = cOre;
    moab::Interface& moab(m_field.get_moab());
    MoFEMFunctionBeginHot;
    if(verb>0) {
      std::ostringstream ss;
      ss << "from: " << bit << std::endl << "to: " << std::endl;
      PetscPrintf(m_field.get_comm(),ss.str().c_str());
    }
    rval = moab.get_adjacencies(
      from_entities,num_netities,to_dimension,false,adj_entities,operation_type
    ); CHKERRQ_MOAB(rval);
    std::vector<BitRefLevel> bit_levels(adj_entities.size());
    rval = moab.tag_get_data(cOre.get_th_RefBitLevel(),adj_entities,&*bit_levels.begin()); CHKERRQ_MOAB(rval);
    std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
    Range::iterator eit = adj_entities.begin();
    //std::cerr << "to:\n";
    for(;eit!=adj_entities.end();b_it++) {
      if(verb>0) {
        RefEntity adj_entiti(m_field.get_basic_entity_data_ptr(),*eit);
        std::ostringstream ss;
        ss << "\t" << adj_entiti << std::endl;
        PetscPrintf(m_field.get_comm(),ss.str().c_str());
      }
      if(!((*b_it)&bit).any() ) {
        eit = adj_entities.erase(eit);
      } else {
        eit++;
      }
    }
    if(b_it!=bit_levels.end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
    }
    MoFEMFunctionReturnHot(0);
  }

}
