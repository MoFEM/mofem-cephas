/** \file GetAdjacancies.cpp
 * \brief Mylti-index containers, data structures and other low-level functions
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

// const static int debug = 1;

PetscErrorCode Core::get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
  
  MoFEMFunctionBeginHot;
  RefEntity from_ref_entiti(basicEntityDataPtr,from_entiti);
  //std::cerr << "from:\n";
  //std::cerr << from_ref_entiti << std::endl;
  rval = moab.get_adjacencies(&from_entiti,1,to_dimension,false,adj_entities); CHKERRQ_MOAB(rval);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  rval = moab.tag_get_data(th_RefBitLevel,adj_entities,&*bit_levels.begin());
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
PetscErrorCode Core::get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
  
  MoFEMFunctionBeginHot;
  RefEntity from_ref_entiti(basicEntityDataPtr,from_entiti);
  //std::cerr << "from:\n";
  //std::cerr << from_ref_entiti << std::endl;
  rval = moab.get_adjacencies(&from_entiti,1,to_dimension,false,adj_entities); CHKERRQ_MOAB(rval);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  rval = moab.tag_get_data(th_RefBitLevel,adj_entities,&*bit_levels.begin());
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
PetscErrorCode Core::get_adjacencies(
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
  ierr = get_adjacencies(bit,from_entities,num_netities,to_dimension,adj_entities,operation_type); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::get_adjacencies(
  const BitRefLevel &bit,
  const EntityHandle *from_entities,
  const int num_netities,
  const int to_dimension,
  Range &adj_entities,
  const int operation_type,
  const int verb
) const {
  
  MoFEMFunctionBeginHot;
  if(verb>0) {
    std::ostringstream ss;
    ss << "from: " << bit << std::endl << "to: " << std::endl;
    PetscPrintf(cOmm,ss.str().c_str());
  }
  rval = moab.get_adjacencies(
    from_entities,num_netities,to_dimension,false,adj_entities,operation_type
  ); CHKERRQ_MOAB(rval);
  std::vector<BitRefLevel> bit_levels(adj_entities.size());
  rval = moab.tag_get_data(th_RefBitLevel,adj_entities,&*bit_levels.begin());
  std::vector<BitRefLevel>::iterator b_it = bit_levels.begin();
  Range::iterator eit = adj_entities.begin();
  //std::cerr << "to:\n";
  for(;eit!=adj_entities.end();b_it++) {
    if(verb>0) {
      RefEntity adj_entiti(basicEntityDataPtr,*eit);
      std::ostringstream ss;
      ss << "\t" << adj_entiti << std::endl;
      PetscPrintf(cOmm,ss.str().c_str());
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
